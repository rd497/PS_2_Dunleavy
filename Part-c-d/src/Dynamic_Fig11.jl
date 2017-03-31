# Script to estimate the acetate and celmass in W3110
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("FluxDriver.jl")
include("Utility.jl")
using GLPK
using PyPlot
include("DataDictionary_noO2.jl")

# setup the time-scale -
time_start = 0.0
time_stop = 12.0
time_step = 0.1
time_array = collect(time_start:time_step:time_stop)
number_of_timesteps = length(time_array)

# Fire up the max cellmass -
data_dictionary = maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 10.5
K_glucose_uptake = 36.0

# initialize the problem -
number_of_external_states = 5
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.5     # 2 acetate
state_array[1,3] = 0.001   # 3 cellmass
state_array[1,4] = 0.001   # 3 formate
state_array[1,5] = 0.001   # 3 ethanol

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1

  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  glucose = state_array[time_step_index,1]
  acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,3]
  formate = state_array[time_step_index,4]
  ethanol = state_array[time_step_index,5]

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

  @show qGlc

  # setup the species bounds -
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  species_bounds_array[81,1] = -qGlc
  species_bounds_array[81,2] = -0.99*qGlc

  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

  # grab the fluxes from the flux_array -
  mu = 1.23*flux_array[24]
  qAcetate_production = flux_array[35]
  qFormate_production = flux_array[41]
  qEthanol_production = flux_array[40]

  @show mu

  # update the external state -
  state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
  state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
  state_array[time_step_index+1,3] = cellmass + mu*cellmass
  state_array[time_step_index+1,4] = formate + qFormate_production*cellmass
  state_array[time_step_index+1,5] = ethanol + qEthanol_production*cellmass

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)
end

glucose = state_array[:,1]
acetate = state_array[:,2]
cellmass = state_array[:,3]
formate = state_array[:,4]
ethanol = state_array[:,5]

clf()
subplot(3,2,1)
plot(time_array,cellmass)
#xlabel("Time [hr]")
ylabel("X [1/hr]")

subplot(3,2,2)
plot(time_array,glucose)
#xlabel("Time [hr]")
ylabel("Glucose [mM]")

subplot(3,2,3)
plot(time_array,acetate)
#xlabel("Time [hr]")
ylabel("Acetate [mM]")

subplot(3,2,4)
plot(time_array,formate)
xlabel("Time [hr]")
ylabel("Formate [mM]")

subplot(3,2,5)
plot(time_array,ethanol)
xlabel("Time [hr]")
ylabel("Ethanol [mM]")
