include("Include.jl")

# load the data dictionary -
#data_dictionary = maximize_acetate_data_dictionary(0,0,0)
#data_dictionary = maximize_atp_data_dictionary(0,0,0)
data_dictionary = maximize_cellmass_data_dictionary(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
