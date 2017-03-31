include("Include.jl")
include("AtomCheck.jl")

# load the data dictionary -
data_dictionary = maximize_acetate_data_dictionary(0,0,0)
#data_dictionary = maximize_formate_data_dictionary(0,0,0)
#data_dictionary = maximize_ethanol_data_dictionary(0,0,0)
#data_dictionary = maximize_atp_data_dictionary(0,0,0)

data_dictionary = maximize_cellmass_data_dictionary(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

#return uptake_array[73:92]

Flux = show_flux_profile(flux_array,0.5,data_dictionary)

#Check if matrix is balanced
path_to_atom_file = "./Atom.txt"

A = generate_atom_matrix("./Atom.txt",data_dictionary)

Atoms = transpose(A)*uptake_array

#return Flux
#return Atoms
return checkAllBalances(path_to_atom_file,data_dictionary)
