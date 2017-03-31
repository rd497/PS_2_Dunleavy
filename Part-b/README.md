## CHEME7770 *E.coli* Core Network Example
This repository holds the flux balance analysis problem (FBA) for the *E.coli* network example from [Palsson laboratory](http://systemsbiology.ucsd.edu/Downloads/EcoliCore) and is described in [EcoSal Chapter 10.2.1 - Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide by Orth, Fleming, and Palsson (2010)](http://www.asmscience.org/content/journal/ecosalplus/10.1128/ecosalplus.10.2.1#backarticlefulltext).

The core *E.coli* network is a small-scale model that can by senior undergraduate and first-year graduate students learning about constraint-based modeling and systems biology. This model has enough reactions and pathways to enable interesting and insightful calculations, but it is also simple enough that the results of such calculations can be understood easily.

### Installation and Requirements
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/CHEME7770-Ecoli-Core-Network-Example.git

or

	$ git clone https://github.com/varnerlab/CHEME7770-Ecoli-Core-Network-Example.git
	
The constraints are encoded the [DataDictionary.jl](https://github.com/varnerlab/CHEME7770-Ecoli-Core-Network-Example/blob/master/src/DataDictionary.jl) file, a basic driver to perform the flux balance analysis calculation is given in [Solve.jl](https://github.com/varnerlab/CHEME7770-Ecoli-Core-Network-Example/blob/master/src/Solve.jl).

### Requirements
The flux balance analysis problem is solved as a [Linear Programming (LP)](https://www.math.ucla.edu/~tom/LP.pdf) problem using the [GLPK solver](https://www.gnu.org/software/glpk/). You can install the [GLPK package for Julia](https://github.com/JuliaOpt/GLPK.jl) using the ``Pkg.add`` command from the REPL:
	
	julia> Pkg.add("GLPK")



 
