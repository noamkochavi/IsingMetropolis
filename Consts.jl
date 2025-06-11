"""
    Consts
    This module defines constants used in the simulation of the Ising model.
"""
module Consts
export dim, z, J, possible_s_sum

const dim = 2 # Dimension of the lattice
const z = 2^dim # Number of nearest neighbors
const J = 1 # Coupling constant
const possible_s_sum = [2:2:z;] # Possible values of the sum of spins relevant for the energy calculation

end