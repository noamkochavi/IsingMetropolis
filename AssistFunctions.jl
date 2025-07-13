"""
    AssistFunctions
    This module provides utility functions for the Ising model simulation.
"""
module AssistFunctions

export ind_mod, moving_average, idx_flat_to_2d, idx_2d_to_flat, get_neighbors, round_idxs

"""
    ind_mod(n, m)
    Modulo operation for 1-based indexing.
"""
ind_mod(n, m) = (n+m-1)%m+1

"""
    moving_average(vs, n)
    Moving average of a vector `vs` with window size `n`.
"""
moving_average(vs, n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

"""
    idx_flat_to_2d(idx, L)
    Convert a flat index to 2D index in a lattice of size LxL.
"""
function idx_flat_to_2d(idx, L)
    return ((idx - 1) ÷ L + 1, ind_mod(idx, L))
end

"""
    idx_2d_to_flat((i, j), L)
    Convert a 2D index (i, j) to a flat index in a lattice of size LxL.
"""
function idx_2d_to_flat((i, j), L)
    return L*(i-1)+j
end

"""
    get_neighbors(i, L)
    Get the indices of the nearest neighbors of a site `i` in a 2D square lattice of size LxL with periodic boundary conditions.
    Optimized for 2D square lattice only!
"""
function get_neighbors(i, L)
    a = (i-1)÷L
    b = ind_mod(i, L)
    return [ind_mod(b+1,L)+a*L, ind_mod(b-1,L)+a*L, ((a-1+L)%L)*L+b, ((a+1+L)%L)*L+b]
end

"""
    show_lattice(lattice)
    Display the lattice as a grayscale image.
    The lattice is expected to be a 1D array of values in {-1, 1}.
"""
function show_lattice(lattice, L)
    bw_array = Gray.((lattice .+ 1) ./ 2)
    img = reshape(bw_array, L, L)
    display(img)
end

"""
    print_lattice(lattice)
    Print the lattice in a human-readable format.
    The lattice is expected to be a 1D array of values in {-1, 1}.
"""
function print_lattice(lattice)
    for i in 1:L
        println(join(map(x -> (x==1) ? "\e[31m↑ \e[0m" : "\e[34m↓ \e[0m", lattice[(i-1)*L+1:i*L] .== 1), ' '))
    end
end

end