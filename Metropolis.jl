"""
    Metropolis
    This module implements the Metropolis algorithm for the Ising model simulation.
    It includes functions for performing Metropolis steps, sweeps, and trials.
"""
module Metropolis

export step_debug, step, sweep, fastersweep, do_sweeps, do_trial, meas_trial

include("Consts.jl")
include("AssistFunctions.jl")
using .Consts
using .AssistFunctions
using Random

"""
    step_debug(lattice, N, β, M, exp_lookup, neighbor_lookup)
    -----------------------------------
    Perform a single Metropolis step on the lattice with debug output.
    This function randomly selects a site, calculates the change in energy (ΔE) for flipping the spin at that site,
    and decides whether to flip the spin based on the Metropolis criterion.
    It prints debug information about the selected site, the change in energy, and whether the spin was flipped.
    This function modifies the lattice in place and updates the magnetization `M`.
    Parameters:
    - `lattice`: The current state of the lattice (array of spins).
    - `N`: The total number of sites in the lattice.
    - `β`: The inverse temperature (1/kT).
    - `M`: The current magnetization of the lattice.
    - `exp_lookup`: Precomputed exponential values for efficiency.
    - `neighbor_lookup`: A lookup table for the nearest neighbors of each site.
    Returns:
    - `M`: The updated magnetization after the step.
"""
function step_debug(lattice, N, β, M, exp_lookup, neighbor_lookup)
    # Get a random index in the lattice, find its neighbors, and calculate the change in energy
    r_idx = rand(1:N)
    nn = neighbor_lookup[r_idx]
    ΔE = 2J*lattice[r_idx]*sum(lattice[nn])

    # If the change in energy is less than or equal to zero, flip the spin
    # Otherwise, decide whether to flip the spin based on the Metropolis criterion
    if ΔE <= 0
        lattice[r_idx] *= (-1)
        println("replaced (-): i=$r_idx, ΔE=$ΔE")
        M += 2*lattice[r_idx]
    else
        # Generate a random number and compare it with the acceptance probability
        # If the random number is less than the acceptance probability, flip the spin
        r = rand()
        A = exp_lookup[-β*ΔE]
        if r < A
            lattice[r_idx] *= (-1)
            M += 2*lattice[r_idx]
            println("replaced (+): i=$r_idx, A=$A, r=$r, ΔE=$ΔE")
        else
            println("not replaced (+): i=$r_idx, A=$A, r=$r, ΔE=$ΔE")
        end
    end
    return M
end

"""
    step(lattice, N, β, M, exp_lookup, neighbor_lookup)
    -----------------------------------
    Perform a single Metropolis step on the lattice without debug output.
    This function randomly selects a site, calculates the change in energy (ΔE) for flipping the spin at that site,
    and decides whether to flip the spin based on the Metropolis criterion.
    This function modifies the lattice in place and updates the magnetization `M`.
    Parameters:
    - `lattice`: The current state of the lattice (array of spins).
    - `N`: The total number of sites in the lattice.
    - `β`: The inverse temperature (1/kT).
    - `M`: The current magnetization of the lattice.
    - `exp_lookup`: Precomputed exponential values for efficiency.
    - `neighbor_lookup`: A lookup table for the nearest neighbors of each site.
    Returns:
    - `M`: The updated magnetization after the step.
"""
function step(lattice, N, β, M, exp_lookup, neighbor_lookup)
    # Get a random index in the lattice, find its neighbors, and calculate the change in energy
    r_idx = rand(1:N)
    nn = neighbor_lookup[r_idx]
    ΔE = 2J*lattice[r_idx]*sum(lattice[nn])

    # If the change in energy is less than or equal to zero, flip the spin
    # Otherwise, decide whether to flip the spin based on the Metropolis criterion
    if ΔE <= 0
        lattice[r_idx] *= (-1)
        M += 2*lattice[r_idx]
    else
        # Generate a random number and compare it with the acceptance probability
        # If the random number is less than the acceptance probability, flip the spin
        r = rand()
        A = exp_lookup[-β*ΔE]
        if r < A
            lattice[r_idx] *= (-1)
            M += 2*lattice[r_idx]
        end
    end
    return M
end

"""
    sweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
    -----------------------------------
    Perform a full sweep of Metropolis steps on the lattice.
    This function iterates over a fixed number of Metropolis steps (equal to the number of sites in the lattice),
    calling the `step` function for each step to update the lattice and magnetization.
    Parameters:
    - `lattice`: The current state of the lattice (array of spins).
    - `N`: The total number of sites in the lattice.
    - `β`: The inverse temperature (1/kT).
    - `M`: The current magnetization of the lattice.
    - `exp_lookup`: Precomputed exponential values for efficiency.
    - `neighbor_lookup`: A lookup table for the nearest neighbors of each site.
    Returns:
    - `M`: The updated magnetization after the sweep.
"""
function sweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
    # Perform N Metropolis steps, updating the lattice and magnetization
    for i in range(1, N)
        # Call the `step` function to perform a single Metropolis step
        M = step(lattice, N, β, M, exp_lookup, neighbor_lookup)
    end
    return M
end

"""
    fastersweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
    -----------------------------------
    Perform an optimized sweep of Metropolis steps on the lattice.
    This function performs a fixed number of Metropolis steps, randomly selecting sites and calculating the change in energy (ΔE) for flipping the spin at those sites.
    It modifies the lattice in place and updates the magnetization `M`.
    The steps are implemented in this function to improve performance by reducing the number of function calls and using precomputed values.
    Parameters:
    - `lattice`: The current state of the lattice (array of spins).
    - `N`: The total number of sites in the lattice.
    - `β`: The inverse temperature (1/kT).  
    - `M`: The current magnetization of the lattice.
    - `exp_lookup`: Precomputed exponential values for efficiency.
    - `neighbor_lookup`: A lookup table for the nearest neighbors of each site.
    Returns:
    - `M`: The updated magnetization after the sweep.
"""
function fastersweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
    # Perform a fixed number of Metropolis steps
    for i in range(1, N)
        # Get a random index in the lattice, find its neighbors, and calculate the change in energy
        r_idx = rand(1:N)
        nn = neighbor_lookup[r_idx]
        ΔE = 2J*lattice[r_idx]*sum(lattice[nn])
        
        # If the change in energy is less than or equal to zero, flip the spin
        # Otherwise, decide whether to flip the spin based on the Metropolis criterion
        if ΔE <= 0 || rand() < exp_lookup[-β*ΔE]
            lattice[r_idx] *= (-1)
            M += 2*lattice[r_idx]
        end
    end
    return M
end

"""
    do_sweeps(lattice, L, β, M0; n=100)
    -----------------------------------
    Perform multiple sweeps of the Metropolis algorithm on the lattice.
    This function precomputes the exponential lookup table and neighbor lookup table for efficiency,
    and then performs `n` sweeps using the `fastersweep` function.
    Parameters:
    - `lattice`: The current state of the lattice (array of spins).
    - `L`: The linear size of the lattice (LxL grid).
    - `β`: The inverse temperature (1/kT).
    - `M0`: The initial magnetization of the lattice.
    - `n`: The number of sweeps to perform (default: 100).
    Returns:
    - `M`: The updated magnetization after all sweeps.
"""
function do_sweeps(lattice, L, β, M0; n=100)
    # Precompute possible exponential values for the Metropolis criterion
    possible_exp_par = -2*β*J.*possible_s_sum
    exp_lookup = Dict(x => ℯ^x for x in possible_exp_par)

    # Compute the total number of sites in the lattice
    N = L^2

    # Precompute the neighbor lookup table for all sites
    neighbor_lookup = Dict(x => get_neighbors(x, L) for x in 1:N)
    
    # Initialize the magnetization
    M = M0

    # Perform `n` sweeps using the `fastersweep` function
    for i in range(1, n)
        M = fastersweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
    end

    return M
end

"""
    do_trial(βs, L, init_mode; n_sweeps=100)
    -----------------------------------
    Perform a trial of the Metropolis algorithm for a range of inverse temperatures (βs).
    This function initializes the lattice based on the specified `init_mode`, performs sweeps for each β,
    and returns the magnetization values after each set of sweeps.
    Parameters:
    - `βs`: A vector of inverse temperatures (1/kT) to iterate over.
    - `L`: The linear size of the lattice (LxL grid).
    - `init_mode`: The initialization mode for the lattice. Can be:
        - `:pos` for all spins up (+1),
        - `:neg` for all spins down (-1),
        - `:rand` for random spins (+1 or -1).
    - `n_sweeps`: The number of sweeps to perform for each β. Can be:
        - An integer (same number of sweeps for all βs),
        - A vector (specific number of sweeps for each β).
    Returns:
    - `M_arr`: A vector of magnetization values after each set of sweeps.
"""
function do_trial(βs, L, init_mode; n_sweeps=100)
    # Compute the total number of sites in the lattice
    N = L^2

    # Initialize the lattice based on the specified mode
    if init_mode == :pos
        lattice = ones(Int8, N) # Stable array with all spins up
    elseif init_mode == :neg
        lattice = -ones(Int8, N) # Stable array with all spins down
    elseif init_mode == :rand
        lattice = 2 .* rand(Bool, N) .- 1 # Random array of spins (+1 or -1)
    else
        throw(ArgumentError("Invalid init_mode. Use :pos, :neg, or :rand."))
    end

    # Ensure `n_sweeps` is valid
    if n_sweeps isa Int
        n_sweeps = fill(n_sweeps, length(βs)) # Use the same number of sweeps for all βs
    elseif !(n_sweeps isa Vector)
        throw(ArgumentError("n_sweeps should be an Int or a Vector."))
    elseif length(n_sweeps) != length(βs)
        throw(ArgumentError("n_sweeps should have the same length as βs."))
    end

    # Initialize the magnetization
    M = sum(lattice)

    # Array to store magnetization values after each set of sweeps
    M_arr = Float64[]

    # Perform sweeps for each β
    for (n, β) in zip(n_sweeps, βs)
        M = do_sweeps(lattice, L, β, M, n=n) # Use L for lattice size
        append!(M_arr, M) # Store the magnetization
    end

    return M_arr
end

"""
    meas_trial(seed, L, Ts, equ_sweeps, meas_sweeps, n_meas)
    --------------------------------------------------------
    Perform a measurement trial for the Ising model using the Metropolis algorithm.
    This function initializes the lattice, equilibrates it for a given number of sweeps,
    and then performs measurements at specified temperatures.

    This function uses the `fastersweep` function to perform the Metropolis steps and
    implements the loop over the sweeps for equilibration and measurements.

    Parameters:
    - `seed`: Random seed for reproducibility.
    - `L`: The linear size of the lattice (LxL grid).
    - `Ts`: A vector of temperatures to iterate over.
    - `equ_sweeps`: Number of sweeps for equilibration. Can be:
        - An integer (same number of sweeps for all temperatures),
        - A vector (specific number of sweeps for each temperature).
    - `meas_sweeps`: Number of sweeps between measurements. Can be:
        - An integer (same number of sweeps for all temperatures),
        - A vector (specific number of sweeps for each temperature).
    - `n_meas`: Number of measurements to perform. Can be:
        - An integer (same number of measurements for all temperatures),
        - A vector (specific number of measurements for each temperature).

    Returns:
    - `T_arr`: An array of temperatures corresponding to each sweep.
    - `M_arr`: An array of magnetization values corresponding to each sweep.
    - `meas_list`: A list of dictionaries containing measurement data (temperature, measurement number, and lattice state).
"""
function meas_trial(seed, L, Ts, equ_sweeps, meas_sweeps, n_meas)
    # Set the random seed for reproducibility
    Random.seed!(seed)

    # Compute the total number of sites in the lattice
    N = L^2

    # Initialize the lattice with random spins (+1 or -1)
    lattice = 2 .* rand(Bool, N) .- 1

    # Initialize the magnetization
    M = sum(lattice)

    # Precompute the neighbor lookup table for all sites
    neighbor_lookup = Dict(x => get_neighbors(x, L) for x in 1:N)

    # Ensure `equ_sweeps` is valid
    if equ_sweeps isa Int
        equ_sweeps = fill(equ_sweeps, length(Ts)) # Use the same number of sweeps for all temperatures
    elseif !(equ_sweeps isa Vector)
        throw(ArgumentError("equ_sweeps should be Int or Vector"))
    elseif length(equ_sweeps) != length(Ts)
        throw(ArgumentError("equ_sweeps should be the same length as Ts"))
    end

    # Ensure `meas_sweeps` is valid
    if meas_sweeps isa Int
        meas_sweeps = fill(meas_sweeps, length(Ts)) # Use the same number of sweeps for all temperatures
    elseif !(meas_sweeps isa Vector)
        throw(ArgumentError("meas_sweeps should be Int or Vector"))
    elseif length(meas_sweeps) != length(Ts)
        throw(ArgumentError("meas_sweeps should be the same length as Ts"))
    end

    # Ensure `n_meas` is valid
    if n_meas isa Int
        n_meas = fill(n_meas, length(Ts)) # Use the same number of measurements for all temperatures
    elseif !(n_meas isa Vector)
        throw(ArgumentError("n_meas should be Int or Vector"))
    elseif length(n_meas) != length(Ts)
        throw(ArgumentError("n_meas should be the same length as Ts"))
    end

    # Calculate the total number of sweeps for all temperatures
    n_sweeps = equ_sweeps .+ meas_sweeps .* (n_meas .- 1)

    # Initialize arrays to store magnetization and temperature values
    M_arr = fill(0.0, sum(n_sweeps)) # Magnetization array
    T_arr = fill(0.0, sum(n_sweeps)) # Temperature array

    # Initialize a list to store measurement data
    meas_list = []

    # Index for storing data in arrays
    arr_idx = 1

    # Iterate over temperatures
    for (Ti, T) in enumerate(Ts)
        β = 1 / T # Compute inverse temperature

        # Precompute possible exponential values for the Metropolis criterion
        possible_exp_par = -2 * β * J .* possible_s_sum
        exp_lookup = Dict(x => ℯ^x for x in possible_exp_par)

        # Perform equilibration sweeps
        for _ in 1:equ_sweeps[Ti]
            M = fastersweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
            M_arr[arr_idx] = M
            T_arr[arr_idx] = T
            arr_idx += 1
        end

        # Store the lattice state after equilibration
        push!(meas_list, Dict("T" => T, "meas_num" => 1, "lattice" => copy(lattice)))

        # Perform measurement sweeps
        for meas_num in 2:n_meas[Ti]
            for _ in 1:meas_sweeps[Ti]
                M = fastersweep(lattice, N, β, M, exp_lookup, neighbor_lookup)
                M_arr[arr_idx] = M
                T_arr[arr_idx] = T
                arr_idx += 1
            end
            # Store the lattice state after each measurement
            push!(meas_list, Dict("T" => T, "meas_num" => meas_num, "lattice" => copy(lattice)))
        end
    end

    return T_arr, M_arr, meas_list
end

end