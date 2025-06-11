include("Metropolis.jl")

using .Metropolis
using CSV, DataFrames

# Initialize parameters
seeds = [9413]  # Random seed(s) for reproducibility
L = 100  # Lattice size

# Temperature array
Ts = vcat([4.0, 3.0], 2.6:-0.02:2.0, 1.8:-0.2:1.0) # Temperatures to simulate
equ_sweeps = fill(10000, length(Ts))  # Array of sweeps for equilibration
meas_sweeps = fill(1000, length(Ts))  # Array of sweeps for measurement
n_meas = fill(10, length(Ts))  # Array of number of measurements to take

# Loop over seeds
for seed in seeds
    save_dir = "meas_trials\\trial_$(seed)_L_$(L)_Ts_$(length(Ts))"  # Directory to save results
    
    # Create directory if it doesn't exist
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    # Save setup parameters to a CSV file
    df = DataFrame((T=Ts, equ_sweeps=equ_sweeps, meas_sweeps=meas_sweeps, n_meas=n_meas))
    CSV.write("$(save_dir)\\setup.csv", df)
    
    # Perform measurements
    T_arr, M_arr, meas_list = meas_trial(seed, L, Ts, equ_sweeps, meas_sweeps, n_meas)
    
    # Save lattice configurations for each measurement
    for meas in meas_list
        lat = round.(reshape(meas["lattice"], L, L))
        df = DataFrame(lat, :auto)
        CSV.write("$(save_dir)\\T_$(meas["T"])_meas_$(meas["meas_num"]).csv", df; writeheader=false)
    end
    
    # Save magnetization measurements to a CSV file
    open("$(save_dir)\\M_meas.csv", "w") do io
        println(io, "T,M")  # Write header
        for i in eachindex(M_arr)
            println(io, join([T_arr[i], M_arr[i]], ','))  # Write temperature and magnetization
        end
    end
end