using JLD
using DSP
using MHDSim
using ProgressMeter


function generate_forcing(N, nϕ, nψ, p_switch, ϕ_min, ϕ_max, ψ_min, ψ_max)

    ϕ = zeros(nϕ, N)
    ψ = zeros(nψ, N)

    for i = 2:N
        if rand() <= p_switch
            ϕ[:, i] = ϕ[:, i-1]
            ψ[:, i] = ψ[:, i-1]
        else
            ϕ[:, i] = (ϕ_max - ϕ_min) .* rand(nϕ, 1) .+ ϕ_min
            ψ[:, i] = (ψ_max - ψ_min) .* rand(nψ, 1) .+ ψ_min
        end

    end


    responsetype = Lowpass(0.3)
    designmethod = Butterworth(1)


    for i = 1:nϕ
        ϕ[i, :] = filtfilt(digitalfilter(responsetype, designmethod), ϕ[i, :])
    end

    for i = 1:nψ
        ψ[i, :] = filtfilt(digitalfilter(responsetype, designmethod), ψ[i, :])
    end

    ϕ, ψ
end

function main()

    if length(ARGS) < 2
        println("Usage: julia generate-trajectories.jl <count> <length>")
        return
    end

    K = parse(Int, ARGS[1])
    T = parse(Float64, ARGS[2])

    Δt = 0.5
    τ = 80.0 # Desired switch time (s)
    n_el = 4 # Number of electrodes
    n_mag = 4 # Number of electromagnets
    σ = 5.0 # Conductivity (S/m)
    ν = 1.0016e-6 # Kinematic viscosity (m²/s)
    ρ = 998.2 # Density (kg/m³)
    Δt = 0.5 # Time step size (s)
    sim_step_per_ctrl = 2 # Number of simulation steps per control step

    nx = 32 # Number of grid points in x-directionj
    ny = 32 # Number of grid points in y-direction


    ϕ_min = 0.0 # Minimum voltage applied to electrodes
    ϕ_max = 10.0 # Maximum voltage applied to electrodes
    ψ_min = -1.0 # Minimum scaling factor of electromagnets
    ψ_max = 1.0 # Maximum scaling factor of electromagnets

    N = Int(T ÷ Δt) # Total number of control steps

    ## Setup Simulation ∼ Takes a while 

    function resolve_fields_path(name::String; base_dir="~/.julia/MHDSimFields", url_base::String)
        base_dir = expanduser(base_dir)
        path = joinpath(base_dir, name)
        if !isdir(path)
            println("Downloading field data '$name'...")
            mkpath(base_dir)
            archive = joinpath(base_dir, "$name.tar.gz")
            url = joinpath(url_base, "$name.tar.gz")
            Downloads.download(url, archive)
            run(`tar -xzf $archive -C $base_dir`)
        end
        return path
    end

    # Download if needed and load into simulation
    path_to_fields = resolve_fields_path(
        "fields_8mm_4el_4mag";
        url_base="https://github.com/aa4cc/MHDSim.jl/releases/download/v1.0"
    )


    # Initialize the simulation data

    data = SimData(path_to_fields, n_el, n_mag, σ, ν, ρ, Δt / sim_step_per_ctrl)

    for k = 1:K

        cache = nothing

        ϕ = zeros(4, N) # Voltages applied to electrodes
        ψ = zeros(4, N) # Scaling factors of the individual coils


        ϕ, ψ = generate_forcing(N, n_el, n_mag, 1 - Δt / τ, ϕ_min, ϕ_max, ψ_min, ψ_max)


        state = SimState(data)

        xs = LinRange(-50e-3, 50e-3, nx)
        ys = LinRange(-50e-3, 50e-3, ny)
        zs = [8e-3]

        peh = PointEvalHandler(data, xs, ys, zs)

        u, v, _ = measure_velocity(state, data, peh)

        X = zeros(2, nx, ny, N + 1)


        @showprogress desc = "Trajectory $k/$K" for t = 1:N

            for _ = 1:sim_step_per_ctrl
                # Simulate the system for sim_step_per_ctrl steps
                state, cache = solve_step(data, state, ϕ[:, t] ./ 10, ψ[:, t], cache)
            end

            u, v, _ = measure_velocity(state, data, peh)

            X[1, :, :, t+1] = reshape(u, nx, ny)
            X[2, :, :, t+1] = reshape(v, nx, ny)

        end

        JLD.save("trajectories/traj-$k.jld", "X", X, "psi", ψ, "phi", ϕ)

    end

end

main()