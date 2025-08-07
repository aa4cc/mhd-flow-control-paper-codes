using MHDSim
using AlternatingMPC
using SparseArrays
using LinearAlgebra
using JLD

include("src/delay-embed.jl")
include("src/vorticity-matrix.jl")

function main()

    # Basic parameters
    n_el = 4 # Number of electrodes
    n_mag = 4 # Number of electromagnets
    σ = 5.0 # Conductivity (S/m)
    Δt = 0.5 # Control timestep (s)
    sim_per_ctrl = 2 # Number of simulation steps per control step
    T = 200.0 # Total simulation time (s)

    nx = 32 # Number of measurement points in the x direction
    ny = 32 # Number of measurement points in the y direction
    ref_mag = 1e-2 # Magnitude of the reference velocity field (m/s)

    Np = 10 # Prediction horizon (number of steps)
    w_el = 1e-2 # Weight for the electric field in the cost function
    w_mag = 1.0 # Weight for the magnetic field in the cost function
    w_err = 1e0 # Weight for the error in the cost function

    paraview_prefix = "paraview/vorticity"


    ## MPC setup

    # Load the Koopman model
    file_dict = JLD.load("../learning/model/koopman_model.jld")
    F = sparse(file_dict["F"])
    G = sparse(file_dict["G"])
    H = sparse(file_dict["H"])
    Ξ = file_dict["Xi"]
    nd = file_dict["nd"]


    function vorticity_profile(x, y, x0, y0)

        η = 2.5e-2

        r2 = ((x - x0)^2 + (y - y0)^2) / η^2

        if r2 >= 1.0
            return NaN
        else
            return 4 / 3 * (1 - r2)
        end

    end

    xs = collect(LinRange(-50e-3, 50e-3, nx)) # x coordinates of the points where the fields will be evaluated
    ys = collect(LinRange(-50e-3, 50e-3, ny)) # y coordinates of the points where the fields will be evaluated
    zs = [8e-3] # z coordinates of the points where the fields will be evaluated

    ref = [] # Reference vector field for the MPC

    p2idx = LinearIndices((nx, ny)) # Mapping from ( x, y) to linear index

    S_rows = Int64[]
    S_cols = Int64[]
    S_vals = Float64[]

    idx = 1

    for i = 1:nx
        for j = 1:ny
            ω = vorticity_profile(xs[i], ys[j], 0.0, 0.0)
            if ~isnan(ω)
                push!(ref, ω)

                push!(S_rows, idx)
                push!(S_cols, p2idx[i, j])
                push!(S_vals, 1.0)

                idx += 1
            end
        end
    end

    N_ref = length(ref) # Number of points where the reference is defined

    ref = vcat(ref...) # Flatten the reference vector field

    S = sparse(S_rows, S_cols, S_vals, N_ref, nx * ny) # Selection matrix for the reference

    Q = w_err / N_ref * sparse(I, N_ref, N_ref) # Cost matrix for the MPC
    R_el = w_el * sparse(I, n_el, n_el) # Cost matrix for the electrodes
    R_mag = w_mag * sparse(I, n_mag, n_mag) # Cost matrix for the electromagnets

    ϕ_min = zeros(n_el) # Minimum electrode signals
    ϕ_max = 10 * ones(n_el) # Maximum electrode signals

    ψ_min = -1.0 * ones(n_mag) # Minimum electromagnet signals
    ψ_max = 1.0 * ones(n_mag) # Maximum electromagnet signals


    E = vorticity_matrix_using_least_squares(nx, ny, xs[2] - xs[1], ys[2] - ys[1]) # Matrix to estimate the vorticity from the velocity field

    problem = AlternatingMPCProblem(F, G, S * E * H, Q, R_el, R_mag, ϕ_min, ϕ_max, ψ_min, ψ_max, Np) # May take a while to initialize


    ## Simulation setup
    # Function to resolve the path to the electric and magnetic fields data.

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

    data = SimData(path_to_fields, n_el, n_mag, σ, Δt / sim_per_ctrl) # Takes a lot of time 

    state = SimState(data) # Initial state of the simulation

    peh = MHDSim.PointEvalHandler(data, xs, ys, zs) # Point evaluation handler for measuring velocity and pressure


    ## Simulation loop

    pvd = paraview_collection(paraview_prefix * ".pvd")

    vtk_grid(paraview_prefix * "-0.vtu", data.dh_v) do vtk
        vtk_point_data(vtk, data.dh_v, state.uₙ)
        vtk_save(vtk)
        pvd[0] = vtk
    end

    ts = 0.0:Δt:T
    ϕs = zeros(n_el, length(ts)) # Electrode signals
    ψs = zeros(n_mag, length(ts)) # Electromagnet signals

    vec_fields = zeros(2, nx, ny, length(ts) + 1) # Vector fields to store the velocity at each point

    ϕ_init = (ϕ_max - ϕ_min) .* rand(4) .+ ϕ_min

    ϕs[:, 1:nd-1] = repeat(ϕ_init, 1, nd - 1) # Initial electrode signals

    Φ_opt = repeat(ϕ_init, 1, Np) # Optimal electrode signals on the prediction horizon

    cache = nothing
    for (k, t) in enumerate(ts)

        if k >= nd
            Y = reshape(vec_fields[:, :, :, k-nd+1:k], :, nd)
            Ỹ = mapreduce(x -> Ξ' * x, hcat, eachcol(Y)) # Last nd of POD reduced measurements

            U = mapreduce(x -> vec(x[1]x[2]'), hcat, zip(eachcol(ψs[:, k-nd+1:k-1]), eachcol(ϕs[:, k-nd+1:k-1]))) # Last nd-1 of virtual control inputs

            z₀ = vec(delay_embed(Ỹ, U, nd))

            Φ_opt, Ψ_opt, _, _ = solve(problem, z₀, Φ_opt, repeat(ref, 1, Np), ϕs[:, k-1], ψs[:, k-1])

            # Extract the first control signals from the optimal solution
            ϕs[:, k] = Φ_opt[:, 1]
            ψs[:, k] = Ψ_opt[:, 1]

        end

        # Run the simulation for two steps

        for i = 1:sim_per_ctrl
            # Solve the simulation step with the current control signals
            state, cache = solve_step(data, state, ϕs[:, k] ./ 10, ψs[:, k], cache)
        end

        vx, vy, _ = measure_velocity(state, data, peh) # Measure velocity at the points

        vec_fields[1, :, :, k+1] = reshape(vx, nx, ny) # Store x-component of the velocity
        vec_fields[2, :, :, k+1] = reshape(vy, nx, ny) # Store y-component of the velocity

        @info "Step $k at time $t: Electrode signals: $(ϕs[:, k]), Electromagnet signals: $(ψs[:, k])"

        vtk_grid(paraview_prefix * "-$t.vtu", data.dh_v) do vtk
            vtk_point_data(vtk, data.dh_v, state.uₙ)
            vtk_save(vtk)
            pvd[k] = vtk
        end


    end

    vtk_save(pvd)

end

main()