using JLD2
using JLD
using SparseArrays
using AlternatingMPC

include("delay_embed.jl")
include("Kalman.jl")

mutable struct Controller
    X::Array{Float64,2}
    Φ::Array{Float64,2}
    Ψ::Array{Float64,2}
    U::Array{Float64,2}
    ref::Array{Float64,2}
    ϕ_init::Array{Float64,2}
    k::Int
    nd::Int
    P::SparseMatrixCSC{Float64,Int}
    problem::AlternatingMPCProblem
    filter::KalmanFilter{Float32}
end


function reference(xs, ys, ref_vec_field, ref_mag)

    nx = length(xs)
    ny = length(ys)

    p2idx = LinearIndices((2, nx, ny)) # Mapping from (component, x, y) to linear index

    S_rows = Int64[]
    S_cols = Int64[]
    S_vals = Float64[]

    idx = 1

    for i = 1:nx
        for j = 1:ny
            if ~isnan(ref_vec_field[1, i, j]) && ~isnan(ref_vec_field[2, i, j])
                push!(ref, ref_mag * ref_vec_field[:, i, j])

                push!(S_rows, idx)
                push!(S_cols, p2idx[1, i, j])
                push!(S_vals, 1.0)

                push!(S_rows, idx + 1)
                push!(S_cols, p2idx[2, i, j])
                push!(S_vals, 1.0)
                idx += 2
            end
        end
    end

    N_ref = length(ref) # Number of points where the reference is defined

    ref = vcat(ref...) # Flatten the reference vector field

    S = sparse(S_rows, S_cols, S_vals, 2 * N_ref, 2 * nx * ny) # Selection matrix for the reference

    return S, ref, N_ref
end


function setup_mpc()

    file_dict = JLD2.load("models/koopman_model.jld")

    H = sparse(file_dict["H"])
    F = sparse(file_dict["F"])
    G = sparse(file_dict["G"])
    Ξ = sparse(file_dict["Xi"])
    nd = file_dict["nd"]


    Np = 10

    nx = size(H, 1)
    npod = size(Ξ, 2)

    ## Kalman Filter

    P_init = Matrix{Float32}(zeros(nx, nx))
    x_init = zeros(Float32, nx)
    Q_kalman = zeros(Float32, nx, nx)
    Q_kalman[(nd-1)*npod+1:nd*npod, (nd-1)*npod+1:nd*npod] .= 1f-8 * Matrix{Float32}(I, npod, npod)
    R_kalman = 5f-5 * Matrix{Float32}(I, 2 * 32 * 32, 2 * 32 * 32)

    filter = KalmanFilter(x_init, P_init, x_init, P_init, convert(Matrix{Float32}, H), convert(Matrix{Float32}, F), convert(Matrix{Float32}, G), Symmetric(R_kalman), Symmetric(Q_kalman), convert(Matrix{Float32}, Ξ), nd, npod)

    ## MPC

    xs = LinRange(-50e-3, 50e-3, 32)
    ys = xs

    ref_vec_field = JLD.load("references/two-vortices.jld")["vec_field"]
    ref_mag = 1e-2

    S, ref, N_ref = reference(xs, ys, ref_vec_field, ref_mag)

    ref = repeat(ref, 1, Np)

    Q = (1e4 / (N_ref)) * sparse(I, 2 * N_ref, 2 * N_ref)

    R_ϕ = 1e-2 * sparse(I, 4, 4)
    R_ψ = 1.0 * sparse(I, 4, 4)

    ϕ_min = 0.0 * ones(4)
    ϕ_max = 10.0 * ones(4)

    ψ_min = -1.0 * ones(4)
    ψ_max = 1.0 * ones(4)

    problem = AlternatingMPCProblem(F, G, S*H, Q, R_ϕ, R_ψ, ϕ_min, ϕ_max, ψ_min, ψ_max, Np)

    ϕ_init = zeros(4, Np)

    ϕ_init[:, 1] = (ϕ_max - ϕ_min) .* rand(4) .+ ϕ_min

    ϕ_init = repeat(ϕ_init[:, 1], 1, Np)

    Φ = zeros(4, 1000)
    Ψ = zeros(4, 1000)

    X = zeros(2 * 32 * 32, 1000)

    Φ = repeat(ϕ_init[:, 1], 1, 1000)


    Controller(X, Φ, Ψ, zeros(16, 1000), ref, ϕ_init, 1, nd, Ξ, problem, filter)

end


function update_controller!(controller::Controller, vec_field::Array{Float64,3})


    nd = controller.nd
    k = copy(controller.k)
    ref = controller.ref #.* reshape(ref_mod.(k*0.5:0.5:(k-1)*0.5+5), 1, 10)

    x = permutedims(vec_field, [1, 3, 2])
    x = x[:, :, end:-1:1]
    x .*= 100e-3 / 1024
    x = reshape(x, 2 * 32 * 32)

    controller.X[:, k] = x

    data_step!(controller.filter, convert.(Float32, x))

    if k > 1
        controller.U[:, k-1] = vec(controller.Ψ[:, k-1] * controller.Φ[:, k-1]')
    end

    if k >= nd


        z0 = convert.(Float64, controller.filter.xtt)

        ϕ_t, ψ_t, J = solve(controller.problem, z0, controller.ϕ_init, ref, controller.Φ[:, k-1], controller.Ψ[:, k-1])


        controller.Φ[:, k] = ϕ_t[:, 1]
        controller.Ψ[:, k] = ψ_t[:, 1]

        controller.ϕ_init = ϕ_t

    end

    controller.k += 1

    time_step!(controller.filter, convert.(Float32, controller.U[:, k]))

    return controller.Φ[:, k], controller.Ψ[:, k]
end



function fi_inv(fi)
    i = ((-2.5 .+ sqrt.(2.5^2 .+ 4 * 5.38 * tan.((1 / 0.889) * abs.(fi)))) ./ (2 * 5.38)) .* (-1) .* sign.(fi)
end









