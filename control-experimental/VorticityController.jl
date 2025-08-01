using JLD2
using SparseArrays
using AlternatingMPC

include("delay_embed.jl")
include("Kalman.jl")
include("vorticity-matrix.jl")

mutable struct VorticityController
    X::Array{Float64,2}
    Φ::Array{Float64,2}
    Ψ::Array{Float64,2}
    U::Array{Float64,2}
    J::Array{Float64}
    ref::Array{Float64,2}
    ϕ_init::Array{Float64,2}
    Np::Int
    k::Int
    nd::Int
    P::SparseMatrixCSC{Float64,Int}
    problem::AlternatingMPCProblem
    filter::KalmanFilter{Float32}
end


function ref_mod(k)
    if k > 203
        return -1.0
    else
        return 1.0
    end
end

function vorticity_profile(x, y, x0, y0)

    η = 2.5e-2

    r2 = ((x - x0)^2 + (y - y0)^2) / η^2

    if r2 >= 1.0
        return NaN
    else
        return 4 / 3 * (1 - r2)
    end

end

function reference(xs, ys)

    nx = length(xs)
    ny = length(ys)

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

    return S, ref, N_ref
end


function setup_mpc()

    file_dict = JLD2.load("models/koopman_model.jld")

    F = sparse(file_dict["F"])
    G = sparse(file_dict["G"])
    H = sparse(file_dict["H"])
    Ξ = sparse(file_dict["Xi"])
    nd = file_dict["nd"]

    nx = size(F, 1)
    npod = size(Ξ, 2)


    ## Kalman Filter

    P_init = 1f-3 * Matrix{Float32}(I, nx, nx)
    x_init = zeros(Float32, nx)
    Q_kalman = zeros(Float32, nx, nx)
    Q_kalman[(nd-1)*npod+1:nd*npod, (nd-1)*npod+1:nd*npod] .= 1f-8 * Matrix{Float32}(I, npod, npod)
    R_kalman = 1f-3 * Matrix{Float32}(I, 2 * 32 * 32, 2 * 32 * 32)

    filter = KalmanFilter(x_init, P_init, x_init, P_init, convert(Matrix{Float32}, F), convert(Matrix{Float32}, G), convert(Matrix{Float32}, H), Symmetric(R_kalman), Symmetric(Q_kalman), convert(Matrix{Float32}, Ξ), nd, npod)

    ## MPC 

    Np = 10

    xs = LinRange(-50e-3, 50e-3, 32)
    ys = xs

    S, ref, N_ref = reference(xs, ys)

    ref = repeat(ref, 1, Np)

    Q = (1e0 / (N_ref)) * sparse(I, N_ref, N_ref)
    R_ϕ = 1e-2 * sparse(I, 4, 4)
    R_ψ = 1.0 * sparse(I, 4, 4)

    ϕ_min = 0.0 * ones(4)
    ϕ_max = 10.0 * ones(4)

    ψ_min = -1.0 * ones(4)
    ψ_max = 1.0 * ones(4)

    E = vorticity_matrix_using_least_squares(32, 32, xs[2] - xs[1], ys[2] - ys[1])


    problem = AlternatingMPCProblem(F, G, S * E * H, Q, R_ϕ, R_ψ, ϕ_min, ϕ_max, ψ_min, ψ_max, Np)

    ϕ_init = zeros(4, Np)

    ϕ_init[:, 1] = (ϕ_max - ϕ_min) .* rand(4) .+ ϕ_min

    ϕ_init = repeat(ϕ_init[:, 1], 1, Np)

    Φ = zeros(4, 1000)
    Ψ = zeros(4, 1000)
    J = zeros(1000)

    X = zeros(2 * 32 * 32, 1000)

    Φ = repeat(ϕ_init[:, 1], 1, 1000)

    VorticityController(X, Φ, Ψ, zeros(16, 1000), J, ref, ϕ_init, Np, 1, nd, Ξ, problem, filter)

end


function update_controller!(controller::VorticityController, vec_field::Array{Float64,3})


    nd = controller.nd
    k = copy(controller.k)
    ref = controller.ref .* reshape(ref_mod.(k:(k-1)+controller.Np), 1, controller.Np)

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

        ϕ_t, ψ_t, J = solve(controller.problem, z0, controller.ϕ_init, ref, controller.Φ[:, k-1], controller.Ψ[:, k-1]; maxiter=30)

        controller.J[k] = J

        @info "J = $J"

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









