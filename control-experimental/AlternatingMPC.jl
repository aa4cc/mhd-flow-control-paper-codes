using OSQP
using LinearAlgebra
using ToeplitzMatrices
using SparseArrays
using Statistics


struct AlternatingMPCProblem
    H::SparseMatrixCSC{Float64, Int64}
    F::SparseMatrixCSC{Float64, Int64}
    G::SparseMatrixCSC{Float64, Int64}
    R̂_ϕ::SparseMatrixCSC{Float64, Int64}
    R̂_ψ::SparseMatrixCSC{Float64, Int64}
    R̄̄_ϕ::SparseMatrixCSC{Float64, Int64}
    R̄̄_ψ::SparseMatrixCSC{Float64, Int64}
    ϕ_min::Vector{Float64}
    ϕ_max::Vector{Float64}
    ψ_min::Vector{Float64}
    ψ_max::Vector{Float64}
    mϕ::Int
    mψ::Int
    N::Int
end



function AlternatingMPCProblem(A::SparseMatrixCSC, B::SparseMatrixCSC, C::SparseMatrixCSC, Q::SparseMatrixCSC, R_ϕ::SparseMatrixCSC, R_ψ::SparseMatrixCSC, ϕ_min, ϕ_max, ψ_min, ψ_max, N)

    n = size(A, 1)
    p = size(C, 1)

    m = size(B, 2)

    mϕ = size(R_ϕ, 1)
    mψ = size(R_ψ, 1)


    Ā = spzeros(Float64, n * N, m * N)
    Ā̄ = spzeros(Float64, n * N, n)

    C̄ = spzeros(Float64, p * N, n * N)

    C̄ = kron(sparse(I, N, N), C)
    Q̄ = kron(sparse(I, N, N), Q)

    R̄_ϕ = kron(sparse(I, N-1, N-1), R_ϕ)
    R̄_ψ = kron(sparse(I, N-1, N-1), R_ψ)

    B̄ = kron(sparse(I, N, N), B)

    tmp = Toeplitz(collect(1:N), vcat(1, zeros(N-1)))

    A_power = copy(A^0)
    
    for i = 1:N
        Ā .+= kron(tmp .== i, A_power*B) 
        A_power = A * A_power # A^i
        Ā̄ .+= kron(sparse(collect(1:N)) .== i, A_power)
    end


    #H₁ = B̄' * Ā' * C̄' * Q̄ * C̄ * Ā * B̄ # Needs to be multipled by B from the both sides
    H =  Ā' * C̄' * Q̄ * C̄ * Ā 

    #G = Ā̄' * C̄' * Q̄ * C̄ * Ā  * B̄ # Needs to be multipled by B from the right

    G = Ā̄' * C̄' * Q̄ * C̄ * Ā

    #F = Q̄ * C̄ * Ā * B̄ # r * Q * y, Needs to be multipled by B from the right

    F = Q̄ * C̄ * Ā

    R̂_ϕ = spzeros(Float64, mϕ, N*mϕ)
    R̂_ϕ[1:mϕ, 1:mϕ] .= R_ϕ

    R̂_ψ = spzeros(Float64, mψ, N*mψ)
    R̂_ψ[1:mψ, 1:mψ] .= R_ψ


    R̄̄_ϕ = kron(spdiagm(N, N, 1 => -ones(N-1), -1 => -ones(N-1)), R_ϕ) + kron(2*sparse(I, N, N),  R_ϕ)

    R̄̄_ϕ[end-mϕ+1:end, end-mϕ+1:end] ./= 2

    R̄̄_ψ = kron(spdiagm(N, N, 1 => -ones(N-1), -1 => -ones(N-1)), R_ψ) + kron(2*sparse(I, N, N),  R_ψ)

    R̄̄_ψ[end-mψ+1:end, end-mψ+1:end] ./= 2

    ϕ_min_big = vec(kron(ones(N,1), ϕ_min))
    ϕ_max_big = vec(kron(ones(N,1), ϕ_max))

    ψ_min_big = vec(kron(ones(N,1), ψ_min))
    ψ_max_big = vec(kron(ones(N,1), ψ_max))



    AlternatingMPCProblem(H, F, G, R̂_ϕ, R̂_ψ, R̄̄_ϕ, R̄̄_ψ, ϕ_min_big, ϕ_max_big, ψ_min_big, ψ_max_big, mϕ, mψ, N)

 
end



function solve_ϕ_step(problem::AlternatingMPCProblem, x0::Vector{Float64}, ψ::Matrix{Float64}, ϕ_prev::Vector{Float64}, ref::Matrix{Float64})

    H = problem.H

    F = problem.F

    G = problem.G

    R̂ = problem.R̂_ϕ

    R̄̄ = problem.R̄̄_ϕ

    N = problem.N
    m = problem.mϕ

    m_big = problem.mψ * problem.mϕ

    B_big = spzeros(m_big*N, m*N)

    for i = 1:N
        B_big[(i-1)*m_big + 1:i*m_big, (i-1)*m + 1:i*m] = kron(sparse(I, problem.mϕ, problem.mϕ), ψ[:, i])
    end

    H = B_big' * H * B_big + R̄̄  

    H = 0.5 * (H + H')

    F =  F * B_big
    G = G * B_big

    b = vcat(G, -R̂, -F)

    q = (vcat(x0, ϕ_prev, vec(ref))' * b)'

    A = sparse(Float64, I, m*N, m*N)

    model = OSQP.Model()

    settings = Dict(:verbose => false)

    OSQP.setup!(model; P=H, q=q, A=A, u=problem.ϕ_max, l=problem.ϕ_min, warm_starting=true, polish=true, settings...)

    results = OSQP.solve!(model)

    if results.info.status != :Solved
       error("OSQP did not solve the problem!")
    end


    ϕ = reshape(results.x, m, N)
    J = results.info.obj_val
    
    ϕ, J
end


function solve_ψ_step(problem::AlternatingMPCProblem, x0::Vector{Float64}, ϕ::Matrix{Float64}, ψ_prev::Vector{Float64}, ref::Matrix{Float64})

    H = problem.H

    F = problem.F

    G = problem.G

    R̂ = problem.R̂_ψ

    R̄̄ = problem.R̄̄_ψ

    N = problem.N
    m = problem.mψ

    m_big = problem.mψ * problem.mϕ

    B_big = spzeros(m_big*N, m*N)

    for i = 1:N
        B_big[(i-1)*m_big + 1:i*m_big, (i-1)*m + 1:i*m] = kron(ϕ[:, i], sparse(I, problem.mψ, problem.mψ))
    end

    H = B_big' * H * B_big +  R̄̄ #+ kron(sparse(I, N, N), 1e-2 * sparse(Float64, I, m, m))

    H = 0.5 * (H + H')

    F = F * B_big
    G = G * B_big

    b = vcat(G, -R̂, -F)

    q = (vcat(x0, ψ_prev, vec(ref))' * b)'

    A = sparse(Float64, I, m*N, m*N)

    model = OSQP.Model()

    settings = Dict(:verbose => false)

    OSQP.setup!(model; P=H, q=q, A=A, u=problem.ψ_max, l=problem.ψ_min, warm_starting=true, polish=true, settings...)

    results = OSQP.solve!(model)

    if results.info.status != :Solved
        error("OSQP did not solve the problem!")
    end

    ψ = reshape(results.x, m, N)
    J = results.info.obj_val

    ψ, J
end


function solve(problem::AlternatingMPCProblem, x0::Vector{Float64}, ϕ::Matrix{Float64}, ref::Matrix{Float64}, ϕ_prev::Vector{Float64}, ψ_prev::Vector{Float64}; maxiter=50, tol=1e-4)

    ψ = nothing
    J = nothing
    for iter = 1:maxiter

        ψ, _ = solve_ψ_step(problem, x0, ϕ, ψ_prev, ref)
        ϕ_next, J = solve_ϕ_step(problem, x0, ψ, ϕ_prev, ref)

        if sqrt.(mean(abs2, ϕ_next - ϕ)) < tol
            ϕ = ϕ_next
            println("Converged in $iter iterations")
            break
        end

        ϕ = ϕ_next

    end

    return ϕ, ψ, J
  
end
