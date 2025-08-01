function delay_embed_block(X, d)
    nx, N = size(X)
    H = zeros(d * nx, N - d + 1)
    for di = 0:d-1
        H[di*nx+1:(di+1)*nx, :] = X[:, (di+1):end-d+1+di]
    end
    return H
end

function delay_embed(X, U, d)

    X_embedded = delay_embed_block(X, d)
    U_embedded = delay_embed_block(U, d - 1)

    Z = vcat(X_embedded, U_embedded)
    return Z

end

function extract_state(Z, nx, d)
    X = Z[(d-1)*nx+1:d*nx, :]
    return X
end

function compute_EDMDc(traj_path, ntrajs, npod, nd, nel, nmag)

    nx = npod # Number of POD modes
    nu = nel * nmag # Number of virtual inputs

    GGt = zeros(nx * nd + nu * (nd - 1) + nu, nx * nd + nu * (nd - 1) + nu)
    XpGt = zeros(nx * nd + nu * (nd - 1), nx * nd + nu * (nd - 1) + nu)

    N = 0

    for i = 1:ntrajs
        traj_file = joinpath(traj_path, "traj-$i.jld")
        file_dict = JLD.load(traj_file)

        X = file_dict["X"]
        ϕ = file_dict["phi"]
        ψ = file_dict["psi"]

        N = size(X, 2)

        U = zeros(nu, N - 1)

        for k = 1:N-1
            U[:, k] = vec(ψ[:, k] * ϕ[:, k]')
        end

        X_lifted = delay_embed(X, U, nd)

        Xp_lifted = X_lifted[:, 2:end]
        X_lifted = X_lifted[:, 1:end-1]

        G_block = vcat(X_lifted, U[:, nd:end])
        GGt += G_block * G_block'
        XpGt += Xp_lifted * G_block'

        N += size(G_block, 2)

    end

    GGt /= N
    XpGt /= N

    GGt = Symmetric(GGt)

    K = XpGt / GGt

    A = K[:, 1:end-nu]
    B = K[:, end-nu+1:end]

    C = zeros(nx, size(A, 1))

    C[:, (nd-1)*nx+1:nd*nx] .= Matrix(I, nx, nx)

    return A, B, C

end