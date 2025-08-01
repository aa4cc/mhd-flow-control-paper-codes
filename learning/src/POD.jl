function compute_POD(traj_path, pod_path, ntraj, pod_dim)
    nx, ny = 32, 32 # Grid size

    XXt = zeros(2 * nx * ny, 2 * nx * ny)

    for i = 1:ntraj

        traj_file = joinpath(traj_path, "traj-$i.jld")
        file_dict = JLD.load(traj_file)
        
        X = file_dict["X"]

        X = reshape(X, 2 * nx * ny, :)

        XXt += X * X'

    end

    λs, Vs = eigen(XXt)

    σs = sqrt.(reverse(λs)) # Sort eigenvalues in descending order
    U = Vs[:, end:-1:1] # Sort eigenvectors accordingly

    Ξ = U[:, 1:pod_dim]

    # Resave trajectories in POD basis
    for i = 1:ntraj
        traj_file = joinpath(traj_path, "traj-$i.jld")
        file_dict = JLD.load(traj_file)
        
        X = file_dict["X"]
        phi = file_dict["phi"]
        psi = file_dict["psi"]

        X = reshape(X, 2 * nx * ny, :)

        X_pod = Ξ' * X

        pod_file = joinpath(pod_path, "traj-$i.jld")
        JLD.save(pod_file, "X", X_pod, "phi", phi, "psi", psi)
    end

    return Ξ
end