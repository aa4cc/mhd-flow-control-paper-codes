function symmetrize_trajectories(traj_path, aug_path, traj_cnt)
    traj_ptr = traj_cnt + 1 # Start from 8 because we already have 7 trajectories as the baseline

    for i = 1:traj_cnt

        traj_file = joinpath(traj_path, "traj-$i.jld")
        file_dict = JLD.load(traj_file)

        X = file_dict["X"]
        _, _, _, N = size(X)

        phi = file_dict["phi"]
        psi = file_dict["psi"]

        aug_file = joinpath(aug_path, "traj-$i.jld")
        JLD.save(aug_file, "X", X, "phi", phi, "psi", psi)


        psi_tmp = copy(psi)
        phi_tmp = copy(phi)

        psi_tmp[1, :], psi_tmp[2, :], psi_tmp[3, :], psi_tmp[4, :] = -psi[2, :], -psi[1, :], -psi[4, :], -psi[3, :]
        phi_tmp[1, :], phi_tmp[2, :], phi_tmp[3, :], phi_tmp[4, :] = phi[1, :], phi[3, :], phi[2, :], phi[4, :]

        X_tmp = copy(X)

        X_tmp[:, 1:16, :, :], X_tmp[:, 17:32, :, :] = X[:, 32:-1:17, :, :], X[:, 16:-1:1, :, :]
        X_tmp[1, :, :, :] = -X_tmp[1, :, :, :]

        aug_file = joinpath(aug_path, "traj-$traj_ptr.jld")
        JLD.save(aug_file, "X", X_tmp, "phi", phi_tmp, "psi", psi_tmp)

        traj_ptr += 1

        for _ = 1:3

            for k = axes(X, 4)
                X[1, :, :, k], X[2, :, :, k] = -rotl90(X[2, :, :, k]), rotl90(X[1, :, :, k])
            end

            psi[1, :], psi[2, :], psi[3, :], psi[4, :] = psi[2, :], psi[4, :], psi[1, :], psi[3, :]
            phi[1, :], phi[2, :], phi[3, :], phi[4, :] = phi[3, :], phi[1, :], phi[4, :], phi[2, :]

            aug_file = joinpath(aug_path, "traj-$traj_ptr.jld")
            JLD.save(aug_file, "X", X, "phi", phi, "psi", psi)

            traj_ptr += 1

            psi_tmp = copy(psi)
            phi_tmp = copy(phi)

            psi_tmp[1, :], psi_tmp[2, :], psi_tmp[3, :], psi_tmp[4, :] = -psi[2, :], -psi[1, :], -psi[4, :], -psi[3, :]
            phi_tmp[1, :], phi_tmp[2, :], phi_tmp[3, :], phi_tmp[4, :] = phi[1, :], phi[3, :], phi[2, :], phi[4, :]

            X_tmp = copy(X)

            X_tmp[:, 1:16, :, :], X_tmp[:, 17:32, :, :] = X[:, 32:-1:17, :, :], X[:, 16:-1:1, :, :]
            X_tmp[1, :, :, :] = -X_tmp[1, :, :, :]

            aug_file = joinpath(aug_path, "traj-$traj_ptr.jld")
            JLD.save(aug_file, "X", X_tmp, "phi", phi_tmp, "psi", psi_tmp) 

            traj_ptr += 1
        end

    end

    return traj_ptr - 1 # Return the total number of trajectories created
end