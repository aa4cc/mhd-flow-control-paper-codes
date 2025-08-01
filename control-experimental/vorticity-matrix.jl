function vorticity_matrix_using_least_squares(nx, ny, Δx, Δy)

    in_idxs = LinearIndices((2, nx, ny))

    out_idxs = LinearIndices((nx, ny))

    DvDx = spzeros(nx * ny, 2 * nx * ny)
    DuDy = spzeros(nx * ny, 2 * nx * ny)

    for i = 1:nx

        for j = 1:ny

            if i > 2 && i < nx - 1

                DvDx[out_idxs[i, j], in_idxs[2, i-2, j]] = -2.0 / (10Δx)
                DvDx[out_idxs[i, j], in_idxs[2, i-1, j]] = -1.0 / (10Δx)
                DvDx[out_idxs[i, j], in_idxs[2, i+1, j]] = 1.0 / (10Δx)
                DvDx[out_idxs[i, j], in_idxs[2, i+2, j]] = 2.0 / (10Δx)



            else # BCs

                if i == 1
                    DvDx[out_idxs[i, j], in_idxs[2, i, j]] = -1.0 / Δx
                    DvDx[out_idxs[i, j], in_idxs[2, i+1, j]] = 1.0 / Δx
                elseif i == nx
                    DvDx[out_idxs[i, j], in_idxs[2, i-1, j]] = -1.0 / Δx
                    DvDx[out_idxs[i, j], in_idxs[2, i, j]] = 1.0 / Δx
                else
                    DvDx[out_idxs[i, j], in_idxs[2, i+1, j]] = 1.0 / (2Δx)
                    DvDx[out_idxs[i, j], in_idxs[2, i-1, j]] = -1.0 / (2Δx)
                end

            end

            if j > 2 && j < ny - 1

                DuDy[out_idxs[i, j], in_idxs[1, i, j-2]] = -2.0 / (10Δy)
                DuDy[out_idxs[i, j], in_idxs[1, i, j-1]] = -1.0 / (10Δy)
                DuDy[out_idxs[i, j], in_idxs[1, i, j+1]] = 1.0 / (10Δy)
                DuDy[out_idxs[i, j], in_idxs[1, i, j+2]] = 2.0 / (10Δy)

            else # BCs

                if j == 1
                    DuDy[out_idxs[i, j], in_idxs[1, i, j]] = -1.0 / Δy
                    DuDy[out_idxs[i, j], in_idxs[1, i, j+1]] = 1.0 / Δy
                elseif j == ny
                    DuDy[out_idxs[i, j], in_idxs[1, i, j-1]] = -1.0 / Δy
                    DuDy[out_idxs[i, j], in_idxs[1, i, j]] = 1.0 / Δy
                else
                    DuDy[out_idxs[i, j], in_idxs[1, i, j+1]] = 1.0 / (2Δy)
                    DuDy[out_idxs[i, j], in_idxs[1, i, j-1]] = -1.0 / (2Δy)
                end

            end


        end

    end


    return DvDx - DuDy

end
