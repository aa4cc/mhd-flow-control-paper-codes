using Makie
using CairoMakie
using JLD
using Statistics
using Interpolations
using LinearAlgebra
using ColorSchemes
using SparseArrays

Δt = 0.5

set_theme!(fontsize=7)

function compute_error(X, ref)

    N = size(X, 4)

    Nx = size(X, 2)
    Ny = size(X, 3)

    err = zeros(N)

    N_points = sum(1 .- isnan.(ref)) ÷ 2


    for t = 1:N

        for i = 1:Nx

            for j = 1:Ny

                if !isnan(ref[1, i, j]) && !isnan(ref[2, i, j])

                    err[t] += (X[1, i, j, t] - ref[1, i, j])^2 + (X[2, i, j, t] - ref[2, i, j])^2

                end

            end
        end

        err[t] = sqrt(err[t] / N_points)

    end


    return err

end


## Vorticity

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

function vorticity_profile(x, y, x0, y0)

    η = 2.5e-2

    r2 = ((x - x0)^2 + (y - y0)^2) / η^2

    if r2 >= 1.0
        return NaN
    else
        return 4 / 3 * (1 - r2)
    end

end



function compute_error(Ω, xs, ys)

    err = zeros(size(Ω, 3))


    for k = 1:size(Ω, 3)

        cnt = 0

        for i = 1:nx

            for j = 1:ny
                val = vorticity_profile(1e-2 * xs[i], 1e-2 * ys[j], 0, 0)

                if ~isnan(val)

                    if k > 200

                        val = -val

                    end
                    err[k] += (Ω[i, j, k] - val)^2
                    cnt += 1
                end

            end

        end

        err[k] = sqrt(err[k] / cnt)

    end

    return err
end



function generate_vorticity_figure(Ω, X, err, ω_ref, xs, ys; size_pt=(100, 100))

    f = Figure(layout=(2, 1), size=size_pt)

    N = 5


    grid = f[1, 1] = GridLayout()
    grid_err = f[2, 1]

    ax = Axis(grid_err, xlabel="time (s)", ylabel="error (s⁻¹)", height=size_pt[2] / 6)

    lines!(ax, 0:Δt:(Δt*(size(Ω, 3)-1)), err)


    xlims!(0, 200)

    yline = vlines!(100.0, color=:red, linestyle=:dash, linewidth=1.5, label="transition")

    axislegend(ax, position=:rt, padding=(5, 5, 0, 0), patchsize=(30, 12))


    knots = LinRange(1, size(Ω, 3), N + 1)
    idxs = Int[Int(round((knots[i] + knots[i+1]) / 2)) for i in 1:N]


    hm = nothing
    hm2 = nothing

    max_val = 1.5
    min_val = -1.5

    ts = []


    for i = 1:N

        ax = nothing

        t = idxs[i]

        ax = Axis(grid[1, i], aspect=1.0, xticks=[-5, 0, 5], yticks=[-5, 0, 5])
        xlims!(xs[1], xs[end])
        ylims!(ys[1], ys[end])
        push!(ts, t * Δt)

        if t * Δt < 100
            Makie.heatmap!(xs, ys, ω_ref, colormap=:balance, colorrange=(-1.5, 1.5), interpolate=true)

        elseif t * Δt == 100.0
            ω_ref_comp = zeros(length(xs), length(ys))

            for i = 1:length(xs)
                for j = 1:length(ys)

                    if i <= length(xs) / 2
                        ω_ref_comp[i, j] = ω_ref[i, j]
                    else
                        ω_ref_comp[i, j] = -ω_ref[i, j]
                    end

                end
            end
            Makie.heatmap!(xs, ys, ω_ref_comp, colormap=:balance, colorrange=(-1.5, 1.5), interpolate=true)
        else
            Makie.heatmap!(xs, ys, -ω_ref, colormap=:balance, colorrange=(-1.5, 1.5), interpolate=true)
        end

        hideydecorations!(ax, ticks=true, grid=false)
        hidexdecorations!(ax, ticks=true, grid=false)

        ax = Axis(grid[2, i], aspect=1.0, xticks=[-5, 0, 5], yticks=[-5, 0, 5])
        xlims!(xs[1], xs[end])
        ylims!(ys[1], ys[end])

        hideydecorations!(ax, ticks=true, grid=false)
        hidexdecorations!(ax, ticks=true, grid=false)

        hm = heatmap!(xs, ys, Ω[:, :, t], colormap=:balance, colorrange=(min_val, max_val))

        ax = nothing

        if i == N
            ax = Axis(grid[3, i], xlabel="position, x (cm)", ylabel="position, y (cm)", xticks=[-5, 0, 5], yticks=[-5, 0, 5], aspect=1.0, yaxisposition=:right)
            xlims!(xs[1], xs[end])
            ylims!(ys[1], ys[end])

            #hideydecorations!(ax, ticks = true, grid = false)

        else
            ax = Axis(grid[3, i], aspect=1.0, xticks=[-5, 0, 5], yticks=[-5, 0, 5])
            xlims!(xs[1], xs[end])
            ylims!(ys[1], ys[end])

            hideydecorations!(ax, ticks=true, grid=false)
            hidexdecorations!(ax, ticks=true, grid=false)
        end

        X_mag = sqrt.(X[1, :, :, t] .^ 2 + X[2, :, :, t] .^ 2)

        U_t = X[1, :, :, t]
        V_t = X[2, :, :, t]


        hm2 = heatmap!(xs, ys, X_mag, colormap=:plasma, colorrange=(0, 1))

        u = scale(interpolate(U_t, BSpline(Constant())), (xs, ys))
        v = scale(interpolate(V_t, BSpline(Constant())), (xs, ys))

        F = (x, y) -> Point2(u(x, y), v(x, y))

        streamplot!(F, xs[1] .. xs[end], ys[1] .. ys[end], stepsize=15e-3, colormap=:grays, color=color = p -> RGBAf(0, 0, 0), linewidth=0.5, arrow_size=3.5, gridsize=(20, 20))


    end

    Label(f[1, 1, TopLeft()], "A", font=:bold, fontsize=12, tellwidth=false, tellheight=false, padding=(0, 20, 10, 0), halign=:right)
    Label(grid_err[1, 1, TopLeft()], "B", font=:bold, fontsize=12, tellwidth=false, tellheight=false, padding=(0, 20, 30, 0), halign=:right)

    Colorbar(grid[1:2, end, Right()], hm, label="vorticity (s⁻¹)", vertical=true, flipaxis=true, size=5, halign=:right, alignmode=Outside(2), height=Relative(0.99))
    Colorbar(grid[3, 1, Bottom()], hm2, label="velocity (cms⁻¹)", vertical=false, flipaxis=false, size=5, alignmode=Outside(0), width=Relative(1.0))

    Label(grid[1, 1, Left()], "Reference\n vorticity", font=:bold, rotation=pi / 2, fontsize=10, padding=(-50, 0, 0, 0), tellwidth=false, tellheight=true)
    Label(grid[2, 1, Left()], "Estimated\n vorticity", font=:bold, rotation=pi / 2, fontsize=10, padding=(-50, 0, 0, 0), tellwidth=false, tellheight=true)
    Label(grid[3, 1, Left()], "Measured\n velocity", font=:bold, rotation=pi / 2, fontsize=10, padding=(-50, 0, 0, 0), tellwidth=false, tellheight=true)

    for (i, t) in enumerate(ts)
        if i == 1
            Label(grid[1, i, Top()], "Time: $t s", font=:bold, fontsize=10, padding=(0, 0, 15, 0), tellwidth=true, tellheight=false)
        else
            Label(grid[1, i, Top()], "$t s", font=:bold, fontsize=10, padding=(0, 0, 15, 0), tellwidth=true, tellheight=false)
        end
    end

    rowgap!(f.layout, 10)


    rowgap!(grid, 3)

    colgap!(grid, 3)

    f
end

##

file_dict = JLD.load("../dataset/experiments/vorticity.jld")


X = file_dict["vec_field"]
X = reshape(X, 2 * 32 * 32, :)


nx = 32
ny = 32

xs = LinRange(-5, 5, nx)
ys = LinRange(-5, 5, ny)

Δx = 1e-2 * (xs[2] - xs[1])
Δy = 1e-2 * (ys[2] - ys[1])

V = vorticity_matrix_using_least_squares(nx, ny, Δx, Δy)

Ω = reshape(V * X, nx, ny, :)

err = zeros(size(Ω, 3))

err = compute_error(Ω, xs, ys)

ω_ref = zeros(length(xs), length(ys))

for i = 1:length(xs)
    for j = 1:length(ys)
        ω_ref[i, j] = vorticity_profile(1e-2 * xs[i], 1e-2 * ys[j], 0, 0)
    end

end


size_mm = (178, 150) # 246, 170 # 110
size_inches = size_mm ./ 25.4
size_pt = 72 .* size_inches

f = generate_vorticity_figure(Ω, 100 * reshape(X, 2, 32, 32, :), err, ω_ref, xs, ys; size_pt=size_pt)

path = joinpath(@__DIR__, "figures/")
mkpath(path)

save(joinpath(path, "vorticity.pdf"), f; pt_per_unit=1)

f
