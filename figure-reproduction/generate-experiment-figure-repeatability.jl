using Makie
using CairoMakie
using JLD
using Statistics
using LinearAlgebra
using ColorSchemes


Δt = 0.5

base_dir = "../dataset/experiments/repeatability/"


function compute_error(X, ref)

    N = size(X, 4)

    Nx = size(X, 2)
    Ny = size(X, 3)

    err = zeros(N)

    N_points = sum(1 .-isnan.(ref)) ÷ 2


    for t = 1:N

        for i =  1:Nx

            for j = 1:Ny

                if !isnan(ref[1, i, j]) && !isnan(ref[2, i, j])

                    err[t] += (X[1, i, j, t] - ref[1, i, j])^2 + (X[2, i, j, t] - ref[2, i, j])^2

                end
    
            end
        end

        err[t] = sqrt(err[t]/N_points)

    end

    return err

end


CairoMakie.activate!(type = "svg")

set_theme!(fontsize = 7)
##

repeats = 5

X = zeros(2, 32, 32, 4)
N_err = 81

err = zeros(N_err, repeats, 4)


# First
for i = 1:repeats
    file_dict = JLD.load(joinpath(base_dir, "two-vortices-$i.jld"))
    vec_field = file_dict["vec_field"]
    ref = file_dict["ref_vec_field"]
    err[:, i, 1] = compute_error(vec_field, ref)
end


# Second

for i = 1:repeats
    file_dict = JLD.load(joinpath(base_dir, "jet-$i.jld"))
    vec_field = file_dict["vec_field"]
    ref = file_dict["ref_vec_field"]
    err[:, i, 2] = compute_error(vec_field, ref)
end


# Third

for i = 1:repeats
    file_dict = JLD.load(joinpath(base_dir, "sides-$i.jld"))
    vec_field = file_dict["vec_field"]
    ref = file_dict["ref_vec_field"]
    err[:, i, 3] = compute_error(vec_field, ref)
end


# Fourth


for i = 1:repeats
    file_dict = JLD.load(joinpath(base_dir, "cross-$i.jld"))
    vec_field = file_dict["vec_field"]
    ref = file_dict["ref_vec_field"]
    err[:, i, 4] = compute_error(vec_field, ref)
end

err *= 100

mean_err = mean(err, dims = 2)
mean_err = dropdims(mean_err, dims=2)

##

size_mm = (178, 80) # 246, 170 # 110
size_inches =  size_mm ./ 25.4
size_pt = 72 .* size_inches


f = Figure(size = size_pt)

ax1 = Axis(f[1, 1], xlabel = "time (s)", ylabel = "error (cm/s)", title="Two Vortices")
xlims!(ax1, 0, 40)
ax2 = Axis(f[1, 2], xlabel = "time (s)", ylabel = "error (cm/s)", title="Jet")
xlims!(ax2, 0, 40)
ax3 = Axis(f[2, 1], xlabel = "time (s)", ylabel = "error (cm/s)", title="Sides")
xlims!(ax3, 0, 40)
ax4 = Axis(f[2, 2], xlabel = "time (s)", ylabel = "error (cm/s)", title="Cross")
xlims!(ax4, 0, 40)


for i = 1:repeats
    lines!(ax1, (0:Δt:(N_err-1)*Δt), err[:, i, 1], color = (:blue, 0.2))
    lines!(ax1, (0:Δt:(N_err-1)*Δt), mean_err[:, 1], color = :blue, linewidth = 1)
    lines!(ax2, (0:Δt:(N_err-1)*Δt), err[:, i, 2], color = (:blue, 0.2))
    lines!(ax2, (0:Δt:(N_err-1)*Δt), mean_err[:, 2], color = :blue, linewidth = 1)
    lines!(ax3, (0:Δt:(N_err-1)*Δt), err[:, i, 3],  color = (:blue, 0.2))
    lines!(ax3, (0:Δt:(N_err-1)*Δt), mean_err[:, 3], color = :blue, linewidth = 1)
    lines!(ax4, (0:Δt:(N_err-1)*Δt), err[:, i, 4], color = (:blue, 0.2)) 
    lines!(ax4, (0:Δt:(N_err-1)*Δt), mean_err[:, 4], color = :blue, linewidth = 1)
end

rowgap!(f.layout, 5)

path = joinpath(@__DIR__, "figures/")
mkpath(path)

save(joinpath(path, "repeatability.pdf"), f; pt_per_unit = 1)