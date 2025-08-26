using Makie
using CairoMakie
using JLD
using Statistics
using Interpolations
using ImageFiltering, ImageContrastAdjustment
using LinearAlgebra
using ColorSchemes
using LaTeXStrings


Δt = 0.5

t = 0:Δt:40

CairoMakie.activate!(type = "svg")

set_theme!(fontsize = 7)

Φ = zeros(4, 81, 4)
Ψ = zeros(4, 81, 4)

# Two Vortices
file_dict = JLD.load("../dataset/experiments/two-vortices.jld")
Φ[:, :, 1] = file_dict["phi"]
Ψ[:, :, 1] = file_dict["psi"]

# Jet
file_dict = JLD.load("../dataset/experiments/jet.jld")
Φ[:, :, 2] = file_dict["phi"]
Ψ[:, :, 2] = file_dict["psi"]

# Sides
file_dict = JLD.load("../dataset/experiments/sides.jld")
Φ[:, :, 3] = file_dict["phi"]
Ψ[:, :, 3] = file_dict["psi"]

# Cross
file_dict = JLD.load("../dataset/experiments/cross.jld")
Φ[:, :, 4] = file_dict["phi"]
Ψ[:, :, 4] = file_dict["psi"]

## 

size_mm = (178, 120) # 246, 170 # (160, 130) # (160, 150)
size_inches =  size_mm ./ 25.4
size_pt = 72 .* size_inches

f = Figure(size = size_pt)

#linestyles = [:solid, :dash, :dot, :dashdot]
linestyles = [:solid, :solid, :solid, :solid]

ax = nothing

labels = [L"\phi_1,\,\psi_1", L"\phi_2,\,\psi_2", L"\phi_3,\,\psi_3", L"\phi_4,\,\psi_4"]

colors = ColorSchemes.okabe_ito

ls = []

for i = 1:4
    ax = Axis(f[i, 1], xlabel = "time (s)")

    if i != 4
        hideydecorations!(ax, ticks = true, grid = false)
        hidexdecorations!(ax, ticks = true, grid = false)
    end
    xlims!(0, 40)
    for j = 1:4
        l = lines!(ax, t, Φ[j, :, i], linestyle = linestyles[j], label = labels[j], alpha = 1.0, color = colors[j], linewidth = 2.0)
        if i == 1
            push!(ls, l)
        end
    end

    ax = Axis(f[i, 2], xlabel = "time (s)", yaxisposition = :right)
    if i != 4
        hideydecorations!(ax, ticks = true, grid = false)
        hidexdecorations!(ax, ticks = true, grid = false)
    end
    xlims!(0, 40)
    for j = 1:4
        l = lines!(ax, t, Ψ[j, :, i],  linestyle =  linestyles[j], label = labels[j], alpha = 1.0, color = colors[j], linewidth = 2.0)
    end
end

Label(f[0, 1], "Electrodes", tellwidth = false, tellheight = true, padding = (0, 0, 5, 0), font = :bold, fontsize = 10)
Label(f[0, 2], "Electromagnets", tellwidth = false, tellheight = true, padding = (0, 0, 5, 0), font = :bold, fontsize = 10)

Label(f[1, 0], "Two Vortices", tellwidth = true, tellheight = false, padding = (0, 5, 0, 0), font = :bold, fontsize = 10, rotation = pi/2)
Label(f[2, 0], "Jet", tellwidth = true, tellheight = false, padding = (0, 5, 0, 0), font = :bold, fontsize = 10, rotation = pi/2)
Label(f[3, 0], "Sides", tellwidth = true, tellheight = false, padding = (0, 5, 0, 0), font = :bold, fontsize = 10, rotation = pi/2)
Label(f[4, 0], "Cross", tellwidth = true, tellheight = false, padding = (0, 5, 0, 0), font = :bold, fontsize = 10, rotation = pi/2)

Legend(f[1:3, 2, Right()], ax, "Signals", framevisible = false, labelsize=10)

colgap!(f.layout, 10)
rowgap!(f.layout, 5)

f

save("figures/control-signals.pdf", f; pt_per_unit = 1)