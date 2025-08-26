using Makie
using CairoMakie
using JLD
using Statistics
using Interpolations
using ImageFiltering, ImageContrastAdjustment
using LinearAlgebra
using ColorSchemes

Δt = 0.5

CairoMakie.activate!(type = "svg")

set_theme!(fontsize = 7)

include("src/velocity-references.jl")

function timelapse(video, α = 0.1)

    out = zeros(size(video, 1), size(video, 2))

    for i in axes(video, 3)
        
        out .= (1 - α) * out .+ α * equalize!(video[:, :, i])
    end

    out
end

function equalize!(img::Array{<:UInt8}; nbins::Int = 256, minval::AbstractFloat = 0.0, maxval::AbstractFloat = 1.0)
    adjust_histogram!(img, Equalization(nbins = nbins, minval = minval, maxval = maxval))
end

function equalize!(img::Array{<:AbstractFloat}; nbins::Int = 256, minval::AbstractFloat = 0.0, maxval::AbstractFloat = 1.0)
    adjust_histogram!(img, Equalization(nbins = nbins, minval = minval, maxval = maxval))
end

function sharpen!(img::AbstractArray{<:AbstractFloat}; amount::Float64 = 3.0, sigma::Float64 = 1.0)
    blurred = imfilter(img, Kernel.gaussian(sigma))
    img .= img .+ amount .* (img .- blurred)
    clamp!(img, 0.0, 1.0)
end

function sharpen(img::Array{T}; amount::T = 3.0, sigma::T = 1.0) where T <: AbstractFloat
    blurred = imfilter(img, Kernel.gaussian(sigma))
    sharpened = img .+ amount .* (img .- blurred)
    clamp!(sharpened, 0.0, 1.0)
    return sharpened
end


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


function generate_frame_streamlines_with_timelapse_and_error(Y, ref, frame, err, xs, ys; max_val = 1.5, size_pt = (100, 100), labels=nothing)

    M = size(Y, 4)

    color = p -> RGBAf(0, 0, 0)

    f = Figure(layout = (5, M), size = size_pt, figure_padding = 5)

  
    hm = nothing 

    for i in 1:M

        Y_mag = sqrt.(Y[1, :, :, i].^2 + Y[2, :, :, i].^2)

        U_t = Y[1, :, :, i]
        V_t = Y[2, :, :, i]

        u = scale(interpolate(U_t, BSpline(Linear())), (xs, ys))
        v = scale(interpolate(V_t, BSpline(Linear())), (xs, ys))
    
        F = (x, y) -> Point2(u(x, y), v(x, y))

        ax1 = Makie.Axis(f[1, i],  xticks = LinearTicks(6), yticks=  LinearTicks(6), aspect = 1.0)

        hideydecorations!(ax1, ticks = true, grid = false)
        hidexdecorations!(ax1, ticks = true, grid = false)

        strength = vec(sqrt.(ref[1, :, :, i] .^ 2 .+ ref[2, :, :, i] .^ 2))
    

        reference!(i) # Plot the reference vector field hardcoded in references.jl

        xlims!(-5, 5)
        ylims!(-5, 5)
    
        ax2 =  Makie.Axis(f[2, i], xticks =  LinearTicks(6), yticks=  LinearTicks(6), aspect = 1.0)

        hideydecorations!(ax2, ticks = true, grid = false)
        hidexdecorations!(ax2, ticks = true, grid = false)

        linkxaxes!(ax1, ax2)
    
        hm = Makie.heatmap!(xs, ys, Y_mag; interpolate=false, colorrange = (0.0, max_val), colormap = :plasma)


        streamplot!(F, xs[1]..xs[end], ys[1]..ys[end], stepsize=15e-3, colormap=:grays, color=color, linewidth=0.5, arrow_size=3.5, gridsize=(20, 20))

        xlims!(-5, 5)
        ylims!(-5, 5)

        ax3 = nothing

        if i == 4
            ax3 = Makie.Axis(f[3, i], xticks = LinearTicks(6), yticks= LinearTicks(6), aspect = 1.0, xlabel = "position, x (cm)",  ylabel = "position, y (cm)", yaxisposition = :right)

        else
            ax3 = Makie.Axis(f[3, i], xticks = LinearTicks(6), yticks= LinearTicks(6), aspect = 1.0)

            hideydecorations!(ax3, ticks = true, grid = false)
            hidexdecorations!(ax3, ticks = true, grid = false)

        end
    

        linkxaxes!(ax1, ax3)
    
        Makie.image!(xs[1]..xs[end], ys[1]..ys[end], rotr90(frame[:, :, i]))
    
        xlims!(-5, 5)
        ylims!(-5, 5)

    end


    Colorbar(f[1:2, 4, Right()], hm, label="velocity magnitude (cm/s)", vertical = true, flipaxis = true, size = 5, halign = :left, alignmode=Outside(-5.0))

    Label(f[1, 1:4, Left()], "Reference", font = :bold,  rotation = pi/2,  fontsize = 10, padding = (0, 0, 0, 0), tellwidth = false, tellheight = true)

    Label(f[2, 1:4, Left()], "Measured\n (at 40 s)", font = :bold,  rotation = pi/2,  fontsize = 10, word_wrap = false, padding = (0, 0, 0, 0),tellwidth = false, tellheight = true)

    Label(f[3, 1:4, Left()], "Long-exposure\n shot\n (20–40 s)", font = :bold,  rotation = pi/2,  fontsize = 10, word_wrap = false,  padding = (0, 0, 0, 0), tellwidth = false, tellheight = true)

    Label(f[1, 1, TopLeft()], "A", font = :bold, fontsize = 12, tellwidth = false, tellheight = false, padding = (0, 5, 5, 0), halign = :right)

    Label(f[1:3, 1, Top()], "Two Vortices", font = :bold, fontsize = 10, tellwidth = false, tellheight = true, padding = (0, 0, 5, 0))
    Label(f[1:3, 2, Top()], "Jet", font = :bold, fontsize = 10, tellwidth = false, tellheight = true, padding = (0, 0, 5, 0))
    Label(f[1:3, 3, Top()], "Sides", font = :bold, fontsize = 10, tellwidth = false, tellheight = true, padding = (0, 0, 5, 0))
    Label(f[1:3, 4, Top()], "Cross", font = :bold, fontsize = 10, tellwidth = false, tellheight = true, padding = (0, 0, 5, 0))


    Label(f[4, 1, TopLeft()], "B", font = :bold, fontsize = 12, tellwidth = false, tellheight = false, padding = (0, 5, 30, 0), halign = :right)


    titles = ["Two Vortices", "Jet", "Sides", "Cross"]
    ax = Axis(f[4, :], xlabel = "time (s)", ylabel = "error (cm/s)", height = size_pt[2] / 6)#, title="Error (cm/s)")

    colors = ColorSchemes.okabe_ito

    #linestyles = [:solid, :dash, :dot, :dashdot]
    linestyles = [:solid, :solid, :solid, :solid]

    for i in axes(err, 2)

        T = (length(err[:, i])-1) * Δt

        lines!(ax, 0:Δt:T, err[:, i], label = titles[i], linestyle = (linestyles[i], :dense), linewidth= 2.0, color = colors[i], alpha = 1.0)

        xlims!(0, T)

    end

    axislegend(ax, position = :rt, padding = (5, 5, 0, 0), patchsize = (35, 12), orientation = :horizontal)

    rowgap!(f.layout, 5)
    colgap!(f.layout, -10)

    f
end


X = zeros(2, 32, 32, 4)
ref = zeros(2, 32, 32, 4)
N_err = 81

err = zeros(N_err, 4)

Φ = zeros(4, N_err, 4)
Ψ = zeros(4, N_err, 4)

frame = zeros(1024, 1024, 4)


# Two Vortices
file_dict = JLD.load("../dataset/experiments/two-vortices.jld")
X1 = file_dict["vec_field"]
X[:, :, :, 1] = X1[:, :, :, end]
ref[:, :, :, 1] = file_dict["ref_vec_field"]

err[:, 1] = compute_error(X1, ref[:, :, :, 1])

file_dict_video = JLD.load("../dataset/experiments/two-vortices-video.jld")
video = Float64.(file_dict_video["video"]) ./ 255.0
frame[:, :, 1] = timelapse(video[:, :, end-500:end], 0.1)


# Jet
file_dict = JLD.load("../dataset/experiments/jet.jld")
X2 = file_dict["vec_field"]
X[:, :, :, 2] = X2[:, :, :, end]
ref[:, :, :, 2] = file_dict["ref_vec_field"]

err[:, 2] = compute_error(X2, ref[:, :, :, 2])

file_dict_video = JLD.load("../dataset/experiments/jet-video.jld")
video = Float64.(file_dict_video["video"]) ./ 255.0
frame[:, :, 2] = timelapse(video[:, :, end-500:end], 0.1)


# Sides
file_dict = JLD.load("../dataset/experiments/sides.jld")
X3 = file_dict["vec_field"]
X[:, :, :, 3] = X3[:, :, :, end]
ref[:, :, :, 3] = file_dict["ref_vec_field"]

err[:, 3] = compute_error(X3, ref[:, :, :, 3])

file_dict_video = JLD.load("../dataset/experiments/sides-video.jld")
video = Float64.(file_dict_video["video"]) ./ 255.0
frame[:, :, 3] = timelapse(video[:, :, end-500:end], 0.1)


# Cross
file_dict = JLD.load("../dataset/experiments/cross.jld")
X4 = file_dict["vec_field"]
X[:, :, :, 4] = X4[:, :, :, end]
ref[:, :, :, 4] = file_dict["ref_vec_field"]

err[:, 4] = compute_error(X4, ref[:, :, :, 4])

file_dict_video = JLD.load("../dataset/experiments/cross-video.jld")
video = Float64.(file_dict_video["video"]) ./ 255.0
frame[:, :, 4] = timelapse(video[:, :, end-500:end], 0.1)


# Sharpen the frames
for i = 1:4
    frame[:, :, i] .= sharpen(frame[:, :, i]; amount = 3.0, sigma = 20.0)
end


# m/s -> cm/s
X *= 100
err *= 100
ref *= 100


max_val = max(maximum(sqrt.(X[1, :, :, :, :].^2 + X[2, :, :, :, :].^2)), 1.0)


xs = LinRange(-5, 5, 32)
ys = LinRange(-5, 5, 32)

##

size_mm = (178, 150) # 246, 170 # (160, 130) # (160, 150)
size_inches =  size_mm ./ 25.4
size_pt = 72 .* size_inches

f = generate_frame_streamlines_with_timelapse_and_error(X, ref, frame, err, xs, ys; size_pt = size_pt, max_val = max_val)

path = joinpath(@__DIR__, "figures/")
mkpath(path)

save(joinpath(path, "velocity.pdf"), f; pt_per_unit = 1)
