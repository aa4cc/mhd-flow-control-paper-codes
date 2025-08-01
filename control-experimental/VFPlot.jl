#=
Plotting vector fields
=#

using GLMakie

function plot_field(vec_field::Array{<:AbstractFloat}, ax::T, fig; video::Bool = true) where {T<:Union{GLMakie.Axis,GLMakie.Scene}}
    scale = (100e-3) / 1024
    vec_field .*= scale
    
    (_, width, height) = size(vec_field)
    xs = LinRange(0, 100e-3, width)
    ys = LinRange(100e-3, 0, height)
    us = vec_field[1, :, :]
    vs = vec_field[2, :, :]

    us = permutedims(us, [2, 1])
    vs = permutedims(vs, [2, 1])
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    us = GLMakie.Observable(us)
    vs = GLMakie.Observable(vs)

    ax.backgroundcolor = "black"
    arr = arrows!(ax, xs, ys, us, vs, arrowsize = 5, lengthscale = 100e-3,
    arrowcolor = strength, linecolor = strength, colorrange=(0, 0.016), normalize=false)
    Colorbar(fig[1, 2], arr)
    if(video)
        function update(vec_field::Array{<:Real})
            
            us_u = permutedims(vec_field[1, :, :].*scale, [2, 1])
            vs_u = permutedims(vec_field[2, :, :].*scale, [2, 1])
            
            # remove this line later ... (this should be done elsewhere == interpolation)
            strength_u = vec(sqrt.(us_u .^ 2 .+ vs_u .^ 2))
            us[] = us_u
            vs[] = vs_u
            arr.arrowcolor = strength_u
            arr.linecolor = strength_u
        end
        return update
    else
        #wait(display(ax))
        arr
    end
end

function plot_field(vec_field::Array{<:AbstractFloat}; plt_scale::Number = 1, video::Bool = true)
    (_, width, height) = size(vec_field)
    xs = LinRange(0, width * 2 * plt_scale, width)
    ys = LinRange(height * 2 * plt_scale, 0, height)
    us = vec_field[1, :, :] .* plt_scale
    vs = vec_field[2, :, :] .* plt_scale
    # should not really be doing this
    #us[abs.(us) .>(8*plt_scale) ] .= 0
    #vs[abs.(vs) .>(8*plt_scale) ] .= 0
    
    us = permutedims(us, [2, 1])
    vs = permutedims(vs, [2, 1])
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    us = GLMakie.Observable(us)
    vs = GLMakie.Observable(vs)
    
    fig = GLMakie.Figure(size=(600, 600))
    ax = GLMakie.Axis(fig[1, 1])

    ax.backgroundcolor = "black"
    arr = arrows!(ax, xs, ys, us, vs, arrowsize = 5, lengthscale = 100e-3,colorrange=(0, 0.016),
    arrowcolor = strength, linecolor = strength)
    if(video)
        function update(vec_field::Array{<:AbstractFloat})
            us_u = permutedims(vec_field[1, :, :] .* plt_scale, [2, 1])
            vs_u = permutedims(vec_field[2, :, :] .* plt_scale, [2, 1])
            # remove this line later ... (this should be done elsewhere == interpolation)
            #us_u[abs.(us_u) .>(8 * plt_scale) ] .= 0
            #vs_u[abs.(vs_u) .>(8 * plt_scale) ] .= 0
            strength_u = vec(sqrt.(us_u .^ 2 .+ vs_u .^ 2))
            us[] = us_u
            vs[] = vs_u
            arr.arrowcolor = strength_u
            arr.linecolor = strength_u
        end
        display(fig)
        return update
    else
        wait(display(fig))
    end
end

function make_gif(v_field::Array{<:AbstractFloat}; framerate=30, filename="measurements/measurement1.mp4")
    scale = (100e-3) / 1024
    
    vec_field = copy(v_field) .*scale
    
    (_, width, height, frames) = size(vec_field)
    xs = LinRange(0, 100e-3, width)
    ys = LinRange(100e-3, 0, height)
    us = vec_field[1, :, :, 1] 
    vs = vec_field[2, :, :, 1] 
    @info "$frames frames to animate"

    us = permutedims(us, [2, 1])
    vs = permutedims(vs, [2, 1])
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    us = GLMakie.Observable(us)
    vs = GLMakie.Observable(vs)
    
    fig = GLMakie.Figure(size=(600, 600))
    ax = GLMakie.Axis(fig[1, 1])
    
    ax.backgroundcolor = "black"
    arr = arrows!(ax, xs, ys, us, vs, arrowsize = 5, lengthscale = 100e-3,
    arrowcolor = strength, linecolor = strength, colorrange=(0, 0.016), normalize=false)
    Colorbar(fig[1, 2], arr)
    @info "starting animation"
    record(fig, filename, 1:frames; framerate=framerate) do i
        us_u = permutedims(vec_field[1, :, :, i], [2, 1])
        vs_u = permutedims(vec_field[2, :, :, i], [2, 1])
        strength_u = vec(sqrt.(us_u .^ 2 .+ vs_u .^ 2))
        us[] = us_u
        vs[] = vs_u
        arr.arrowcolor = strength_u
        arr.linecolor = strength_u
    end
    @info "animation done"
end

function animate_trajectory(Y; framerate=30, filename="tmp.mp4")
 
    (_ ,width, height, N) = size(Y)
    max_val = maximum(sqrt.(Y[1,:, :, :].^2 + Y[2,:, :, :].^2))
    xs = LinRange(0, width * 2, width)
    ys = LinRange(height * 2, 0, height)
 
    y0 = copy(Y[:, :, :, 1])
 
    us = copy(y0[1,:,:])
    vs = copy(y0[2,:,:])
    us = permutedims(us, [2, 1])
    vs = permutedims(vs, [2, 1])
    strength = vec(sqrt.(us .^ 2 + vs .^ 2))
 
 
    us_obs = Observable(y0[1,:,:])
    vs_obs = Observable(y0[2,:,:])
 
 
    f = Figure(size = (800, 800))
    GLMakie.Axis(f[1, 1], backgroundcolor = "black")
 
    arr = arrows!(xs, ys, us_obs, vs_obs, arrowsize = 6, lengthscale = 1,
    arrowcolor = strength, linecolor = strength, colorrange = (0.0, max_val), colormap= :jet, normalize = true)
 
    record(f, filename, 1:N;
        framerate=framerate) do frame
 
        yi = copy(Y[:, :, :, frame])
 
        us = yi[1, :, :]
        vs = yi[2, :, :]
        us = permutedims(us, [2, 1])
        vs = permutedims(vs, [2, 1])
 
 
        strength = vec(sqrt.(us .^ 2 + vs .^ 2))
 
        us_obs[] = us
        vs_obs[] = vs
 
        arr.linecolor = strength
        arr.arrowcolor = strength
    end
end

function get_flow_line(vec_field, x, y)
    # didnt test this
    (_, width, height) = size(vec_field)
    xs = Float32[]
    ys = Float32[]
    cx, cy = x, y
    visited = Set([[cx, cy]])
    push!(xs, cx)
    push!(ys, cy)
    while((0 < cx <= width) && (0 < cy <= height))
        cx += vec_field[1, cx, cy]
        cy += vec_field[2, cx, cy]
        if([cx, cy] in visited)
            break
        end 
        push!(xs, cx)
        push!(ys, cy)
        push!(visited, [cx, cy])
    end
    return xs, ys
end


if abspath(PROGRAM_FILE) == @__FILE__
    using Images
    #include("PIV.jl")
    # -------------load images ------------
    #img_a = permutedims(convert(Array{Float32}, load("img_0_test.png")), [2, 1])
    img_a = convert(Array{Float32}, load("img_0.png"))
    img_b = convert(Array{Float32}, load("img_1.png"))
    
    imgs = [img_a, img_b, img_a];
    # ------------- initialize PIV --------------
    data = PIVData(1024, 1024, 32)
    vec_field = zeros(2, 32, 32)
    # ------------ initialize PLOT --------------
    compute_vector_field!(data, vec_field, img_a, img_b)
    update_vf = plot_field(vec_field, video=false)
    #=
    for i = 1:50
        for j = 1:length(imgs) - 1
            compute_vector_field!(data, vec_field, imgs[j], imgs[j + 1])
            update(vec_field)
            sleep(1/10)
        end
    end
    =#
end