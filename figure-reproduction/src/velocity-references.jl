using Makie.GeometryBasics


function generate_annulus(xc, yc, ri, ro)

    outer = Point2f[]
    inner = Point2f[]

    for θ = -π:0.1:π
        push!(outer, Point2f(ro * cos(θ) + xc, ro * sin(θ)+yc)) 
        push!(inner, Point2f(ri * cos(θ)+xc, ri * sin(θ)+yc))
    end

    annulus = GeometryBasics.Polygon(inner, [outer])
end


function rectangle_with_arrow(x, y, a, b, α; γ = 0, orientation = :forward)
    color = Makie.categorical_colors(:plasma, 256)[end]
    M = [cos(α) -sin(α); sin(α) cos(α)]
    ps = [[x + a/2, y + b/2], [x + a/2, y - b/2], [x - a/2, y - b/2], [x - a/2, y + b/2]]
    ps = [M * x for x in ps]
    ps = Point2f.(ps)
    poly!(ps, color = color, strokecolor = :black, strokewidth = 1.0)


    start_point = (ps[1] + ps[2]) / 2 
    end_point = (ps[3] + ps[4]) / 2

    τ = (end_point - start_point) 

    end_point = end_point - (1-γ) * τ
    start_point = start_point + (1-γ) * τ


    lines!([start_point, end_point] , linewidth=1.5, color = :black)

    if orientation == :forward
        τ = (end_point - start_point) / norm(end_point - start_point)
        β = atan(τ[2], τ[1])
        scatter!(end_point, markersize=7, color = :black, marker=:rtriangle, rotation = β)
    elseif orientation == :backward

        τ = (start_point - end_point) / norm(start_point - end_point)
        β = atan(τ[2], τ[1])
        scatter!(start_point, markersize=7, color = :black, marker=:rtriangle, rotation = β)
    else 
        error("Invalid orientation") 
    end





end


function annulus_with_arrow(xc, yc, α,r1, r2; orientation = :forward, θ_offset = 0)
    r = (r1 + r2) / 2
    g = θ -> [xc + r * cos(θ), yc + r * sin(θ)]
    θ = 2π:-0.01:pi/4
    θ = θ .+ θ_offset

    M = [cos(α) -sin(α); sin(α) cos(α)]
    ps = g.(θ)
    rot = (x) -> M * x
    ps = rot.(ps)
    ps = [Point2(x...) for x in ps]
    color = Makie.categorical_colors(:plasma, 256)[end]
    annulus = generate_annulus(M*[xc, yc]..., 1f0, 2.5f0)
    poly!(annulus, color = color, strokecolor = :black, strokewidth = 1.0)
    lines!(ps, linewidth=1.5, color = :black)
    if orientation == :forward
        τ = (ps[end] - ps[end-1]) / norm(ps[end] - ps[end-1])

        β = atan(τ[2], τ[1])

        scatter!(ps[end], markersize=7, color = :black, marker=:rtriangle, rotation = β)

    elseif orientation == :backward
        τ = (ps[1] - ps[2]) / norm(ps[1] - ps[2])
        β = atan(τ[2], τ[1])
        scatter!(ps[1], markersize=7, color = :black, marker=:rtriangle, rotation =  β )
    else 
        error("Invalid orientation") 
    end

end

function reference!(i)

    if i == 1
        annulus_with_arrow(-2.5, 0, pi/4, 1, 2.5, orientation = :forward, θ_offset = -pi/4)
        annulus_with_arrow(2.5, 0, pi/4, 1, 2.5, orientation = :backward, θ_offset = pi)
    elseif i == 2
        rectangle_with_arrow(0, 0, sqrt(2) * 10, 2.0, pi/4, γ = 0.7)
    elseif i == 3
        rectangle_with_arrow(3.0, 1.5, 4.0, 2.0, 0, γ = 0.7, orientation = :backward)
        rectangle_with_arrow(-3.0, -1.5, 4.0, 2.0, 0, γ = 0.7)
    elseif i == 4
        rectangle_with_arrow(3.0, 0.0, 4.0, 2.0, 0, γ = 0.7, orientation = :backward)
        rectangle_with_arrow(-3.0, 0.0, 4.0, 2.0, 0, γ = 0.7, orientation = :forward)
        rectangle_with_arrow(3.0, 0.0, 4.0, 2.0, pi/2, γ = 0.7, orientation = :forward)
        rectangle_with_arrow(-3.0, 0.0, 4.0, 2.0, pi/2, γ = 0.7, orientation = :backward)

    end
end
