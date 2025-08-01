using BenchmarkTools
using FFTW
using Images
using Statistics
using Base.Threads
using AppleAccelerate

abstract type PIVType end

struct Default <: PIVType end
struct Overlayed <: PIVType end


struct PIVData{C <: PIVType}
    fft::AbstractFFTs.Plan
    ifft::AbstractFFTs.Plan
    window_size::Tuple{Int, Int}
    overlay::Tuple{Int, Int}
    rows::Int
    cols::Int
    sub_pixel_fnc::Union{Function, Nothing}
end

const OverlayedPIV = PIVData{Overlayed}
const DefaultPIV = PIVData{Default}

#=
Macro to turn on or off threads based on a flag
=#
macro dothreads(flag, expr)
    quote
        if $(flag)
            Base.Threads.@threads $expr
        else
            $expr
        end
    end |> esc
end
# FFTW is most effective when run on single thread -> multiple (single thread)FFTW calculations on different threads is fast 
FFTW.set_num_threads(1)

function PIVData(width::Int, height::Int, window_size::T; sub_pixel_function::String = "") where {T<:Union{Int, Tuple{Int, Int}}}
    # Inits PIV, FFTW finds best FFT algorithm suitable for current hardware
    if(typeof(window_size) == Int)
        window_size = (window_size, window_size)
    end
    @assert mod(width, window_size[2]) ≈ 0 && mod(height, window_size[1]) ≈ 0

    frame1 = randn(Float32, width, height)
    frame2 = randn(Float32, width, height)

    rows = height ÷ window_size[1]
    cols = width ÷ window_size[2]

    window1 = view(frame1, 1:window_size[1], 1:window_size[2]) 
    window2 = view(frame2, 1:window_size[1], 1:window_size[2])
    rfft_p = plan_rfft(window1, flags=FFTW.ESTIMATE)
    fft_a = rfft_p*window1
    fft_b = rfft_p*window2
    cross_corr = conj.(fft_a) .* fft_b
    ifft_p = plan_irfft(cross_corr,window_size[1], flags=FFTW.ESTIMATE)
    sub_pixel_fnc = nothing
    if(sub_pixel_function == "PARAB")
        sub_pixel_fnc = parabolic_peak_fit
    end

    PIVData{Default}(rfft_p, ifft_p, window_size, (0, 0), rows, cols, sub_pixel_fnc)
end

function OverlayedPIVData(width::Int, height::Int, window_size::T, overlay::T; sub_pixel_function::String = "") where {T<:Union{Int, Tuple{Int, Int}}}
    # Inits PIV, FFTW finds best FFT algorithm suitable for current hardware
    if(typeof(window_size) == Int)
        window_size = (window_size, window_size)
    end
    if(typeof(overlay) == Int)
        overlay = (overlay, overlay)
    end
    @assert mod(width, window_size[2]) ≈ 0 && mod(height, window_size[1]) ≈ 0
    @assert all(window_size .> overlay)
    @assert all(overlay .> 0)

    frame1 = randn(Float32, width, height)
    frame2 = randn(Float32, width, height)

    rows = height ÷ window_size[1]
    cols = width ÷ window_size[2]

    rows = UInt16((height - window_size[1]) / overlay[1]) + 1
    cols = UInt16((width - window_size[2])/ overlay[2]) + 1

    window1 = view(frame1, 1:window_size[1], 1:window_size[2]) 
    window2 = view(frame2, 1:window_size[1], 1:window_size[2])
    rfft_p = plan_rfft(window1, flags=FFTW.ESTIMATE)
    fft_a = rfft_p*window1
    fft_b = rfft_p*window2
    cross_corr = conj.(fft_a) .* fft_b
    ifft_p = plan_irfft(cross_corr,window_size[1], flags=FFTW.ESTIMATE)
    sub_pixel_fnc = nothing
    if(sub_pixel_function == "PARAB")
        sub_pixel_fnc = parabolic_peak_fit
    end

    PIVData{Overlayed}(rfft_p, ifft_p, window_size, overlay, rows, cols, sub_pixel_fnc)
end

function compute_vector_field!(data::DefaultPIV, vector_field::Array{<:Real}, frame1::Array{<:Real}, frame2::Array{<:Real}; threaded::Bool=true, fps::Int = 30)
    # computes vector field in place  stores it in vector_field
    window_size, rows, cols = data.window_size, data.rows, data.cols
    @dothreads threaded for idx in CartesianIndices((rows, cols))
        r = idx[1] - 1
        c = idx[2] - 1
        @inbounds window1 = view(frame1, r*window_size[1] + 1 : (r+1)*window_size[1], c*window_size[2] + 1: (c+1)*window_size[2])
        @inbounds window2 = view(frame2, r*window_size[1] + 1 : (r+1)*window_size[1], c*window_size[2] + 1: (c+1)*window_size[2])
        #substract mean
        #window1 .= window1 .- mean(window1, dims=:) # This is for normed PIV, propably is not needed (should add division by std)
        #window2 .= window2 .- mean(window2, dims=:)
        fft1 = conj.(data.fft*window1)
        fft2 = data.fft*window2
        cross_corr = fft1 .* fft2
        result = fftshift(data.ifft*cross_corr)
        
        res = argmax(result)
        if(!isnothing(data.sub_pixel_fnc))
            row, col = parabolic_peak_fit(result, res[1], res[2])
        else
            row = res[1] - 1
            col = res[2] - 1
        end

        @inbounds vector_field[1, r+1, c+1] = col - window_size[2] ÷ 2
        @inbounds vector_field[2, r+1, c+1] = window_size[1] ÷ 2 - row
    end
    vector_field .*= fps
end


function compute_vector_field!(data::OverlayedPIV, vector_field::Array{<:Real}, frame1::Array{<:Real}, frame2::Array{<:Real}; threaded::Bool = true, fps::Int = 30)
    # computes vector field in place  stores it in vector_field
    # uses overlaying windows
    window_size, rows, cols, overlay = data.window_size, data.rows, data.cols, data.overlay
    or, oc = overlay
    
    @dothreads threaded for idx in CartesianIndices((rows, cols))
        r = idx[1] - 1
        c = idx[2] - 1
        @inbounds window1 = view(frame1, r*or + 1 : ((r*or)+window_size[1]), c*oc + 1: ((c*oc)+ window_size[2]))
        @inbounds window2 = view(frame2, r*or + 1 : ((r*or)+window_size[1]), c*oc + 1: ((c*oc)+ window_size[2]))

        fft1 = conj.(data.fft*window1)
        fft2 = data.fft*window2
        cross_corr = fft1 .* fft2
        result = fftshift(data.ifft*cross_corr)
        
        res = argmax(result)
        if(!isnothing(data.sub_pixel_fnc))
            row, col = parabolic_peak_fit(result, res[1], res[2])
        else
            row = res[1] - 1
            col = res[2] - 1
        end

        @inbounds vector_field[1, r+1, c+1] = col - window_size[2] ÷ 2
        @inbounds vector_field[2, r+1, c+1] = window_size[1] ÷ 2 - row
    end
    vector_field .*= fps
end

function _err(vec_field::Array{<:Real}, r::Int, c::Int, U::Array{<:Real}, V::Array{<:Real}, ϵ::Number) :: Float32
    # used for outlier rejection
    med_u = median(U)
    med_v = median(V)
    fluct_u = vec_field[1, r, c] - med_u
    fluct_v = vec_field[2, r, c] - med_v
    res_u = U .- med_u
    res_v = V .- med_v
    med_res_u = median(abs.(res_u))
    med_res_v = median(abs.(res_v))
    norm_fluct_u = abs(fluct_u/(med_res_u+ϵ))
    norm_fluct_v = abs(fluct_v/(med_res_v+ϵ))
    sqrt(norm_fluct_u^2 + norm_fluct_v^2)
end

function validate(vec_field::Array{<:Real}; threshold::Number=2, threaded::Bool = true)
    # Outlier rejection algorithm 
    b = 1 # neighbourhood radius
    ϵ = 0.1
    n_neighbours = b==1 ? 8 : 24
    (_, rows, cols) = size(vec_field)
    vec_field_new = copy(vec_field)
    # zero out corners
    for i = 1:4
        r = 1
        c = 1
        if(i in [3, 4])
            r = rows
        end
        if(i in [2, 4])
            c = cols
        end
        vec_field[1, r, c] = 0.0
        vec_field[2, r, c] = 0.0
    end
    @dothreads threaded for i in CartesianIndices((rows - 2b, cols - 2b))
        r = i[1] + b
        c = i[2] + b
        @inbounds U, V = vec_field[1, r-b:r+b, c-b:c+b], vec_field[2, r-b:r+b, c-b:c+b]
        Uf, Vf = U[:], V[:]
        @inbounds U = vcat(Uf[1:(2b+1)b+b], Uf[(2b+1)b+b+2:end])
        @inbounds V = vcat(Vf[1:(2b+1)b+b], Vf[(2b+1)b+b+2:end])
        med_u = median(U)
        med_v = median(V)
        fluct_u = vec_field[1, r, c] - med_u
        fluct_v = vec_field[2, r, c] - med_v
        res_u = U .- med_u
        res_v = V .- med_v
        med_res_u = median(abs.(res_u))
        med_res_v = median(abs.(res_v))
        norm_fluct_u = abs(fluct_u/(med_res_u+ϵ))
        norm_fluct_v = abs(fluct_v/(med_res_v+ϵ))
        if (sqrt(norm_fluct_u^2 + norm_fluct_v^2) > threshold)
            vec_field_new[1, r, c] = sum(U) / n_neighbours
            vec_field_new[2, r, c] = sum(V) / n_neighbours
        end
    end
    #edges - rows
    n_neighbours = b==1 ? 5 : 14
    @dothreads threaded for i = 1 + b:rows - b
        @inbounds U1 = vec_field[1, i-b:i+b, 1:1+b]
        @inbounds U2 = vec_field[1, i-b:i+b, end-b:end]
        @inbounds V1 = vec_field[2, i-b:i+b, 1:1+b]
        @inbounds V2 = vec_field[2, i-b:i+b, end-b:end]
        U1f, V1f = U1[:], V1[:]
        U2f, V2f = U2[:], V2[:]
        @inbounds U1 = vcat(U1f[1 : b], U1f[b+2:end])
        @inbounds V1 = vcat(V1f[1 : b], V1f[b+2:end])
        @inbounds U2 = vcat(U2f[1 : end-b-1], U2f[end-b+1:end])
        @inbounds V2 = vcat(V2f[1 : end-b-1], V2f[end-b+1:end])
        if(_err(vec_field, i, 1, U1, V1, ϵ) > threshold)
            vec_field_new[1, i, 1] = sum(U1) / n_neighbours
            vec_field_new[2, i, 1] = sum(V1) / n_neighbours
        end
        if(_err(vec_field, i, cols, U2, V2, ϵ) > threshold)
            vec_field_new[1, i, end] = sum(U2) / n_neighbours
            vec_field_new[2, i, end] = sum(V2) / n_neighbours
        end  
    end
    #edges - cols
    @dothreads threaded for i = 1 + b:cols - b
        @inbounds U1 = vec_field[1, 1:1+b, i-b:i+b]
        @inbounds U2 = vec_field[1, end-b:end, i-b:i+b]
        @inbounds V1 = vec_field[2, 1:1+b, i-b:i+b]
        @inbounds V2 = vec_field[2, end-b:end, i-b:i+b]
        U1f, V1f = U1[:], V1[:]
        U2f, V2f = U2[:], V2[:]
        @inbounds U1 = vcat(U1f[1 : 2*b], U1f[2*b+2:end])
        @inbounds V1 = vcat(V1f[1 : 2*b], V1f[2*b+2:end])
        @inbounds U2 = vcat(U2f[1 : 2*b+1], U2f[2*b+3:end])
        @inbounds V2 = vcat(V2f[1 : 2*b+1], V2f[2*b+3:end])
        if(_err(vec_field, 1, i, U1, V1, ϵ) > threshold)
            vec_field_new[1, 1, i] = sum(U1) / n_neighbours
            vec_field_new[2, 1, i] = sum(V1) / n_neighbours
        end
        if(_err(vec_field, rows, i, U2, V2, ϵ) > threshold)
            vec_field_new[1, rows, i] = sum(U2) / n_neighbours
            vec_field_new[2, rows, i] = sum(V2) / n_neighbours
        end  
    end
    return vec_field_new
end

function parabolic_peak_fit(cross_corr::Matrix{<:Number}, row::Number, col::Number)
    # sub pixel precision fitting
    cr = size(cross_corr)
    if(row == 1 || row == cr[1]) || (col == 1 || col == cr[2])
        return row - 1, col - 1
    end
    r = row + ((cross_corr[row - 1, col] - cross_corr[row + 1, col])/(2cross_corr[row - 1, col] - 4cross_corr[row, col] + 2cross_corr[row + 1, col]))
    c = col + ((cross_corr[row, col - 1] - cross_corr[row, col + 1])/(2cross_corr[row, col - 1] - 4cross_corr[row, col] + 2cross_corr[row, col + 1]))
    return r - 1, c - 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Images
    
    img_a = convert(Array{Float32}, load("img_0.png"))
    img_b = convert(Array{Float32}, load("img_1.png"))
    #data = OverlayedPIVData(1024, 1024, 32, 16)
    data = PIVData(1024, 1024, 32)
    

    vec_field = zeros(2, 32, 32)
    compute_vector_field!(data, vec_field, img_a, img_b)
    @benchmark validate(vec_field)

end



