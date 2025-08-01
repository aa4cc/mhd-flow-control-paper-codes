using ImageContrastAdjustment, Statistics

function equalize!(img::Array{<:AbstractFloat}; nbins::Int = 256, minval::AbstractFloat = 0.0, maxval::AbstractFloat = 1.0)
    adjust_histogram!(img, Equalization(nbins = nbins, minval = minval, maxval = maxval))
end

function contrast_stretch!(img::Array{<:AbstractFloat}; t::AbstractFloat = 0.5, slope::Number = 1.0)
    adjust_histogram!(img, ContrastStretching(t = t, slope = slope))
end

function linear_stretch!(img::Array{<:AbstractFloat}; src_minval::AbstractFloat=0.1, src_maxval::AbstractFloat=0.9, dst_minval::AbstractFloat=0.05, dst_maxval::AbstractFloat=0.95)
    adjust_histogram!(img, LinearStretching(src_minval=src_minval, src_maxval=src_maxval, dst_minval=dst_minval, dst_maxval=dst_maxval))
end

function equalize!(img::Array{<:UInt8}; nbins::Int = 256, minval::AbstractFloat = 0.0, maxval::AbstractFloat = 1.0)
    adjust_histogram!(img, Equalization(nbins = nbins, minval = minval, maxval = maxval))
end

function contrast_stretch!(img::Array{<:UInt8}; t::AbstractFloat = 0.5, slope::Number = 1.0)
    adjust_histogram!(img, ContrastStretching(t = t, slope = slope))
end

function linear_stretch!(img::Array{<:UInt8}; src_minval::AbstractFloat=0.1, src_maxval::AbstractFloat=0.9, dst_minval::AbstractFloat=0.05, dst_maxval::AbstractFloat=0.95)
    adjust_histogram!(img, LinearStretching(src_minval=src_minval, src_maxval=src_maxval, dst_minval=dst_minval, dst_maxval=dst_maxval))
end


function intensity_cap!(img::Array{<:AbstractFloat}; n::Number = 1)
    std = Statistics.std(img)
    median = Statistics.median(img)
    max_intensity = median + n * std
    img[img .> max_intensity] .= max_intensity
end