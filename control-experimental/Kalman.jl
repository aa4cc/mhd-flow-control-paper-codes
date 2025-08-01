using JLD
using LinearAlgebra


mutable struct KalmanFilter{T<:Real} 
    xtt::Vector{T}
    Ptt::Matrix{T}
    xt1t::Vector{T}
    Pt1t::Matrix{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    R::Symmetric{T, Matrix{T}}
    Q::Symmetric{T, Matrix{T}}
    P::Matrix{T}
    nd::Int
    npod::Int
end

function data_step!(data::KalmanFilter{T}, y::Vector{T}) where T

    xtt1 = data.xt1t
    Ptt1 = data.Pt1t

    C = data.C
    R = data.R
    P = data.P
    nd = data.nd
    npod = data.npod

    K_inv = Symmetric(P * Ptt1[(nd - 1) * npod + 1:nd * npod, (nd - 1) * npod + 1:nd * npod] * P' + R)

    L = (Ptt1 * C') / K_inv
    data.xtt = xtt1 + L*(y - C*xtt1)
    data.Ptt = Ptt1 - L*C*Ptt1
end

function time_step!(data::KalmanFilter{T}, u::Vector{T}) where T

    xtt = data.xtt
    Ptt = data.Ptt
    Q = data.Q
    A = data.A
    B = data.B

    data.xt1t = A*xtt + B*u
    data.Pt1t = A*Ptt*A' + Q
end
