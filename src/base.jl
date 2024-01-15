export Dimension, Interval, AbstractQuadrule, QuadratureRule

const Dimension = Integer

struct Interval{T<:Real}
    a::T
    b::T
end

abstract type AbstractQuadrule{D}
end

struct QuadratureRule{D,V,W} <: AbstractQuadrule{D}
    x::V
    w::W
    function QuadratureRule(x::V, w::W) where {V,W}
        D = ndims(x)
        @assert length(x) == length(w) "Number of quadrature points and weights are not consistent."
        return new{D,V,W}(x,w)
    end
end