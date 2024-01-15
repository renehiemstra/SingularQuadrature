export GaussRule, Legendre, Lobatto

abstract type AbstractGaussrule <: AbstractQuadrule{1} end
struct Legendre <: AbstractGaussrule end
struct Lobatto <: AbstractGaussrule end

Base.length(Q::AbstractQuadrule) = length(Q.x)
Base.eltype(Q::AbstractQuadrule) = eltype(Q.x)
Base.ndims(Q::AbstractQuadrule{D}) where D = D

struct GaussRule{S,T<:Real, V<:AbstractVector{T}} <: AbstractQuadrule{1}
    x::V
    w::V
    function GaussRule{S}(x::V, w::V) where {S, V}
        @assert length(x)==length(w) "Number of points are not consistent with number of weights."
        T = eltype(w)
        return new{S,T,V}(x,w)
    end
end

function GaussRule(::Type{Legendre}, n::Dimension)
    x, w = gausslegendre(n)
    return GaussRule{Legendre}(x, w)
end

function GaussRule(::Type{Lobatto}, n::Dimension)
    x, w = gausslobatto(n)
    return GaussRule{Lobatto}(x, w)
end

GaussRule(n::Dimension) = GaussRule(Legendre, n)

function GaussRule(method, n::Dimension, I::Interval)
    Q = GaussRule(method, n)
    return affine_transform!(Q, Interval(-1.0,1.0), I)
end

function affine_transform!(Q::AbstractQuadrule{1}, I::Interval, J::Interval)
    LI = I.b - I.a
    LJ = J.b - J.a
    s = LJ / LI
    for i in 1:length(Q)
        alpha = (Q.x[i] - I.a) / LI
        Q.x[i] = (1.0 - alpha) * J.a + alpha * J.b
        Q.w[i] *= s
    end
    return Q
end

function affine_transform(Q::AbstractQuadrule{1}, I::Interval, J::Interval)
    return affine_transform!(copy(Q), I, J)
end

function QuadratureRule(Q::TensorProduct{D,<:AbstractQuadrule{1}}) where D
    x = CartesianProduct(q -> q.x, Q)
    w = KroneckerProduct(q -> q.w, Q, reverse=true)
    return QuadratureRule(x, w)
end

function QuadratureRule(method::Type{<:AbstractGaussrule}, n::Dimension, cube::Cube{D}) where D
    return QuadratureRule(TensorProduct(I -> GaussRule(method, n, I), cube.I))
end
