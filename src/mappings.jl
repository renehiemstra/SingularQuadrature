export Polyhedron, Cube, Pyramid, pullback

abstract type Polyhedron{D} end

Base.ndims(::Polyhedron{D}) where D = D

struct Cube{D,T<:Real} <: Polyhedron{D}
    I::TensorProduct{D, Interval{T}}
end

struct Pyramid{D, X<:TensorProduct, T<:Real} <: Polyhedron{D}
    base::X
    z::T
    function Pyramid(base::X, z::T) where {X,T}
        M = length(base)
        @assert base isa TensorProduct{M,<:Interval{T}} "Base should be of type TensorProduct{D, <:Interval{T}}"
        return new{M+1,X,T}(base, z)
    end
end

pullback(b::Polyhedron) = Cube(TensorProduct(k -> Interval(0.0,1.0), 1:ndims(b)))
