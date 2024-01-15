using Test, SafeTestsets, Revise

@safetestset "Gauss quadrature" begin

using LinearAlgebra, SingularQuadrature, CartesianProducts, KroneckerProducts

@testset "Gauss rules" begin

    @test_throws AssertionError GaussRule{Legendre}([-1.0,0.0,1.0],[1.0,1.0])

    @test length(GaussRule(3))==3
    @test length(GaussRule(Legendre, 3))==3
    @test length(GaussRule(Lobatto, 3))==3

    Q = GaussRule(Lobatto, 3, Interval(0.0, 3.0))
    @test Q.x ≈ [0.0,1.5,3.0]
    @test sum(Q.w) ≈ 3.0

    # A Gauss Legendre rule of n points can exactly integrate polynomials
    # up to degree 2n-1
    Q = GaussRule(Legendre, 3, Interval(0.0, 2.0))
    f(x,p) = x.^p
    for k in 0:5
        @test dot(Q.w, f(Q.x, k)) ≈ 2.0^(k+1) / (k+1)
    end
    @test !(dot(Q.w, f(Q.x, 6)) ≈ 2.0^(7) / 7)

    # A Gauss Lobatto rule of n points can exactly integrate polynomials
    # up to degree 2n-3
    Q = GaussRule(Lobatto, 3, Interval(0.0, 2.0))
    f(x,p) = x.^p
    for k in 0:3
        @test dot(Q.w, f(Q.x, k)) ≈ 2.0^(k+1) / (k+1)
    end
    @test !(dot(Q.w, f(Q.x, 4)) ≈ 2.0^(5) / 5)

    # quadrature rules can be affinely mapped to a new interval
    Q = SingularQuadrature.affine_transform!(GaussRule(Lobatto, 3), Interval(-1.0,1.0),  Interval(0.0,3.0))
    @test Q.x ≈ [0.0,1.5,3.0]
    @test sum(Q.w) ≈ 3.0
end

@testset "Tensor product Gauss rules" begin
    Ω = Interval(0.0, 1.0) ⨱ Interval(0.0, 1.0)
    Q = QuadratureRule(TensorProduct(I -> GaussRule(Legendre, 3, I), Ω))
    @test sum(Q.w) ≈ 1.0
end

@testset "Integration over hyper cubes" begin
    I = Interval(0.0,1.0)
    Q = QuadratureRule(Legendre, 3, Cube(I ⨷ I ⨷ I ⨷ I ⨷ I))
    @test sum(Q.w) ≈ 1.0
end

end # safetestset