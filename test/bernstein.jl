using Test, SafeTestsets, Revise

@safetestset "Testing Bernstein polynomials" begin

using SingularQuadrature, CartesianProducts

@testset "Single Bernstein polynomial evaluation" begin
    I = Interval(1.0, 3.0)

    # test basis function at x=0.0
    x = I.a
    @test bernstein(I,2,0,x) ≈ 1.0
    @test bernstein(I,2,1,x) ≈ 0.0
    @test bernstein(I,2,2,x) ≈ 0.0

    # test basis function at x=1.0
    x = I.b
    @test bernstein(I,2,0,x) ≈ 0.0
    @test bernstein(I,2,1,x) ≈ 0.0
    @test bernstein(I,2,2,x) ≈ 1.0

    # test basis function at x=2.0
    x = I.a + (I.b-I.a)/2
    @test bernstein(I,2,0,x) ≈ 0.25
    @test bernstein(I,2,1,x) ≈ 0.5
    @test bernstein(I,2,2,x) ≈ 0.25
end

@testset "Basis evaluation single point" begin
    # test partition of unity different orders
    # for single point evaluation
    I = Interval(1.0, 3.0)
    X = LinRange(I.a, I.b, 4)
    for p in 2:4
        for x in X
            B = bernstein(I, 2, x)
            @test sum(B) .≈ 1
            @test all(B .>= 0.0)
        end
    end
end

@testset "Basis evaluation multiple points" begin
    # test partition of unity different orders
    # for multiple point evaluation
    I = Interval(1.0, 3.0)
    x = LinRange(I.a, I.b, 10)
    for p in 1:4
        B = bernstein(I, p, x)
        @test all(sum(B, dims=1) .≈ 1)
        @test all(B .>= 0.0)
    end
end

# @testset "Multidimensional Bernstein functions evaluated at single point" begin
#     # test partition of unity different orders
#     # for multiple point evaluation
#     I = Interval(1.0, 3.0)
#     partition = CartesianProduct(I, I, I)
#     element = first(Elements(partition))

#     B = Bernstein(element, 3, [1.0,1.0,1.0])
#     @test sum(B) ≈ 1.0
#     @test sum(B[1]) ≈ 1.0

#     B = Bernstein(element, 3, [2.0,1.0,1.0])
#     @test B[1] ≈ B[4] && B[2] ≈ B[3]
#     @test sum(B[1:4]) ≈ 1.0
#     @test sum(B) ≈ 1.0

#     B = Bernstein(element, 3, [1.5,1.0,1.0])
#     @test B[1] > B[4] && B[2] > B[3]
#     @test sum(B[1:4]) ≈ 1.0
#     @test sum(B) ≈ 1.0
# end

# @testset "Multidimensional Bernstein functions evaluated grid points" begin
#     # test partition of unity different orders
#     # for multiple point evaluation
#     I = Interval(1.0, 3.0)
#     partition = CartesianProduct(I, I, I)
#     element = first(Elements(partition))

#     x = CartesianProduct([1.5,2.0], [1.5,2.0], [1.5,2.0])

#     B = Bernstein(element, 3, x)
#     @test all(sum(B,dims=1) .≈ 1.0)
# end

end