using Test, SafeTestsets, Revise

@safetestset "Elementary Shapes" begin

using SingularQuadrature, CartesianProducts

I = Interval(0.0,1.0)

@testset "Cubes" begin
    @test Cube(I ⨷ I ⨷ I) isa Cube{3}
    @test Cube(I ⨷ I ⨷ I ⨷ I) isa Cube{4}
    @test pullback(Cube(I ⨷ I ⨷ I)) isa Cube{3}
end

@testset "Pyramids" begin
    @test Pyramid(I ⨷ I, 0.0) isa Pyramid{3}
    @test Pyramid(I ⨷ I ⨷ I, 0.0) isa Pyramid{4}
    @test pullback(Pyramid(I ⨷ I, 0.0)) isa Cube{3}
end


end