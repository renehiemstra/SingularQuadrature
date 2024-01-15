using Test

tests = [
            "base",
            "gaussquad",
            "bernstein",
            "mappings"
        ]

@testset "SingularQuadrature" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
