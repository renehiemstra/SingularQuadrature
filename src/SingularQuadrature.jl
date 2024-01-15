module SingularQuadrature

using FastGaussQuadrature, CartesianProducts, KroneckerProducts

include("base.jl")
include("mappings.jl")
include("bernstein.jl")
include("gaussquad.jl")

end # module SingularQuadrature