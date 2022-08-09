using PrimitiveCollisions
using StaticArrays
using Test

@testset "PrimitiveCollisions.jl" begin
    include("inversion.jl")
    include("collision.jl")
    include("closest_pair.jl")
end
