module PrimitiveCollisions

using StaticArrays
using LinearAlgebra
using IfElse

export AbstractShape, AbstractPolygon, Circle, Point, Rect
include("shapes.jl")
export State, CollisionData, invert, check_collision
include("collisions.jl")
export closest_pair
include("closest_pair.jl")

include("utils.jl")
include("prettyprinting.jl")

import PrecompileTools

PrecompileTools.@compile_workload begin
    c1 = Circle(1.0)
    r1 = Rect(1.0, 1.0)

    s = State(SVector{2}(2.5, 2.5), π / 4.0)
    closest_pair(c1, r1, s)

    r1 = Rect(2.0, 1.0)
    c1 = Circle(2.0)

    s = State(SVector{2}(0.5, 0.5), π / 3.0)
    col = check_collision(c1, r1, s)
end

end
