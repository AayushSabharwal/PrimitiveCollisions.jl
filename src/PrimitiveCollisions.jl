module PrimitiveCollisions

using StaticArrays
using LinearAlgebra

export AbstractShape, AbstractPolygon, Circle, Point, Rect
include("shapes.jl")
export State, CollisionData, invert, check_collision
include("collisions.jl")

include("utils.jl")

end
