module PrimitiveCollisions

using StaticArrays
using LinearAlgebra

export AbstractShape, AbstractPolygon, Circle, Point, Rect
include("shapes.jl")
export State, CollisionData, invert, check_collision
include("collisions.jl")
export closest_pair
include("closest_pair.jl")

include("utils.jl")

end
