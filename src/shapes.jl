"""
```julia
abstract type AbstractShape{D} end
```

Abstract supertype for all shapes in `D`-dimensions
"""
abstract type AbstractShape{D} end

"""
```julia
abstract type AbstractPolygon{D} <: AbstractShape{D} end
```

Abstract supertype for all polygons in `D`-dimensions. Non-polygons
can be spheres, etc.
"""
abstract type AbstractPolygon{D} <: AbstractShape{D} end

"""
```julia
struct Circle{R<:Real} <: AbstractShape{2}
```

A 2-dimensional circle, which stores its `radius`.
"""
struct Circle{R<:Real} <: AbstractShape{2}
    radius::R

    function Circle{R}(radius::R) where {R}
        @assert radius >= zero(R)
        return new{R}(radius)
    end
end

"""
```julia
function Point(R::Type{<:Real})
```

A simple wrapper over a circle with zero radius, used to represent a point.
"""
Point(R::Type{<:Real}) = Circle{R}(zero(R))

"""
```julia
struct Rect{R<:Real} <: AbstractPolygon{2}
```

A 2-dimensional rectangle which stores its half-extents (`half_ext`) in both dimensions.
"""
struct Rect{R<:Real} <: AbstractPolygon{2}
    half_ext::SVector{2,R}
end
