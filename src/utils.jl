@inline sqnorm(pt::SVector) = dot(pt, pt)

@inline function point_line_segment_projection(
    p::SVector{2,F}, u::SVector{2,F}, v::SVector{2,F}
) where {F}
    t = clamp(point_line_projection_parameter(p, u, v), zero(F), one(F))
    return u + t * (v - u)
end

@inline function point_line_projection(
    p::SVector{2,F}, u::SVector{2,F}, v::SVector{2,F}
) where {F}
    return u + point_line_projection_parameter(p, u, v) * (v - u)
end

@inline function point_line_projection_parameter(
    p::SVector{2,F}, u::SVector{2,F}, v::SVector{2,F}
) where {F}
    return dot(p - u, v - u) / sqnorm(v - u)
end
