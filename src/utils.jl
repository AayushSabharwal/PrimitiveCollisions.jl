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

"""
```julia
function point_rect_projection(a::Rect, state::State{R})
```

Return the projection of point `state.rel_pos` on each edge of the [`Rect`](@ref) `a`. Order
used for edges: `(left, top, right, bottom)`.
"""
@inline function point_rect_projection(a::Rect, p::SVector{2,R}) where {R}
    return (
        SVector{2}(R(-a.half_ext[1]), clamp(p[2], -a.half_ext[2], a.half_ext[2])),
        SVector{2}(clamp(p[1], -a.half_ext[1], a.half_ext[1]), R(a.half_ext[2])),
        SVector{2}(R(a.half_ext[1]), clamp(p[2], -a.half_ext[2], a.half_ext[2])),
        SVector{2}(clamp(p[1], -a.half_ext[1], a.half_ext[1]), R(-a.half_ext[2])),
    )
end

"""
```julia
function rect_points(b::Rect, state::State{R})
```

Return the vertices of rectangle `b` with [`State`](@ref) `state`.
"""
@inline function rect_points(b::Rect, state::State{R}) where {R}
    b_center = state.rel_pos
    s, c = sincos(state.rel_rot)
    rot_mat = SMatrix{2,2}(c, s, -s, c)
    return (
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(-one(R), -one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(-one(R), one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(one(R), one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(one(R), -one(R))),
    )
end
