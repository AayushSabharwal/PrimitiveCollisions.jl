"""
```julia
function closest_pair(a::AbstractShape{2}, b::AbstractShape{2}, state::State{R})
```

Given two shapes `a` and `b` and a `state` describing the position and rotation of `b`
relative to `a`, return the points on the boundary of `a` and `b` respectively which
are closest to each other. Only valid for shapes that do not intersect.
"""
function closest_pair(a::AbstractShape{2}, b::AbstractShape{2}, state::State{R}) where {R}
    istate = invert(state)
    return reverse(closest_pair(b, a, istate))
end

function closest_pair(a::Circle, b::Circle, state::State{R}) where {R}
    dist = norm(state.rel_pos)
    return (state.rel_pos / dist * a.radius, -state.rel_pos / dist * b.radius)
end

function closest_pair(a::Rect, b::Circle, state::State{R}) where {R}
    points = point_rect_projection(a, state.rel_pos)
    sq_distances = sqnorm.(points .- (state.rel_pos,))
    closest_i = argmin(sq_distances)
    δ = (points[closest_i] .- state.rel_pos)
    return (points[closest_i], δ / norm(δ) * b.radius)
end

function closest_pair(a::Rect, b::Rect, state::State{R}) where {R}
    pa, pb = rect_closest_pair_util(a, b, state)
    istate = invert(state)
    qb, qa = rect_closest_pair_util(b, a, istate)

    if sqnorm(pa - pb) < sqnorm(qa - qb)
        s, c = sincos(-state.rel_rot)
        rot_mat = SMatrix{2,2}(c, s, -s, c)
        return pa, rot_mat * (pb - state.rel_pos)
    else
        s, c = sincos(-istate.rel_rot)
        rot_mat = SMatrix{2,2}(c, s, -s, c)
        return rot_mat * (qa - istate.rel_pos), qb
    end
end

function rect_closest_pair_util(a::Rect, b::Rect, state::State{R}) where {R}
    b_points = rect_points(b, state)

    projections = point_rect_projection.((a,), b_points)
    distances = Tuple(sqnorm.(projections[i] .- (b_points[i],)) for i in 1:4)
    closest_is = argmin.(distances)
    closest_closest_i = argmin(getindex.(distances, closest_is))
    best_index = closest_is[closest_closest_i]

    return projections[closest_closest_i][best_index], b_points[closest_closest_i]
end