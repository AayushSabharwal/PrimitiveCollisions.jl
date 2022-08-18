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
    closest_i = ifelseargmin(sq_distances)
    closest_point = my_getindex(points, closest_i)
    δ = (closest_point .- state.rel_pos)
    return (closest_point, δ / norm(δ) * b.radius)
end

function closest_pair(a::Rect, b::Rect, state::State{R}) where {R}
    a_points = rect_points(a, State(zero(SVector{2,R}), zero(R)))
    b_points = rect_points(b, state)
    function projection_util(x, p, q)
        y = point_line_segment_projection(x, p, q)
        return (sqnorm(y - x), x, y)
    end
    dist_1 = Tuple(projection_util(a_points[i], b_points[j], b_points[mod1(j+1, 4)]) for i in 1:4, j in 1:4)
    dist_2 = Tuple(projection_util(b_points[i], a_points[j], a_points[mod1(j+1, 4)]) for i in 1:4, j in 1:4)
    closest_1 = argmin(x -> x[1], dist_1)
    closest_2 = argmin(x -> x[1], dist_2)
    @show closest_1 closest_2
    if closest_1[1] < closest_2[1]
        s, c = sincos(state.rel_rot)
        return closest_1[2], SMatrix{2,2}(c, -s, s, c) * (closest_1[3] - state.rel_pos)
    else
        s, c = sincos(state.rel_rot)
        return closest_2[3], SMatrix{2,2}(c, -s, s, c) * (closest_2[2] - state.rel_pos)
    end
end

