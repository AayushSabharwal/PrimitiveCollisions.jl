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
    pb, pa = closest_pair(b, a, istate)
    s, c = sincos(istate.rel_rot)
    rot_mat = SMatrix{2,2}(c, -s, s, c)
    return rot_mat * (pa - istate.rel_pos), rot_mat * (pb - istate.rel_pos)
end

function closest_pair(a::Circle, b::Circle, state::State{R}) where {R}
    dist = norm(state.rel_pos)
    return state.rel_pos / dist * a.radius, state.rel_pos / dist * (dist - b.radius)
end

function closest_pair(a::Rect, b::Circle, state::State{R}) where {R}
    points = point_rect_projection(a, state.rel_pos)
    sq_distances = sqnorm.(points .- (state.rel_pos,))
    closest_i = argmin(sq_distances)
    closest_point = points[closest_i]
    δ = (state.rel_pos .- closest_point)
    return closest_point,
    closest_point + δ / norm(δ) * (sqrt(sq_distances[closest_i]) - b.radius)
end

function closest_pair(a::Rect, b::Rect, state::State{R}) where {R}
    if isapprox(floor(2 * state.rel_rot / π), 2 * state.rel_rot / π)
        b_hext = if isapprox(floor(state.rel_rot / π), state.rel_rot / π)
            b.half_ext
        else
            reverse(b.half_ext)
        end
        bx_ext = SVector{2}(-b_hext[1], b_hext[1]) .+ state.rel_pos[1]
        by_ext = SVector{2}(-b_hext[2], b_hext[2]) .+ state.rel_pos[2]

        x_int = any(-a.half_ext[1] .<= bx_ext .<= a.half_ext[1])
        y_int = any(-a.half_ext[2] .<= by_ext .<= a.half_ext[2])

        x_common = SVector{2}(max(bx_ext[1], -a.half_ext[1]), min(bx_ext[2], a.half_ext[1]))
        y_common = SVector{2}(max(by_ext[1], -a.half_ext[2]), min(by_ext[2], a.half_ext[2]))

        x_len = x_common[2] - x_common[1]
        y_len = y_common[2] - y_common[1]

        if x_int || y_int
            if (x_int && !y_int) || (x_int && y_int && x_len >= y_len)
                midpoint = (x_common[1] + x_common[2]) / 2
                return (
                    SVector{2}(
                        midpoint,
                        state.rel_pos[2] > zero(R) ? a.half_ext[2] : -a.half_ext[2],
                    ),
                    SVector{2}(
                        midpoint,
                        state.rel_pos[2] -
                        (state.rel_pos[2] > zero(R) ? b_hext[2] : -b_hext[2]),
                    ),
                )
            else
                midpoint = (y_common[1] + y_common[2]) / 2
                return (
                    SVector{2}(
                        state.rel_pos[1] > zero(R) ? a.half_ext[1] : -a.half_ext[1],
                        midpoint,
                    ),
                    SVector{2}(
                        state.rel_pos[1] -
                        (state.rel_pos[1] > zero(R) ? b_hext[1] : (-b_hext[1])),
                        midpoint,
                    ),
                )
            end
        end
    end
    a_points = rect_points(a, State(zero(SVector{2,R}), zero(R)))
    b_points = rect_points(b, state)
    function projection_util(x, p, q)
        y = point_line_segment_projection(x, p, q)
        return (sqnorm(y - x), x, y)
    end
    dist_1 = Tuple(
        projection_util(a_points[i], b_points[j], b_points[mod1(j + 1, 4)]) for i in 1:4,
        j in 1:4
    )
    dist_2 = Tuple(
        projection_util(b_points[i], a_points[j], a_points[mod1(j + 1, 4)]) for i in 1:4,
        j in 1:4
    )
    closest_1 = argmin(x -> x[1], dist_1)
    closest_2 = argmin(x -> x[1], dist_2)
    if closest_1[1] < closest_2[1]
        return closest_1[2], closest_1[3]
    else
        return closest_2[3], closest_2[2]
    end
end
