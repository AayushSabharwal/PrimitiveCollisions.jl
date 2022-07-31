"""
```julia
struct State{R<:Real}
```

Represents the current state of spatial configuration between two shapes (`A` and `B`) in
2-dimensional space. One of the shapes (`A`) is used as the frame of reference, thus residing
at the origin with no rotation.

# Fields
- `rel_pos::SVector{2,R}`: The position of `B` with respect to `A`.
- `rel_rot::R`: The rotation of `B` with respect to `A`. Rotations are counter-clockwise
  measured from the positive x-axis, in radians.
"""
struct State{R<:Real}
    rel_pos::SVector{2,R}
    rel_rot::R
end

State(rel_pos::SVector{2,R}, rel_rot::Irrational) where {R} = State(rel_pos, R(rel_rot))

"""
```julia
function invert(state::State{R}) where {R}
```

Given a `State` containing the configuration of `B` with respect to `A`, returns the state
containing the configuration of `A` with respect to `B`.
"""
function invert(state::State{R}) where {R}
    s, c = sincos(state.rel_rot)
    return State{R}(SMatrix{2,2,R}(c, -s, s, c) * -state.rel_pos, -state.rel_rot)
end

"""
```julia
struct CollisionData{R<:Real}
```

Struct which stores the collision information about two bodies in 2-dimensional space.

# Fields
- `separation::R`: If the bodies do not intersect, this is the shortest distance between any
  two points on their boundaries. If the bodies intersect, this is the negative of the minimum
  displacement magnitude required to separate them.
- `direction::SVector{2,R}`: Unit vector. If the bodies do not intersect, this is the
  direction of the minimum distance vector from body `A` to body `B`. If the bodies intersect,
  this is the direction to displace body `B` to separate them.
"""
struct CollisionData{R<:Real}
    separation::R
    direction::SVector{2,R}
end

"""
```julia
function invert(coldata::CollisionData{R}, current_state::State{R}) where {R}
```

Given collision data in the `current_state` frame of reference, transform it
into the opposite frame of reference.
"""
function invert(coldata::CollisionData{R}, current_state::State{R}) where {R}
    s, c = sincos(current_state.rel_rot)
    return CollisionData{R}(
        coldata.separation,
        SMatrix{2,2,R}(c, -s, s, c) * -coldata.direction,
    )
end

"""
```julia
function check_collision(a::AbstractShape{2}, b::AbstractShape{2}, state::State{R}) where {R}
```

Given two shapes `a` and `b` and a `state` describing the position and rotation of `b`
relative to `a`, return the [`CollisionData`](@ref) describing their separation.
"""
function check_collision(a::AbstractShape{2}, b::AbstractShape{2}, state::State{R}) where {R}
    istate = invert(state)
    invert(check_collision(b, a, istate), istate)
end

function check_collision(a::Circle, b::Circle, state::State{R}) where {R}
    # circles are easy
    center_dist = norm(state.rel_pos)
    separation = center_dist - a.radius - b.radius
    return CollisionData{R}(separation, state.rel_pos / center_dist)
end

# making rect rotation-invariant simplifies things, since circle already doesn't care
@inbounds function check_collision(a::Rect, b::Circle, state::State{R}) where {R}
    # check if center of circle inside rectangle
    if all(abs.(state.rel_pos) .<= a.half_ext)
        δ = a.half_ext - abs.(state.rel_pos)
        if δ[1] <= δ[2]
            return CollisionData{R}(
                -δ[1],
                state.rel_pos[1] < zero(R) ? SVector{2}(-one(R), zero(R)) :
                SVector{2}(one(R), zero(R)),
            )
        else
            return CollisionData{R}(
                -δ[2],
                state.rel_pos[2] < zero(R) ? SVector{2}(zero(R), -one(R)) :
                SVector{2}(zero(R), one(R)),
            )
        end
    end

    # closest point on each side of the rectangle to the circle
    points = (
        point_line_segment_projection(
            state.rel_pos,
            SVector{2}(R(-a.half_ext[1]), R(-a.half_ext[2])),
            SVector{2}(R(-a.half_ext[1]), R(a.half_ext[2])),
        ),
        point_line_segment_projection(
            state.rel_pos,
            SVector{2}(R(-a.half_ext[1]), R(a.half_ext[2])),
            SVector{2}(R(a.half_ext[1]), R(a.half_ext[2])),
        ),
        point_line_segment_projection(
            state.rel_pos,
            SVector{2}(R(a.half_ext[1]), R(a.half_ext[2])),
            SVector{2}(R(a.half_ext[1]), R(-a.half_ext[2])),
        ),
        point_line_segment_projection(
            state.rel_pos,
            SVector{2}(R(a.half_ext[1]), R(-a.half_ext[2])),
            SVector{2}(R(-a.half_ext[1]), R(-a.half_ext[2])),
        ),
    )

    # distances (possibly negative) from point on rectangle to border of circle
    distances = norm.(points .- (state.rel_pos,)) .- b.radius
    closest_i = argmin(distances)

    return CollisionData{R}(
        distances[closest_i],
        (state.rel_pos - points[closest_i]) * b.radius / (distances[closest_i] + b.radius),
    )

end

function check_collision(a::Rect, b::Rect, state::State{R}) where {R}
    a_separation, a_axis = rect_rect_collision_util(a, b, state)
    istate = invert(state)
    b_separation, b_axis = rect_rect_collision_util(b, a, istate)

    a_coldata = CollisionData{R}(a_separation, a_axis)
    b_coldata = invert(CollisionData{R}(b_separation, b_axis), istate)

    if a_separation > zero(R) && b_separation > zero(R)
        return a_separation < b_separation ? a_coldata : b_coldata
    elseif a_separation > zero(R)
        return a_coldata
    elseif b_separation > zero(R)
        return b_coldata
    else
        return a_separation < b_separation ? b_coldata : a_coldata
    end
end

# use separating axis theorem, projecting points of b onto axes of a
# only 2 axes, and we already know the projection of a.
@inbounds function rect_rect_collision_util(a::Rect, b::Rect, state::State{R}) where {R}
    b_center = state.rel_pos

    s, c = sincos(state.rel_rot)
    rot_mat = SMatrix{2,2}(c, s, -s, c)
    b_points = (
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(-one(R), -one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(-one(R), one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(one(R), one(R))),
        b_center + rot_mat * (R.(b.half_ext) .* SVector{2}(one(R), -one(R))),
    )

    # x axis
    # are the bodies separate on this axis
    separate_1 = true
    # what is the distance between them
    separation_distance_1 = R(Inf)
    # positive or negative axis
    separation_axis_1 = zero(SVector{2,R})
    for i = 1:4
        # don't need to actually project it, just take the x coordinate
        t = b_points[i][1]
        # outside
        if t < -a.half_ext[1]
            distance = -a.half_ext[1] - t
            if separate_1 && distance < separation_distance_1
                separation_distance_1 = distance
                separation_axis_1 = SVector{2,R}(-one(R), zero(R))
            end
        # inside
        elseif -a.half_ext[1] <= t <= a.half_ext[1]
            distance = min(t + a.half_ext[1], a.half_ext[1] - t)
            axis = if t < zero(R)
                SVector{2}(-one(R), zero(R))
            elseif t == zero(R)
                if b_center[1] > zero(R)
                    SVector{2}(one(R), zero(R))
                else
                    SVector{2}(-one(R), zero(R))
                end
            else
                SVector{2}(one(R), zero(R))
            end
            if separate_1
                separate_1 = false
                separation_distance_1 = distance
                separation_axis_1 = axis
            elseif distance > separation_distance_1
                separation_distance_1 = distance
                separation_axis_1 = axis
            end
        else
            # outside
            distance = t - a.half_ext[1]
            if separate_1 && distance < separation_distance_1
                separation_distance_1 = distance
                separation_axis_1 = SVector{2,R}(one(R), zero(R))
            end
        end
    end

    separate_2 = true
    separation_distance_2 = R(Inf)
    separation_axis_2 = zero(SVector{2,R})
    for i = 1:4
        t = b_points[i][2]
        if t < -a.half_ext[2]
            distance = -a.half_ext[2] - t
            if separate_2 && distance < separation_distance_2
                separation_distance_2 = distance
                separation_axis_2 = SVector{2}(zero(R), -one(R))
            end
        elseif -a.half_ext[2] <= t <= a.half_ext[2]
            distance = min(t + a.half_ext[2], a.half_ext[2] - t)
            axis = if t < zero(R)
                SVector{2}(zero(R), -one(R))
            elseif t == zero(R)
                if b_center[2] > zero(R)
                    SVector{2}(zero(R), -one(R))
                else
                    SVector{2}(zero(R), one(R))
                end
            else
                SVector{2}(zero(R), one(R))
            end
            distance = t + a.half_ext[2]
            if separate_2
                separate_2 = false
                separation_distance_2 = distance
                separation_axis_2 = axis
            elseif distance > separation_distance_2
                separation_distance_2 = distance
                separation_axis_2 = axis
            end
        else
            distance = t - a.half_ext[2]
            if separate_2 && distance < separation_distance_2
                separation_distance_2 = distance
                separation_axis_2 = SVector{2}(zero(R), one(R))
            end
        end
    end

    if separate_1 && separate_2
        if separation_distance_1 < separation_distance_2
            return separation_distance_1, separation_axis_1
        else
            return separation_distance_2, separation_axis_2
        end
    elseif separate_1
        return separation_distance_1, separation_axis_1
    elseif separate_2
        return separation_distance_2, separation_axis_2
    else
        if separation_distance_1 < separation_distance_2
            return -separation_distance_1, separation_axis_1
        else
            return -separation_distance_2, separation_axis_2
        end
    end
end
