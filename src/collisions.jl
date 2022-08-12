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
        coldata.separation, SMatrix{2,2,R}(c, -s, s, c) * -coldata.direction
    )
end

"""
```julia
function check_collision(a::AbstractShape{2}, b::AbstractShape{2}, state::State{R}) where {R}
```

Given two shapes `a` and `b` and a `state` describing the position and rotation of `b`
relative to `a`, return the [`CollisionData`](@ref) describing their separation.
"""
function check_collision(
    a::AbstractShape{2}, b::AbstractShape{2}, state::State{R}
) where {R}
    istate = invert(state)
    return invert(check_collision(b, a, istate), istate)
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
    function center_inside_rect()
        δ = a.half_ext - abs.(state.rel_pos)
        return IfElse.ifelse(
            δ[1] <= δ[2],
            CollisionData{R}(
                -δ[1],
                SVector{2}(
                    IfElse.ifelse(state.rel_pos[1] < zero(R), -one(R), one(R)), zero(R)
                ),
            ),
            CollisionData{R}(
                -δ[2],
                SVector{2}(
                    zero(R), IfElse.ifelse(state.rel_pos[2] < zero(R), -one(R), one(R))
                ),
            ),
        )
    end

    function center_outside_rect()
        # closest point on each side of the rectangle to the circle
        points = point_rect_projection(a, state.rel_pos)

        # distances (possibly negative) from point on rectangle to border of circle
        distances = norm.(points .- (state.rel_pos,)) .- b.radius
        closest_i = argmin(distances)

        return CollisionData{R}(
            distances[closest_i],
            (state.rel_pos - points[closest_i]) / (distances[closest_i] + b.radius),
        )
    end

    return IfElse.ifelse(
        # * because all doesn't work right with SVectors of symbols:
        # https://github.com/SciML/ModelingToolkit.jl/issues/1742
        prod(abs.(state.rel_pos) .<= a.half_ext),
        center_inside_rect(),
        center_outside_rect(),
    )
end

function check_collision(a::Rect, b::Rect, state::State{R}) where {R}
    a_separation, a_axis = rect_rect_collision_util(a, b, state)
    istate = invert(state)
    b_separation, b_axis = rect_rect_collision_util(b, a, istate)

    a_coldata = CollisionData{R}(a_separation, a_axis)
    b_coldata = invert(CollisionData{R}(b_separation, b_axis), istate)

    return IfElse.ifelse(
        a_separation > zero(R) && b_separation > zero(R),
        a_separation < b_separation ? a_coldata : b_coldata,
        IfElse.ifelse(
            a_separation > zero(R),
            a_coldata,
            IfElse.ifelse(
                b_separation > zero(R),
                b_coldata,
                a_separation < b_separation ? b_coldata : a_coldata,
            ),
        ),
    )
end

# use separating axis theorem, projecting points of b onto axes of a
# only 2 axes, and we already know the projection of a.
@inbounds function rect_rect_collision_util(a::Rect, b::Rect, state::State{R}) where {R}
    b_center = state.rel_pos
    b_points = rect_points(b, state)

    # x axis
    # are the bodies separate on this axis
    separate_1 = true
    # what is the distance between them
    separation_distance_1 = R(Inf)
    # positive or negative axis
    separation_axis_1 = zero(SVector{2,R})

    function unit_svector(R::Type, index; unit = one(R))
        return IfElse.ifelse(
            index == 1, SVector{2}(unit, zero(R)), SVector{2}(zero(R), unit)
        )
    end

    function outside_left(t, separate, separation_distance, separation_axis, index)
        distance = -a.half_ext[index] - t
        return IfElse.ifelse(
            separate && distance < separation_distance,
            (separate, distance, -unit_svector(R, index)),
            (separate, separation_distance, separation_axis),
        )
    end

    function inside(t, separate, separation_distance, separation_axis, index)
        distance = min(t + a.half_ext[index], a.half_ext[index] - t)
        axis = IfElse.ifelse(
            t < zero(R),
            -unit_svector(R, index),
            IfElse.ifelse(
                t == zero(R),
                unit_svector(
                    R,
                    index;
                    unit = IfElse.ifelse(b_center[index] > zero(R), one(R), -one(R)),
                ),
                unit_svector(R, index),
            ),
        )
        return IfElse.ifelse(
            separate,
            (false, distance, axis),
            IfElse.ifelse(
                distance > separation_distance,
                (separate, distance, axis),
                (separate, separation_distance, separation_axis),
            ),
        )
    end

    function outside_right(t, separate, separation_distance, separation_axis, index)
        distance = t - a.half_ext[index]
        return IfElse.ifelse(
            separate && distance < separation_distance,
            (separate, distance, unit_svector(R, index)),
            (separate, separation_distance, separation_axis),
        )
    end

    function projection_separation(t, separate, separation_distance, separation_axis, index)
        return IfElse.ifelse(
            t < -a.half_ext[index],
            outside_left(t, separate, separation_distance, separation_axis, index),
            IfElse.ifelse(
                -a.half_ext[index] <= t <= a.half_ext[index],
                inside(t, separate, separation_distance, separation_axis, index),
                outside_right(t, separate, separation_distance, separation_axis, index),
            ),
        )
    end

    for i in 1:4
        # don't need to actually project it, just take the x coordinate
        t = b_points[i][1]
        separate_1, separation_distance_1, separation_axis_1 = projection_separation(
            t, separate_1, separation_distance_1, separation_axis_1, 1
        )
    end

    separate_2 = true
    separation_distance_2 = R(Inf)
    separation_axis_2 = zero(SVector{2,R})
    for i in 1:4
        t = b_points[i][2]
        separate_2, separation_distance_2, separation_axis_2 = projection_separation(
            t, separate_2, separation_distance_2, separation_axis_2, 2
        )
    end

    return IfElse.ifelse(
        separate_1 && separate_2,
        IfElse.ifelse(
            separation_distance_1 < separation_distance_2,
            (separation_distance_1, separation_axis_1),
            (separation_distance_2, separation_axis_2),
        ),
        IfElse.ifelse(
            separate_1,
            (separation_distance_1, separation_axis_1),
            IfElse.ifelse(
                separate_2,
                (separation_distance_2, separation_axis_2),
                IfElse.ifelse(
                    separation_distance_1 < separation_distance_2,
                    (-separation_distance_1, separation_axis_1),
                    (-separation_distance_2, separation_axis_2),
                ),
            ),
        ),
    )
end
