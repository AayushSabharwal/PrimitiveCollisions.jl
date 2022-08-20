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
        closest_i = ifelseargmin(distances)
        closest_dist = if closest_i isa Integer
            distances[closest_i]
        else
            minimum(distances)
        end
        return CollisionData{R}(
            closest_dist,
            (state.rel_pos - my_getindex(points, closest_i)) / (closest_dist + b.radius),
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
    b_points = rect_points(b, state)

    x_proj = SVector{2}(minimum(x -> x[1], b_points), maximum(x -> x[1], b_points))
    y_proj = SVector{2}(minimum(x -> x[2], b_points), maximum(x -> x[2], b_points))

    x_intersect = any(-a.half_ext[1] .<= x_proj .<= a.half_ext[1])
    y_intersect = any(-a.half_ext[2] .<= y_proj .<= a.half_ext[2])

    function separation_and_axis(proj, ext)
        temp_1 = abs(proj[2] + ext)
        temp_2 = abs(ext - proj[1])
        return IfElse.ifelse(temp_1 < temp_2, (temp_1, -one(R)), (temp_2, one(R)))
    end

    x_dist, x_dir = separation_and_axis(x_proj, a.half_ext[1])
    x_dist = IfElse.ifelse(x_intersect, -x_dist, x_dist)
    x_ax = SVector{2}(x_dir, zero(R))
    y_dist, y_dir = separation_and_axis(y_proj, a.half_ext[2])
    y_dist = IfElse.ifelse(y_intersect, -y_dist, y_dist)
    y_ax = SVector{2}(zero(R), y_dir)

    return IfElse.ifelse(
        xor(x_intersect, y_intersect),
        IfElse.ifelse(x_intersect, (x_dist, x_ax), (y_dist, y_ax)),
        IfElse.ifelse(abs(x_dist) < abs(y_dist), (x_dist, x_ax), (y_dist, y_ax)),
    )
end
