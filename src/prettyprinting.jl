function Base.show(io::IO, ::MIME"text/plain", r::Rect)
    print(io, "Rect with half-extents $(r.half_ext)")
end

function Base.show(io::IO, c::Circle{F}) where {F}
    if iszero(c.radius)
        print(io, "Point($F)")
    else
        print(io, "Circle{$F}($(c.radius))")
    end
end

function Base.show(io::IO, c::Circle)
    print(io, "Circle with radius $(c.radius)")
end

function Base.show(io::IO, state::State)
    print(io, "State with relative position $(state.rel_pos) and relative rotation $(state.rel_rot) rad")
end

function Base.show(io::IO, col::CollisionData)
    print(io, "Collision with separation distance $(col.separation) in direction $(col.direction)")
end
