@testset "Inversion" begin
    @testset "State inversion" begin
        s = State(SVector{2}(1.0, 0.0), π / 3.0)
        s2 = invert(s)
        @test s2.rel_pos ≈ SVector{2}(-0.5, √3.0 / 2.0)
        @test s2.rel_rot ≈ -s.rel_rot

        s3 = invert(s2)
        @test s3.rel_pos ≈ s.rel_pos
        @test s3.rel_rot ≈ s.rel_rot
    end

    @testset "CollisionData inversion" begin
        s = State(SVector{2}(1.0, 1.0), π / 4.0)
        c = CollisionData(1.0, SVector{2}(1.0, 0.0))
        c2 = invert(c, s)
        @test c2.separation == c.separation
        @test c2.direction ≈ SVector{2}(-1.0 / √2.0, 1.0 / √2.0)
    end

    @testset "Closest Pair Inversion" begin
        c1 = Circle(1.0)
        r1 = Rect(1.0, 1.0)

        s = State(SVector{2}(2.5, 2.5), π / 4.0)
        pa, pb = closest_pair(c1, r1, s)
        @test pa ≈ SVector{2}(1.0 / √2.0, 1.0 / √2.0)
        @test pb ≈ SVector{2}(2.5 - 1.0 / √2.0, 2.5 - 1.0 / √2.0)
    end
end
