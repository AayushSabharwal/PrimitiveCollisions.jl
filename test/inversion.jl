@testset "Inversion" begin
    @testset "State inversion" begin
        s = State(SVector{2}(1.0, 0.0), π / 3)
        s2 = invert(s)
        @test all(s2.rel_pos .≈ SVector{2}(-0.5, √3.0 / 2.0))
        @test s2.rel_rot ≈ -s.rel_rot

        s3 = invert(s2)
        @test all(isapprox.(s3.rel_pos, s.rel_pos, atol = 1e-7))
        @test s3.rel_rot ≈ s.rel_rot
    end

    @testset "CollisionData inversion" begin
        s = State(SVector{2}(1.0, 1.0), π / 4.0)
        c = CollisionData(1.0, SVector{2}(1.0, 0.0))
        c2 = invert(c, s)
        @test c2.separation == c.separation
        @test all(c2.direction .≈ SVector{2}(-1.0 / √2.0, 1.0 / √2.0))
    end
end
