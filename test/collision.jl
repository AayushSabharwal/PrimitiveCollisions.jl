@testset "Collision" begin
    @testset "Circle-Circle" begin
        c1 = Circle(1.0)
        c2 = Circle(2.0)

        s = State(SVector{2}(1.0, 1.0), π / 6.0)
        col = check_collision(c1, c2, s)
        @test col.separation ≈ √2.0 - 3.0
        @test all(col.direction .≈ 1.0 / √2.0)

        s = State(SVector{2}(3.0, 3.0), π / 4.0)
        col = check_collision(c1, c2, s)
        @test col.separation ≈ 3.0√2.0 - 3.0
        @test all(col.direction .≈ 1.0 / √2.0)
    end

    @testset "Circle-Rect" begin
        r1 = Rect(2.0, 1.0)
        c1 = Circle(1.0)

        s = State(SVector{2}(0.5, 0.5), π / 3.0)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ -0.5
        @test all(col.direction .≈ SVector{2}(0.0, 1.0))

        s = State(SVector{2}(1.9, 0.0), 0.0)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ -0.1
        @test all(col.direction .≈ SVector{2}(1.0, 0.0))

        s = State(SVector{2}(-0.5, -0.5), π / 2.0)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ -0.5
        @test all(col.direction .≈ SVector{2}(0.0, -1.0))

        s = State(SVector{2}(-1.9, 0.0), π / 1.2)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ -0.1
        @test all(col.direction .≈ SVector{2}(-1.0, 0.0))

        s = State(SVector{2}(1.5, 1.5), π)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ -0.5
        @test all(col.direction .≈ SVector{2}(0.0, 1.0))

        s = State(SVector{2}(2.2, 1.2), π + 1.3)
        col = check_collision(r1, c1, s)
        @test col.separation ≈ 0.2√2.0 - 1.0
        @test all(col.direction .≈ 1.0 / √2.0)
    end

    @testset "Rect-Rect" begin
        r1 = Rect(1.0, 2.0)
        r2 = Rect(0.5, 1.0)

        s = State(SVector{2}(0.0, 0.0), 0.0)
        col = check_collision(r1, r2, s)
        @test col.separation ≈ 0.5
        @test all(abs.(col.direction) .≈ SVector{2}(1.0, 0.0))

        s = State(SVector{2}(1.0, 2.0), 0.0)
        col = check_collision(r1, r2, s)
        @test col.separation ≈ -0.5
        @test all(col.direction .≈ SVector{2}(1.0, 0.0))

        s = State(SVector{2}(1.0, 2.0), π / 4)
        col = check_collision(r1, r2, s)
        @test col.separation ≈ -0.5
        @test all(col.direction .≈ 1.0 / √2.0)
    end
end
