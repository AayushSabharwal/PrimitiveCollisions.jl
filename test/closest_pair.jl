@testset "Closest Pair" begin
    @testset "Circle-Circle" begin
        c1 = Circle(1.0)
        c2 = Circle(2.0)

        s = State(SVector{2}(3.0, 3.0), π / 6.0)
        pa, pb = closest_pair(c1, c2, s)
        @test all(pa .≈ (1.0 / √2.0, 1.0 / √2.0))
        @test all(pb .≈ (-√2.0, -√2.0))
    end

    @testset "Circle-Rect" begin
        r1 = Rect(2.0, 1.0)
        c1 = Circle(1.0)

        s = State(SVector{2}(3.0, 3.0), π / 3.0)
        pa, pb = closest_pair(r1, c1, s)
        @test all(pa .≈ (2.0, 1.0))
        @test all(pb .≈ (-1.0 / √5.0, -2.0 / √5.0))

        s = State(SVector{2}(3.5, 0.0), 0.0)
        pa, pb = closest_pair(r1, c1, s)
        @test all(pa .≈ (2.0, 0.0))
        @test all(pb .≈ (-1.0, 0.0))
    end

    @testset "Rect-Rect" begin
        r1 = Rect(1.0, 2.0)
        r2 = Rect(0.5, 1.0)

        s = State(SVector{2}(1.5, 3.0), 0.0)
        pa, pb = closest_pair(r1, r2, s)
        @test all(pa .≈ (1.0, 2.0))
        @test all(pb .≈ (-0.5, -1.0))

        s = State(SVector{2}(50.0, 50.0), π / 4)
        pa, pb = closest_pair(r1, r2, s)
        @test all(pa .≈ (1.0, 2.0))
        @test all(pb .≈ (-0.5, 1.0 / √2.0))
    end
end
