using WENO4
using Test

@testset "WENO4" begin
    x = 0:0.2:20
    y = sin.(x)
    @test interpolate_weno4(x, x, y) ≈ y  atol=1e-10
    f(x) = 2.0 + 3.5 * x - 6.9 * x^2
    y = f.(x)
    @test interpolate_weno4(x, x, y) ≈ y  atol=1e-10
    # Constant edges
    @test interpolate_weno4([-1.0], x, y)[1] == y[1]
    @test interpolate_weno4([21.0], x, y)[1] == y[end]
    # Extrapolated edges, exact solution
    @test interpolate_weno4([-10.0], x, y, extrapolate=true)[1] ≈ f(-10) atol=1e-10
    @test interpolate_weno4([50.0], x, y, extrapolate=true)[1] ≈ f(50.) atol=1e-10
    # Zig-zag profile with enough points to fit polynomial, interpolation should
    # never be beyond extremes
    x = [0., 0.5, 0.75, 1., 1+1e-10, 1.25, 1.5, 2]
    y = [0, 0.2, 0.75, 1., -1, -0.75, -0.2, 0]
    new_x = 0.8:0.01:1.2
    new_y = interpolate_weno4(new_x, x, y)
    @test maximum(new_y) <= maximum(y)
    @test minimum(new_y) >= minimum(y)
    @test WENO4.binary_search(a, 1.25) == findfirst(x -> x>1.25, a) - 1
    @test WENO4.weno4_q(1, 1, 2, 3, 4, 1, 2.5, 3.5, 4.5) == (1.0, 1.5)
    @test WENO4.weno4_q(1, 1, 2, 3, 4, 0, 0, 0, 0) == (0.0, 0.0)
    @test WENO4.weno4_β(1, 2, 3, 4, 0, 0, 0, 0.5) == (0.0, 1.0)
    @test WENO4.weno4_β(1, 2, 3, 4, 0, 0, 0, 0) == (0.0, 0.0)
end
