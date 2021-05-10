using FRGRealTime
using Test
using FiniteDifferences
@testset "flowpp" begin
    @test FRGRealTime.flowpp(0.1, 1.0, 5.0, -2.0, 10.0)==central_fdm(5, 1)(k->FRGRealTime.loopfunpp(0.1, 1.0, k, -2.0, 10.0),5.0)
    @test FRGRealTime.flowpp(0.5, 1.0, 5.0, -2.0, 10.0)==central_fdm(5, 1)(k->FRGRealTime.loopfunpp(0.5, 1.0, k, -2.0, 10.0),5.0)
    @test FRGRealTime.flowpp(1.1, 1.0, 5.0, -2.0, 10.0)==central_fdm(5, 1)(k->FRGRealTime.loopfunpp(1.1, 1.0, k, -2.0, 10.0),5.0)
    @test FRGRealTime.flowpp(10.0, 1.0, 5.0, -2.0, 10.0)≈central_fdm(5, 1)(k->FRGRealTime.loopfunpp(10.0, 1.0, k, -2.0, 10.0),5.0)
end

FRGRealTime.loopfunpm(1e-12,1e-11,1.0,0.5,2.0)
using Plots


plot(x->FRGRealTime.loopfunpm(x,1e-10,1.0,0.5,2.0),1e-12,1e-10)
