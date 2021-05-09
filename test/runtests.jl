using FRGRealTime
using Test
using FiniteDifferences
@testset "flowpp" begin
    @test FRGRealTime.flowpp(0.1, 1.0, 5.0, -2.0, 10.0)==central_fdm(5, 1)(k->FRGRealTime.loopfunpp(0.1, 1.0, k, -2.0, 10.0),5.0)
    # Write your tests here.
end