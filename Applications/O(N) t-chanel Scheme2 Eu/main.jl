using Distributed
# addprocs(2)
# nprocs()
@everywhere using FRGRealTime
@everywhere import FRGRealTime.ThresholdFunctions as TF
@everywhere using DifferentialEquations
@everywhere using HCubature
@everywhere using Dierckx
using LaTeXStrings
using Plots
using BenchmarkTools
@everywhere path2="/home/tyy/Documents/FRGRealTime.jl/Applications/O(N) t-chanel Scheme2 Eu/"
@everywhere include(joinpath(path2,"struct.jl"))
@everywhere include(joinpath(path2,"flow.jl"))

include(joinpath(path2,"solver.jl"))

@everywhere config_spec.dtmax=50
@time fourpointdataeu=tchanelSolveFourPointParallel(
    0.1,
    150.0, lpakrang,
    x -> sol1(x)[1],
    x -> sol1(x)[2],
    RelambdaflowEutest,
    6400*2,
    u0 = [6*8.0],
    config=config_spec,
)

writedlm(joinpath(path2,"Data/k=100.dat"),hcat(fourpointdataeu...))

# yaxis = :log,
plot(fourpointdataeu[1],-fourpointdataeu[2],label =L"\mathrm{Im}\lambda",title="k=200MeV",xlabel=L"p_0")
