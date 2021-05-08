module FRGRealTime
using HCubature,
    QuadGK, DelimitedFiles, FastGaussQuadrature, Dierckx, DoubleFloats

export Epi

include("loop.jl")
include("flow.jl")


end
