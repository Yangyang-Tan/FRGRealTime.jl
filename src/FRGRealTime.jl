module FRGRealTime
using HCubature,
    QuadGK, DelimitedFiles, FastGaussQuadrature, Dierckx, DoubleFloats

export Epi
include("artifacts/flowfun.jl")

include("loop.jl")
include("flow.jl")
include("vertex.jl")


end
