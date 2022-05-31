module FRGRealTime
using HCubature,
    QuadGK,
    DelimitedFiles,
    FastGaussQuadrature,
    Dierckx,
    DoubleFloats,
    Roots,
    DifferentialEquations,
    FiniteDifferences

export Eb, RGscale, ODEConfig, SpecSolution, LambdaGridini, LambdaGrid,reverse2!,Interpolate_Linear,Interpolate_Log_Linear
include("artifacts/math.jl")

using .MathFuns

include("artifacts/flowfun.jl")
include("artifacts/integrate.jl")
include("artifacts/struct.jl")

include("GPU/GPUflowfun.jl")
include("GPU/GPUflow.jl")

include("loop.jl")
include("flow.jl")
include("ThresholdFunction.jl")
include("Im.jl")
include("Re.jl")
include("vertexIm.jl")
include("vertexImtest.jl")
include("vertexImdq0.jl")
include("vertexImdq0test.jl")
include("vertexRe.jl")
include("vertexzero.jl")
include("zerosolver.jl")

end
