#load package
using Revise
using FRGRealTime
import FRGRealTime.ThresholdFunctions as TF
using CUDA
device!(1)
using DifferentialEquations
using HCubature
using Dierckx, Interpolations, LinearAlgebra
using Plots
# using BenchmarkTools, SharedArrays, LoopVectorization, ThreadsX, Referenceables
# using Tullio

include(joinpath(path, "flowGPU.jl"))
config_lpa = ODEConfig(dtmax = 100.0)
config_spec = ODEConfig(
    dtmax = 50.0,
    progress = false,
    dense = false,
    save_on = false,
    save_start = false,
    save_end = true,
    atol = 1e-5,
    rtol = 1e-5,
)



ϵ5 = 0.5
ϵ5GPU = 5.0f0
lpakrang = RGscale(800.0f0, 1.0f0)
lpakrangGPU = RGscale(800.0f0, 1.0f0)
ImlambdaUV = 1e-7
ImlambdaUVGPU = 1f-7
RelambdaUV = 3 * 8.0
RelambdaUVGPU = 3 * 8.0f0
RelambdaAUVGPU = 3 * 8.0f0
RelambdaBUVGPU = 1.0f0


function init_fourpoint_2d(p0grid::AbstractArray{T}, q0grid::AbstractArray{T}, m2) where {T}
    N1 = length(p0grid)
    N2 = length(q0grid)
    u = zeros(T, N1, N2 + 1, 2)
    for I in CartesianIndices((N1, N2 + 1))
        u[I, 1] = 0.0
        # u[I,1] = -1f-7
        u[I, 2] = RelambdaAUVGPU
    end

    for I in CartesianIndices((N1, 1))
        u[I, 1] = 0.0
        u[I, 2] = p0grid[I[1]]^2 - m2
    end
    #   for I in CartesianIndices((N1, margin))
    #     u[I,1] = 0.0f0
    #     u[I,2] = 0.0f0
    #   end
    return u
end


function init_fourpoint_2d_InvRe(p0grid::AbstractArray{T}, q0grid::AbstractArray{T}, m2) where {T}
    N1 = length(p0grid)
    N2 = length(q0grid)
    u = zeros(T, N1, N2 + 1, 2)
    for I in CartesianIndices((N1, N2 + 1))
        u[I, 1] = 0.0
        # u[I,1] = -1f-7
        u[I, 2] = 1/RelambdaAUVGPU
    end

    for I in CartesianIndices((N1, 1))
        u[I, 1] = 0.0
        u[I, 2] = p0grid[I[1]]^2 - m2
    end
    #   for I in CartesianIndices((N1, margin))
    #     u[I,1] = 0.0f0
    #     u[I,2] = 0.0f0
    #   end
    return u
end



function init_fourpoint_2d_AB(p0grid::AbstractArray{T}, q0grid::AbstractArray{T}, m2) where {T}
    N1 = length(p0grid)
    N2 = length(q0grid)
    u = zeros(T, N1, N2 + 1, 3)
    for I in CartesianIndices((N1, N2 + 1))
        u[I, 1] = 0.0
        # u[I,1] = -1f-7
        u[I, 2] = 0.0
        u[I, 3] = RelambdaUV
    end

    for I in CartesianIndices((N1, 1))
        u[I, 1] = 0.0
        u[I, 2] = p0grid[I[1]]^2 - m2
    end
    #   for I in CartesianIndices((N1, margin))
    #     u[I,1] = 0.0f0
    #     u[I,2] = 0.0f0
    #   end
    return u
end
