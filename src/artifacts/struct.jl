using CUDA
struct RGscale{T<:Number}
    IR::T
    UV::T
    kspan::Tuple{T,T}
end
RGscale(UV, IR) = RGscale(IR, UV, (UV, IR))



mutable struct ODEConfig
    atol::AbstractFloat
    rtol::AbstractFloat
    dtmax::AbstractFloat
    dtmin::AbstractFloat
    maxiters::Int
    progress::Bool
    adaptive::Bool
    dense::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm
end


function ODEConfig(;
    atol = 1e-14,
    rtol = 1e-14,
    dtmax = 0.1,
    dtmin = 1e-14,
    maxiters = 10^5,
    progress = true,
    adaptive = true,
    dense = true,
    save_on = true,
    save_start = true,
    save_end = true,
    alg = Tsit5(),
)
    ODEConfig(
        atol,
        rtol,
        dtmax,
        dtmin,
        maxiters,
        progress,
        adaptive,
        dense,
        save_on,
        save_start,
        save_end,
        alg,
    )
end


struct SpecSolution
    sol::Any
    u0::AbstractArray
    ImΓ2ini::AbstractArray
    ReΓ2ini::AbstractArray
    k::AbstractArray
    m2::AbstractArray
    λ0::AbstractArray
    p0::AbstractArray
    q0::AbstractArray
    ImΓ2::AbstractArray
    ReΓ2::AbstractArray
    spec::AbstractArray
end


mutable struct LambdaGrid
    Imlambdap0q0::CuArray
    Imlambdap0mq0::CuArray
    Imlambdap0mp0::CuArray
    Imlambdaq0mq0::CuArray
    Imlambdap0Ek::CuArray
    Imlambdap0mEk::CuArray
    Imlambdaq0Ek::CuArray
    Imlambdaq0mEk::CuArray
    Relambdap0q0::CuArray
    Relambdap0p0::CuArray
    Relambdaq0q0::CuArray
    Relambdap0Ek::CuArray
    Relambdaq0Ek::CuArray
end

# function LambdaGridini(T::DataType, dim::Int)
#     LambdaGrid(
#         CuArray{T}(undef, dim, dim),
#         CuArray{T}(undef, dim, dim),
#         CuArray{T}(undef, dim),
#         CuArray{T}(undef, 1, dim),
#         CuArray{T}(undef, dim),
#         CuArray{T}(undef, dim),
#         CuArray{T}(undef, 1, dim),
#         CuArray{T}(undef, 1, dim),
#         CuArray{T}(undef, dim, dim),
#         CuArray{T}(undef, dim),
#         CuArray{T}(undef, 1, dim),
#         CuArray{T}(undef, dim),
#         CuArray{T}(undef, 1, dim),
#     )
# end




function LambdaGridini(T::DataType, dim::Int)
    LambdaGrid(
        CuArray{T}(undef, dim, dim),
        CuArray{T}(undef, dim, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
        CuArray{T}(undef, dim),
    )
end




function myreverseKenel2(
    Imlambdap0q0,
    Imlambdap0mq0,
    Relambdap0q0,
    ImBin,
    ReBin,
    lengthy::Int,
)
    tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    I = CartesianIndices(ImBin)
    @inbounds if tid <= length(I)
        i, j = Tuple(I[tid])
        Imlambdap0q0[i, j] = ImBin[i, j]
        Imlambdap0mq0[i, j] = ImBin[i, lengthy+1-j]
        Relambdap0q0[i, j] = ReBin[i, j]
    end
    return
end

function Interpolate_Linear(Ek::T, dq0, q0min,lengthy) where {T}
    idxEk = max(min((Ek - q0min) / dq0 + 1, lengthy), 1)
    idxmEk = max(min((-Ek - q0min) / dq0 + 1, lengthy), 1)
    idxEk1 = floor(Int32, idxEk)
    idxEk2 = ceil(Int32, idxEk)
    idxmEk1 = floor(Int32, idxmEk)
    idxmEk2 = ceil(Int32, idxmEk)
    local a1::T
    local a2::T
    local a3::T
    local a4::T
    if idxEk1 == idxEk2
        a1 = 1
        a2 = 0
        a3 = 1
        a4 = 0
        return (idxEk1,idxEk2, idxmEk1, idxmEk2, a1, a2, a3, a4)
    else
        a1 = (idxEk - idxEk2) / (idxEk1 - idxEk2)
        a2 = (-idxEk + idxEk1) / (idxEk1 - idxEk2)
        a3 = (idxmEk - idxmEk2) / (idxmEk1 - idxmEk2)
        a4 = (-idxmEk + idxmEk1) / (idxmEk1 - idxmEk2)
        return (idxEk1,idxEk2, idxmEk1, idxmEk2, a1, a2, a3, a4)
    end
end

function Interpolate_Log_Linear(Ek::T, dq0_log, q0min,lengthy) where {T}
    gridlength=(lengthy-1)/2
    idxEk = max(min((log(Ek)-log(q0min)) / dq0_log + 1, gridlength), 1)+gridlength+1
    idxmEk = 2*gridlength+2-idxEk
    idxEk1 = floor(Int32, idxEk)
    idxEk2 = ceil(Int32, idxEk)
    idxmEk1 = floor(Int32, idxmEk)
    idxmEk2 = ceil(Int32, idxmEk)
    local a1::T
    local a2::T
    local a3::T
    local a4::T
    if idxEk1 == idxEk2
        a1 = 1
        a2 = 0
        a3 = 1
        a4 = 0
        return (idxEk1,idxEk2, idxmEk1, idxmEk2, a1, a2, a3, a4)
    else
        a1 = (idxEk - idxEk2) / (idxEk1 - idxEk2)
        a2 = (-idxEk + idxEk1) / (idxEk1 - idxEk2)
        a3 = (idxmEk - idxmEk2) / (idxmEk1 - idxmEk2)
        a4 = (-idxmEk + idxmEk1) / (idxmEk1 - idxmEk2)
        return (idxEk1,idxEk2, idxmEk1, idxmEk2, a1, a2, a3, a4)
    end
end


function myreverseKenel1(
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    ImBin,
    ReBin,
    lengthy::Int,
    dq0,
    q0min,
    Ek,
    Interpfun=Interpolate_Linear,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    idxEk1, idxEk2, idxmEk1, idxmEk2, a1, a2, a3, a4 = Interpfun(Ek, dq0, q0min,lengthy)
    @inbounds if i <= lengthy
        Imlambdap0mp0[i] = ImBin[i, lengthy+1-i]
        Imlambdaq0mq0[i] = ImBin[i, lengthy+1-i]
        Imlambdap0Ek[i] = ImBin[i, idxEk2] * a2 + ImBin[i, idxEk1] * a1
        Imlambdap0mEk[i] = ImBin[i, idxmEk2] * a4 + ImBin[i, idxmEk1] * a3
        Imlambdaq0Ek[i] = ImBin[i, idxEk2] * a2 + ImBin[i, idxEk1] * a1
        Imlambdaq0mEk[i] = ImBin[i, idxmEk2] * a4 + ImBin[i, idxmEk1] * a3
        Relambdap0p0[i] = ReBin[i, i]
        Relambdaq0q0[i] = ReBin[i, i]
        Relambdap0Ek[i] = ReBin[i, idxEk2] * a2 + ReBin[i, idxEk1] * a1
        Relambdaq0Ek[i] = ReBin[i, idxEk2] * a2 + ReBin[i, idxEk1] * a1
    end
    return
end




function reverse2!(
    A::LambdaGrid,
    ImBin,
    ReBin;
    lengthy::Int,
    dq0,
    q0min,
    Ek,
    threads::Int = 1000,
    Interpfun=Interpolate_Linear,
)
    @cuda threads = threads blocks = cld(length(ImBin), threads) myreverseKenel2(
        A.Imlambdap0q0,
        A.Imlambdap0mq0,
        A.Relambdap0q0,
        ImBin,
        ReBin,
        lengthy,
    )
    @cuda threads = threads blocks = cld(lengthy, threads) myreverseKenel1(
        A.Imlambdap0mp0,
        A.Imlambdaq0mq0,
        A.Imlambdap0Ek,
        A.Imlambdap0mEk,
        A.Imlambdaq0Ek,
        A.Imlambdaq0mEk,
        A.Relambdap0p0,
        A.Relambdaq0q0,
        A.Relambdap0Ek,
        A.Relambdaq0Ek,
        ImBin,
        ReBin,
        lengthy,
        dq0,
        q0min,
        Ek,
        Interpfun,
    )
end
