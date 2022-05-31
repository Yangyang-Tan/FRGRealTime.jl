struct RGscale{T<:Number}
    IR::T
    UV::T
    kspan::Tuple{T,T}
end
RGscale(UV,IR)=RGscale(IR,UV,(UV,IR))




mutable struct ODEConfig
    atol::AbstractFloat
    rtol::AbstractFloat
    dtmax::AbstractFloat
    maxiters::Int
    progress::Bool
    adaptive::Bool
    dense::Bool
    save_on::Bool
    save_start::Bool
    save_end::Bool
end



function ODEConfig(;
    atol = 1e-14,
    rtol = 1e-14,
    dtmax = 0.1,
    maxiters = 10^5,
    progress = true,
    adaptive = true,
    dense = true,
    save_on = true,
    save_start = true,
    save_end = true,
)
    ODEConfig(
        atol,
        rtol,
        dtmax,
        maxiters,
        progress,
        adaptive,
        dense,
        save_on,
        save_start,
        save_end,
    )
end
