function tchanelSolveFourPointGPU(
    k0::T,
    Temper::T,
    krang::RGscale,
    ϵ = 0.1;
    u0::AbstractArray{T},
    Nf = 4,
    config::ODEConfig,
) where {T<:Number}
    function flow_tgpu(du, u, p, k)
        Ek = 0.9 * k
        spacing = p
        Imlambda = @view u[:, :, 1]
        Relambda = @view u[:, :, 2]
        idx = max(min(round(Int64, spacing * Ek), length(q0gridGPU)), 1)
        println("idx=", idx)
        Imlambdap0Ek = Imlambda[:, idx] * 0.0f0
        Relambdap0Ek = Relambda[:, idx]
        Imlambdaq0Ek = Imlambdap0Ek'
        Relambdaq0Ek = Relambdap0Ek'
        Imlambdap0p0 = diag(Imlambda)
        Imlambdaq0q0 = diag(Imlambda)' * 0.0f0
        Relambdap0p0 = diag(Relambda)
        Relambdaq0q0 = diag(Relambda)'
        println("k=", k)
        du[:, :, 1] .=
            Imlambdaflow3.(
                p0gridGPU,
                q0gridGPU,
                k,
                0.0f0,
                Temper,
                Imlambda,
                Relambda,
                Imlambdap0p0,
                Imlambdaq0q0,
                Imlambdap0Ek,
                Imlambdaq0Ek,
                Relambdap0p0,
                Relambdaq0q0,
                Relambdap0Ek,
                Relambdaq0Ek,
                ϵ,
                Nf,
            )
        du[:, :, 2] .=
            Relambdaflow3.(
                p0gridGPU,
                q0gridGPU,
                k,
                0.0f0,
                Temper,
                Imlambda,
                Relambda,
                Imlambdap0p0,
                Imlambdaq0q0,
                Imlambdap0Ek,
                Imlambdaq0Ek,
                Relambdap0p0,
                Relambdaq0q0,
                Relambdap0Ek,
                Relambdaq0Ek,
                ϵ,
                Nf,
            )
    end
    alg = Tsit5()
    p = 1 / (p0gridGPU[2] - p0gridGPU[1])#p0 initial
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, k0), p)
    soltemp = solve(
        prob,
        alg,
        dtmax = config.dtmax,
        abstol = config.atol,
        reltol = config.rtol,
        calck = false,
        # dt=0.0000001,
        dtmin = config.dtmin,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
        progress = config.progress,
    )
end



function tchanelSolveTwoPointGPU(
    k0::T,
    Temper::T,
    krang::RGscale,
    ϵ = 0.1;
    u0::AbstractArray{T},
    Nf = 4,
    config::ODEConfig,
) where {T<:Number}
    function flow_twopoint_gpu(du, u, p, k)
        Ek = Eb(k, massfun(k))
        spacing = p
        # prepare u & du
        ImGamma2 = @view u[:, 1, 1]
        ReGamma2 = @view u[:, 1, 2]
        dkImGamma2 = @view du[:, 1, 1]
        dkReGamma2 = @view du[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]
        #end
        idx = min(round(Int64, spacing * Ek) + 1, length(q0gridGPU))
        println("idx=", idx)
        Imlambdap0Ek = Imlambda[:, idx]
        Relambdap0Ek = Relambda[:, idx]
        Imlambdaq0Ek = Imlambdap0Ek'
        Relambdaq0Ek = Relambdap0Ek'
        Imlambdap0p0 = diag(Imlambda)
        Imlambdaq0q0 = diag(Imlambda)'
        Relambdap0p0 = diag(Relambda)
        Relambdaq0q0 = diag(Relambda)'
        println("k=", k)
        dkImlambda .=
            Imlambdaflow3.(
                p0gridGPU,
                q0gridGPU,
                k,
                0.0f0,
                Temper,
                Imlambda,
                Relambda,
                Imlambdap0p0,
                Imlambdaq0q0,
                Imlambdap0Ek,
                Imlambdaq0Ek,
                Relambdap0p0,
                Relambdaq0q0,
                Relambdap0Ek,
                Relambdaq0Ek,
                ϵ,
                Nf,
            )
        dkRelambda .=
            Relambdaflow3.(
                p0gridGPU,
                q0gridGPU,
                k,
                0.0f0,
                Temper,
                Imlambda,
                Relambda,
                Imlambdap0p0,
                Imlambdaq0q0,
                Imlambdap0Ek,
                Imlambdaq0Ek,
                Relambdap0p0,
                Relambdaq0q0,
                Relambdap0Ek,
                Relambdaq0Ek,
                ϵ,
                Nf,
            )
        dkImGamma2 = twopointflow.(k, m2, Temper, lambda, Imlambdap0Ek, Nf)
        dkReGamma2 = twopointflow.(k, m2, Temper, lambda, Relambdap0Ek, Nf)
    end
    alg = Tsit5()
    p = 1 / (p0gridGPU[2] - p0gridGPU[1])#p0 initial
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, k0), p)
    soltemp = solve(
        prob,
        alg,
        dtmax = config.dtmax,
        abstol = config.atol,
        reltol = config.rtol,
        calck = false,
        # dt=0.0000001,
        dtmin = config.dtmin,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
        progress = config.progress,
    )
end


solself = tchanelSolveFourPointCPU(
    798.0,
    145.0,
    lpakrang,
    ϵ5,
    u0 = init_fourpoint_2d(p0grid, q0grid),
    Nf = 4.0,
    config = config_spec,
)
