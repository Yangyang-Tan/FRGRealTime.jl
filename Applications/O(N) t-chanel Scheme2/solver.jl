function tchanelSolveTwoPointParallel3(
    p0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D,Interpolations.Extrapolation},
    lamdafun::Union{Function,Spline1D,Interpolations.Extrapolation},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int,
    ϵ = 0.1;
    u0::AbstractArray{T},
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
    end
    alg = Vern9()
    p = krang.IR
    prob = ODEProblem(flow_t, u0, krang.kspan, p)
    kList = range(krang.IR, stop = krang.UV, length = numberOfParameters)
    kListperm=sortperm(randn(length(kList)))
    kListperm2=sortperm(kListperm)
    kList=kList[kListperm]
    parameterChange =
        (prob, i, repeat) -> remake(
            prob,
            tspan = (krang.UV, kList[i]),
            p = Epi(kList[i], massfun(kList[i])),
        )
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = (sol, i) -> (sol[end], false),
    )
    soltemp = solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
        abstol = config.atol,
        reltol = config.rtol,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
    )
    imdata = hcat(soltemp...)'[:, 1]
    redata = hcat(soltemp...)'[:, 2]
    kList=kList[kListperm2]
    imdata=imdata[kListperm2]
    redata=redata[kListperm2]
    imfun = Spline1D(kList, imdata)
    refun = Spline1D(kList, redata)
    imgamm2 =
        -hcubature(
            x -> twopointflow2(p0, x[1], massfun(x[1]), Temper, imfun),
            [krang.IR],
            [krang.UV],
            atol = 1e-8,
            rtol = 1e-8,
            initdiv = 400,
        )[1]
    regamm2 =
        -hcubature(
            x ->
                twopointflow2(p0, x[1], massfun(x[1]), Temper, refun) +
                twopointflow2(p0, x[1], massfun(x[1]), Temper, lamdafun, 2),
            [krang.IR],
            [krang.UV],
            atol = 1e-8,
            rtol = 1e-8,
            initdiv = 400,
        )[1] + p0^2 - massfun(krang.UV)
    return [imgamm2, regamm2]
end





function tchanelZeroSolve(
    Temper::T,
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    krang::RGscale,
    ϵ = 0.1;
    u0::AbstractArray{T},
    config::ODEConfig,
) where {T<:Number}
    function flow_Zero(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        du[1] = Imlambdaflowfun(
            0.0,
            Epi(k, -u[3]),
            k,
            -u[3],
            u[2],
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            0.0,
            Epi(k, -u[3]),
            k,
            -u[3],
            u[2],
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[3] =
            twopointflow2(0.0, k, -u[3], Temper, x -> u[2], 4) +
            twopointflow2(0.0, k, -u[3], Temper, x -> u[2], 2)
        # println("k=", k, " ", "m=", -u[3])
    end
    prob = ODEProblem(flow_Zero, u0, krang.kspan)
    alg = Vern9()
    solve(
        prob,
        alg,
        abstol = 1e-12,
        reltol = 1e-12,
        adaptive = config.adaptive,
    )
end




function tchanelLPASolve(
    Temper,
    krang::RGscale,
    massflowfun,
    lambdaflowfun;
    u0,
    config::ODEConfig,
    Nf
)
    function flow_LPA(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = massflowfun(k, u[1], p, u[2],Nf)
        du[2] = lambdaflowfun(k, u[1], p, u[2],Nf)
    end
    alg = Vern9()
    prob = ODEProblem(flow_LPA, u0, krang.kspan, Temper)
    solve(
        prob,
        alg,
        progress = config.progress,
        abstol = config.atol,
        reltol = config.rtol,
        adaptive = config.adaptive,
        dtmax = config.dtmax
    )
end


function tchanelSolveTwoPointParallel4(
    p0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D,Interpolations.Extrapolation},
    lamdafun::Union{Function,Spline1D,Interpolations.Extrapolation},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int,
    ϵ = 0.1;
    u0::AbstractArray{T},
    Nf,
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
            Nf
        )
    end
    alg = config.alg
    p = krang.IR
    prob = ODEProblem(flow_t, u0, krang.kspan, p)
    kList = range(start=krang.IR, stop = krang.UV, length = numberOfParameters)
    kListperm=sortperm(randn(length(kList)))
    kListperm2=sortperm(kListperm)
    kList=kList[kListperm]
    parameterChange =
        (prob, i, repeat) -> remake(
            prob,
            tspan = (krang.UV, kList[i]),
            p = Epi(kList[i], massfun(kList[i])),
        )
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = (sol, i) -> (sol[end], false),
    )
    soltemp = solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
        abstol = config.atol,
        reltol = config.rtol,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
    )
    imdata = hcat(soltemp...)'[:, 1]
    redata = hcat(soltemp...)'[:, 2]
    kList=kList[kListperm2]
    imdata=imdata[kListperm2]
    redata=redata[kListperm2]
    imfun = Spline1D(kList, imdata)
    refun = Spline1D(kList, redata)
    function flow_twopoint(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] =twopointflow2(p0, k, massfun(k), Temper, imfun, Nf)
        du[2] =
            twopointflow2(p0, k, massfun(k), Temper, refun, Nf) +
            twopointflow2(p0, k, massfun(k), Temper, lamdafun, 2)
    end
    prob2 =
        ODEProblem(flow_twopoint, [0.0,p0^2 - massfun(krang.UV)], krang.kspan)
    solgamm2 = solve(
        prob2,
        Tsit5(),
        abstol = 1e-12,
        reltol = 1e-12,
        dtmax = 50.0,
        adaptive = config.adaptive,
        dense = false,
        save_on = false,
        save_start = false,
        save_end = true,
    )
    return solgamm2.u[end]
    # regamm2 =
    #     -hcubature(
    #         x ->
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, refun) +
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, lamdafun,2),
    #         [krang.IR],
    #         [krang.UV],
    #         atol = 1e-6,
    #         rtol = 1e-6,
    #         initdiv = 400,
    #     )[1] + p0^2 - massfun(krang.UV)
    # return [imgamm2, regamm2]
end

function tchanelSolveTwoPointParallel5(
    p0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D,Interpolations.Extrapolation},
    lamdafun::Union{Function,Spline1D,Interpolations.Extrapolation},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    tlist::AbstractArray{T},
    ϵ = 0.1;
    u0::AbstractArray{T},
    Nf,
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
            Nf
        )
    end
    alg = config.alg
    p = krang.IR
    prob = ODEProblem(flow_t, u0, krang.kspan, p)
    kListperm=sortperm(randn(length(tlist)))
    kListperm2=sortperm(kListperm)
    kList=tlist[kListperm]
    parameterChange =
        (prob, i, repeat) -> remake(
            prob,
            tspan = (krang.UV, kList[i]),
            p = Epi(kList[i], massfun(kList[i])),
        )
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = (sol, i) -> (sol[end], false),
    )
    soltemp = solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = length(kList),
        abstol = config.atol,
        reltol = config.rtol,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
    )
    imdata = hcat(soltemp...)'[:, 1]
    redata = hcat(soltemp...)'[:, 2]
    kList=kList[kListperm2]
    imdata=imdata[kListperm2]
    redata=redata[kListperm2]
    imfun = Spline1D(kList, imdata)
    refun = Spline1D(kList, redata)
    function flow_twopoint(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] =twopointflow2(p0, k, massfun(k), Temper, imfun, Nf)
        du[2] =
            twopointflow2(p0, k, massfun(k), Temper, refun, Nf) +
            twopointflow2(p0, k, massfun(k), Temper, lamdafun, 2)
    end
    prob2 =
        ODEProblem(flow_twopoint, [0.0,p0^2 - massfun(krang.UV)], krang.kspan)
    solgamm2 = solve(
        prob2,
        Tsit5(),
        abstol = 1e-12,
        reltol = 1e-12,
        dtmax = 0.5,
        adaptive = config.adaptive,
        dense = false,
        save_on = false,
        save_start = false,
        save_end = true,
    )
    return solgamm2.u[end]
    # regamm2 =
    #     -hcubature(
    #         x ->
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, refun) +
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, lamdafun,2),
    #         [krang.IR],
    #         [krang.UV],
    #         atol = 1e-6,
    #         rtol = 1e-6,
    #         initdiv = 400,
    #     )[1] + p0^2 - massfun(krang.UV)
    # return [imgamm2, regamm2]
end




function tchanelSolveTwoPointParalleltest(
    p0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D},
    lamdafun::Union{Function,Spline1D},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int,
    ϵ = 0.1;
    u0::AbstractArray{T},
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            p0,
            p,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
    end
    alg = Vern9()
    p = krang.IR
    prob = ODEProblem(flow_t, u0, krang.kspan, p)
    kList = range(krang.IR, stop = krang.UV, length = numberOfParameters)
    parameterChange =
        (prob, i, repeat) -> remake(
            prob,
            tspan = (krang.UV, kList[i]),
            p = Epi(kList[i], massfun(kList[i])),
        )
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = (sol, i) -> (sol[end], false),
    )
    soltemp = solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
        abstol = config.atol,
        reltol = config.rtol,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
    )
    imdata = hcat(soltemp...)'[:, 1]
    redata = hcat(soltemp...)'[:, 2]
    imfun = Spline1D(kList, imdata)
    refun = Spline1D(kList, redata)
    function flow_twopoint(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] =
            twopointflow2(p0, k, massfun(k), Temper, refun, 4) +
            twopointflow2(p0, k, massfun(k), Temper, lamdafun, 2)
    end
    prob2 =
        ODEProblem(flow_twopoint, [p0^2 - massfun(krang.UV)], krang.kspan)
    solregamm2 = solve(
        prob2,
        alg,
        abstol = 1e-10,
        reltol = 1e-10,
        dtmax = 10.0,
        adaptive = config.adaptive,
    )
    # regamm2 =
    #     -hcubature(
    #         x ->
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, refun) +
    #             twopointflow2(p0, x[1], massfun(x[1]), Temper, lamdafun,2),
    #         [krang.IR],
    #         [krang.UV],
    #         atol = 1e-6,
    #         rtol = 1e-6,
    #         initdiv = 400,
    #     )[1] + p0^2 - massfun(krang.UV)
    # return [imgamm2, regamm2]
end



function tchanelSolveFourPointParallel(
    k0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D},
    lamdafun::Union{Function,Spline1D},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int,
    ϵ = 0.1;
    u0::AbstractArray{T},
    Nf=4,
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        du[1] = Imlambdaflowfun(
            p,
            Eb(k0,massfun(k0)),
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            p,
            Eb(k0,massfun(k0)),
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
            Nf
        )
    end
    alg = Tsit5()
    p = 0.0#p0 initial
    prob = ODEProblem(flow_t, u0, (krang.UV, k0), p)
    p0List = range(0.0, stop = 400.0, length = numberOfParameters)
    parameterChange = (prob, i, repeat) -> remake(prob, p = p0List[i])
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = (sol, i) -> (sol[end], false),
    )
    soltemp = solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
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
    imdata = hcat(soltemp...)'[:, 1]
    redata = hcat(soltemp...)'[:, 2]
    return [p0List, imdata, redata]
end




function tchanelSolveFourPointParallelSimple(
    q0::T,
    Temper::T,
    krang::RGscale,
    massfun::Union{Function,Spline1D},
    lamdafun::Union{Function,Spline1D},
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    ϵ = 0.1;
    u0::AbstractArray{T},
    config::ODEConfig,
) where {T<:Number}
    function flow_t(du, u, p, k)
        du[1] = 0.0*Imlambdaflowfun(
            0.0,
            q0,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
        du[2] = Relambdaflowfun(
            0.0,
            q0,
            k,
            massfun(k),
            lamdafun(k),
            Temper,
            u[1],
            u[2],
            ϵ,
        )
    end
    alg = Vern9()
    prob = ODEProblem(flow_t, u0, (krang.UV, krang.IR))
    soltemp = solve(
        prob,
        alg,
        dtmax = config.dtmax,
        abstol = config.atol,
        reltol = config.rtol,
        calck = false,
        # dt=0.0000001,
        dtmin = 1e-12,
        adaptive = config.adaptive,
        progress = config.progress,
    )
    # imdata = hcat(soltemp...)'[:, 1]
    # redata = hcat(soltemp...)'[:, 2]
    # return [soltemp.t,imdata, redata]
end
