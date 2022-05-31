function tchanelSolve(
    p0,
    T,
    IRScale,
    UVScale,
    massfun,
    Imlambdaflowfun,
    Relambdaflowfun,
    ϵ = 0.1;
    u0,
    progress = true,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...,
)
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(p, k, massfun(k), T, u[1], u[2], ϵ)
        du[2] = Relambdaflowfun(p, k, massfun(k), T, u[1], u[2], ϵ)
    end
    alg = Tsit5()
    kspan = (UVScale, IRScale)
    prob = ODEProblem(flow_t, u0, kspan, p0)
    solve(
        prob,
        alg,
        progress = progress,
        atol = atol,
        rtol = rtol,
        adaptive = adaptive,
        dtmax = dtmax,
        kwargs...,
    )
end



function tchanelSolveParallel(
    T,
    IRScale,
    UVScale,
    massfun,
    Imlambdaflowfun,
    Relambdaflowfun,
    numberOfParameters,
    ϵ = 0.1;
    u0,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...,
)
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(p, k, massfun(k), T, u[1], u[2], ϵ)
        du[2] = Relambdaflowfun(p, k, massfun(k), T, u[1], u[2], ϵ)
    end
    alg = Tsit5()
    kspan = (UVScale, IRScale)
    p = 1.0
    prob = ODEProblem(flow_t, u0, kspan, p)
    parameterList = range(1.0, stop = 500.0, length = numberOfParameters)
    # function parameterChange(prob, i, repeat)
    #     prob.p = parameterList[i]
    #     prob
    # end
    parameterChange = (prob, i, repeat) -> remake(prob, p = parameterList[i])
    ensembleProb = EnsembleProblem(prob, prob_func = parameterChange)
    solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
        atol = atol,
        rtol = rtol,
        adaptive = adaptive,
        dtmax = dtmax,
        kwargs...,
    )
end


function tchanelSolveTwoPointParallel(
    T,
    IRScale,
    UVScale,
    massfun,
    lamdafun,
    Imlambdaflowfun,
    Relambdaflowfun,
    numberOfParameters,
    ϵ = 0.1;
    u0,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...,
)
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(p, k, massfun(k), lamdafun(k), T, u[1], u[2], ϵ)
        du[2] = Relambdaflowfun(p, k, massfun(k), lamdafun(k), T, u[1], u[2], ϵ)
    end
    alg = Tsit5()
    kspan = (UVScale, IRScale)
    p = 0.0
    prob = ODEProblem(flow_t, u0, kspan, p)
    parameterList = range(0, stop = 400.0, length = numberOfParameters)
    # function parameterChange(prob, i, repeat)
    #     prob.p = parameterList[i]
    #     prob
    # end
    parameterChange = (prob, i, repeat) -> remake(prob, p = parameterList[i])
    function imrefun(sol, i)
        p0 = parameterList[i]
        imgamm2 =
            -hcubature(
                x -> twopointflow2(p0, x[1], massfun(x[1]), T, y -> sol(y)[1]),
                [1.0],
                [800.0],
                atol = 1e-6,
                rtol = 1e-6,
            )[1]
        regamm2 =
            hcubature(
                x -> twopointflow2(p0, x[1], massfun(x[1]), T, y -> sol(y)[2]),
                [1.0],
                [800.0],
                atol = 1e-6,
                rtol = 1e-6,
            )[1] + p0^2 - massfun(800.0)
        return ([imgamm2, regamm2], false)
    end
    ensembleProb = EnsembleProblem(
        prob,
        prob_func = parameterChange,
        output_func = imrefun,
    )
    solve(
        ensembleProb,
        alg,
        EnsembleThreads(),
        trajectories = numberOfParameters,
        atol = atol,
        rtol = rtol,
        adaptive = adaptive,
        dtmax = dtmax,
        kwargs...,
    )
end

function tchanelSolveTwoPointParallel3(
    p0::T,
    Temper::T,
    krang::RGscale,
    massfun::Function,
    lamdafun::Function,
    Imlambdaflowfun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int,
    ϵ = 0.1;
    u0::AbstractArray{T},
    config::ODEConfig
) where {T<:Number}
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] =
            Imlambdaflowfun(p0, p, k, massfun(k), lamdafun(k), Temper, u[1], u[2], ϵ)
        du[2] =
            Relambdaflowfun(p0, p, k, massfun(k), lamdafun(k), Temper, u[1], u[2], ϵ)
    end
    alg = Tsit5()
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
        EnsembleSerial(),
        trajectories = numberOfParameters,
        atol = config.atol,
        rtol = config.rtol,
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
    imgamm2 =
        -hcubature(
            x -> twopointflow2(p0, x[1], massfun(x[1]), Temper, imfun),
            [krang.IR],
            [krang.UV],
            atol = 1e-6,
            rtol = 1e-6,
            initdiv = 100,
        )[1]
    regamm2 =
        hcubature(
            x -> twopointflow2(p0, x[1], massfun(x[1]), Temper, refun),
            [krang.IR],
            [krang.UV],
            atol = 1e-6,
            rtol = 1e-6,
            initdiv = 100,
        )[1] + p0^2 - massfun(krang.UV)
    return [imgamm2, regamm2]
end



function tchanelSolveFourPointParallel(
    k0::T,
    Temper::T,
    krang::RGscale,
    massfun::Function,
    lamdafun::Function,
    Relambdaflowfun::Function,
    numberOfParameters::Int;
    u0::AbstractArray{T},
    config::ODEConfig
) where {T<:Number}
    function flow_t(du, u, p, k)
        du[1] =Relambdaflowfun(p, Epi(k0, massfun(k0)), k, massfun(k), lamdafun(k), Temper, u[1])
    end
    alg = Tsit5()
    p = 0.0#p0 initial
    prob = ODEProblem(flow_t, u0, (krang.UV,k0), p)
    p0List = range(0.0, stop = 1600.0, length = numberOfParameters)
    parameterChange =
        (prob, i, repeat) -> remake(
            prob,
            p = p0List[i],
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
        dtmax=config.dtmax,
        atol = config.atol,
        rtol = config.rtol,
        adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
    )
    return [p0List,hcat(soltemp...)[1,:]]
end
