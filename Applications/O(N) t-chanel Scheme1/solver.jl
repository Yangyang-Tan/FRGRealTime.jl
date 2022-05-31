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
    Relambdaflowfun, numberOfParameters,
    ϵ = 0.1;
    u0,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...
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
    Relambdaflowfun, numberOfParameters,
    ϵ = 0.1;
    u0,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...
)
    function flow_t(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = Imlambdaflowfun(p, k, massfun(k),lamdafun(k), T, u[1], u[2], ϵ)
        du[2] = Relambdaflowfun(p, k, massfun(k),lamdafun(k), T, u[1], u[2], ϵ)
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
    function imrefun(sol,i)
        p0=parameterList[i];
        imgamm2=-hcubature(
            x -> twopointflow2(p0, x[1], massfun(x[1]), T, y->sol(y)[1]),
            [1.0],
            [800.0],
            atol = 1e-6,
            rtol = 1e-6,
        )[1];
        regamm2=hcubature(
            x -> twopointflow2(p0, x[1], massfun(x[1]), T,  y->sol(y)[2]),
            [1.0],
            [800.0],
            atol = 1e-6,
            rtol = 1e-6,
        )[1]  +p0^2 - massfun(800.0);
        return ([imgamm2,regamm2],false)
    end
    ensembleProb = EnsembleProblem(prob, prob_func = parameterChange,output_func=imrefun)
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
