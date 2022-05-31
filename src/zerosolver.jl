function ZeroLPASolve(
    T,
    Npi,
    IRScale,
    UVScale;
    u0,
    progress = true,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...,
)
    function flow_LPA(du, u, h, p, k)
        h2(x) = h(p, x)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = lam4piflow(k, -u[2], p, Npi, u[1])
        du[2] = propReLPAflow(k, p, Npi, h2)
    end
    h(p, k) = zero(u0)
    alg = MethodOfSteps(RK4())
    kspan = (UVScale, IRScale)
    prob = DDEProblem(flow_LPA, u0, h, kspan, T)
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


function ZeroSolve(
    T,
    Npi,
    IRScale,
    UVScale;
    u0,
    atolint=1e-6,
    rtolint=1e-6,
    progress = true,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.5,
)
    function flow_Zero(du, u, h, p, k)
        h2(x) = h(p, x)
        #display(plot(x->-h2(x)[2],k,Λ))
        du[1] = lam4piflow(k, -u[2], p, Npi, u[1])
        du[2] = propReZeroflow(
            k,
            p,
            Npi,
            h2,
            UVScale,
            atol = atolint,
            rtol = rtolint,
        )
        println("k=",k," ","m=",-u[2])
    end
    kspan = (UVScale, IRScale)
    h(p, k) = zero(u0)
    prob = DDEProblem(flow_Zero, u0, h, kspan, T)
    alg = MethodOfSteps(RK4())
    solve(
        prob,
        alg,
        progress = progress,
        atol = atol,
        rtol = rtol,
        adaptive = adaptive,
        dtmax = dtmax,
    )
end
