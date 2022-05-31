using Distributed
addprocs(2)
nprocs()
@everywhere using FRGRealTime
@everywhere using FRGRealTime.ThresholdFunctions
@everywhere using DifferentialEquations
@everywhere using HCubature
using Plots
@everywhere include("flow.jl")
include("solver.jl")


@everywhere function tchanelLPASolve(
    T,
    IRScale,
    UVScale,
    massflowfun,
    lambdaflowfun;
    u0,
    progress = true,
    atol = 1e-14,
    rtol = 1e-14,
    adaptive = true,
    dtmax = 0.01,
    kwargs...,
)
    function flow_LPA(du, u, p, k)
        #display(plot(x->-h2(x)[2],k,Λ))
        #println(k, " m2=", -u[2], " λ=", u[1])
        du[1] = massflowfun(k, u[1], p, u[2])
        du[2] = lambdaflowfun(k, u[1], p, u[2])
    end
    alg = Tsit5()
    kspan = (UVScale, IRScale)
    prob = ODEProblem(flow_LPA, u0, kspan, T)
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





@everywhere sol1 = tchanelLPASolve(
    150.0,
    1.0,
    800.0,
    massflow2,
    lambdaflow2;
    u0 = [-116585.3, 6 * 8.0],
)



@everywhere sol1 = tchanelLPASolve(
    145.0,
    1.0,
    1000.0,
    massflow1,
    lambdaflow1;
    u0 = [-146585.3, 5 * 7.0],
)


plot(k -> sol1(k)[1], 1.0, 800.0)


@time tsol1 = tchanelSolve(
    2000.0,
    145.0,
    1.0,
    1000.0,
    x -> sol1(x)[1],
    Imlambdaflow1,
    Relambdaflow1,
    0.1,
    u0 = [-0.0001, 5 * 7.0],
    dtmax = 0.5,
)


@time tsolparallel1 = tchanelSolveParallel(
    145.0,
    1.0,
    1000.0,
    x -> sol1(x)[1],
    Imlambdaflow1,
    Relambdaflow1,
    4000,
    0.2,
    u0 = [-0.0001, 5 * 7.0],
    dtmax = 0.05,
)




@time tsolTwoPointparallel1 = tchanelSolveTwoPointParallel(
    150.0,
    1.0,
    800.0,
    x -> sol1(x)[1],
    x -> sol1(x)[2],
    Imlambdaflow2,
    Relambdaflow2,
    100,
    0.05,
    u0 = [-0.00001, 6 * 8.0],
    dtmax = 0.04,
)

imvec=hcat((tsolTwoPointparallel1.u)...)'[:,1]
revec=hcat((tsolTwoPointparallel1.u)...)'[:,2]

lambda_omega = range(0.0, stop = 400.0, length = 100)


plot(lambda_omega,revec)

plot!(lambda_omega,abs.(revec),yaxis=:log)


plot(lambda_omega,-imvec./(imvec.^2+revec.^2),yaxis=:log)
tsolparallel1[3](40.0)

plot(k -> sol1(k)[1], 1.0, 800.0)


plot(k -> -tsol1(k)[1], 1, 1000.0)


plot(k -> tsol1(k)[2], 1.0, 800.0, yaxis = :log)
