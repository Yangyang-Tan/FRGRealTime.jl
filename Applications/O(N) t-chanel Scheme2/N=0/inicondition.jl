@everywhere include(joinpath(path,"ini.jl"))
# @everywhere Tself=Tself+0.101
@everywhere include(joinpath(path,"Zero.jl"))

@everywhere function getlambda(α)
    solNf = tchanelLPASolve(
        Tlpa,
        lpakrang,
        massflow2,
        lambdaflow2;
        u0 = [-116585.9, α*6 * 8.0],
        config = config_lpa,
        Nf=0
    )
    solNf.u[end][1]
end
getlambda(5.6417581372728721)
@everywhere lambdaα=5.6417581372728721
solNf0 = tchanelLPASolve(
    Tlpa,
    lpakrang,
    massflow2,
    lambdaflow2;
    u0 = [-116585.9, lambdaα*6 * 8.0],
    config = config_lpa,
    Nf=0
)
plot(x->solNf0(x)[1],1,800)

solself=tchanelZeroSelfSolve2(
    Tlpa,
    lpakrang,
    flowslefZeroNf,
    [-116585.9, lambdaα * 6 * 8.0],
    0.0,
)


solselft=solself.t
solselfu1=hcat(solself.u...)[1,:]
solselfu2=hcat(solself.u...)[2,:]

-solself[end][1]

@eval @everywhere solselft = $solselft
@eval @everywhere solselfu1 = $solselfu1
@eval @everywhere solselfu2 = $solselfu2

@everywhere include(joinpath(path,"solver.jl"))

@everywhere tempsolmassfun=Spline1D(reverse(solselft),-reverse(solselfu1))
@everywhere tempsollambdafun=Spline1D(reverse(solselft),reverse(solselfu2))

@everywhere solmassfun = CubicSplineInterpolation(0.0:0.5:800.0,tempsolmassfun.(0.0:0.5:800.0))
@everywhere sollambdafun = CubicSplineInterpolation(0.0:0.5:800.0,tempsollambdafun.(0.0:0.5:800.0))

writedlm(
    joinpath(path, "N=0/realtime_self_N=0_k.dat"),
    solselft,
)

writedlm(
    joinpath(path, "N=0/realtime_self_N=0_m2k.dat"),
    -solselfu1,
)
writedlm(
    joinpath(path, "N=0/realtime_self_N=0_lambdak.dat"),
    solselfu2,
)
writedlm(
    joinpath(path, "N=0/realtime_self_N=0_α.dat"),
    lambdaα,
)


massselfT=collect(145:1.0:400.0)
massslpaT=collect(145:2.0:400.0)
massselfdata=pmap(massselfT) do temp
    tchanelZeroSelfSolve(
        temp,
        lpakrang,
        flowslefZeroNf,
        [-116585.9, lambdaα * 6 * 8.0],
        0.0,
    )
end

masslpadata=map(massslpaT) do temp
    tchanelLPASolve(
        temp,
        lpakrang,
        massflow2,
        lambdaflow2;
        u0 = [-116585.9, lambdaα*6 * 8.0],
        config = config_lpa,
        Nf=0
    )[end][1]
end

plot(massselfT,sqrt.(massselfdata))
plot(massslpaT,sqrt.(masslpadata))

writedlm(
    joinpath(path, "N=0/massselfT_N=0.dat"),
    massselfT,
)

writedlm(
    joinpath(path, "N=0/massself_N=0.dat"),
    sqrt.(massselfdata),
)

writedlm(
    joinpath(path, "N=0/masslpaT_N=0.dat"),
    massselfT,
)

writedlm(
    joinpath(path, "N=0/masslpa_N=0.dat"),
    sqrt.(massselfdata),
)
