@everywhere include(joinpath(path,"ini.jl"))
# @everywhere Tself=Tself+0.101
@everywhere include(joinpath(path,"Zero.jl"))
@everywhere include(joinpath(path,"biselect.jl"))
@everywhere lpakrang = RGscale(800.0, 0.1)
@everywhere config_self = ODEConfig(
    dtmax = 50.0,
    rtol = 1e-6,
    atol = 1e-6,
    progress = false,
    dense = false,
    save_on = false,
    save_start = false,
    save_end = true,
    alg=Tsit5(),
)
@everywhere Nf_ini=1
@everywhere function getlambda(α)
    solNf = tchanelLPASolve(
        Tlpa,
        lpakrang,
        massflow2,
        lambdaflow2;
        u0 = [mass2_UV, α*lambda_UV],
        config = config_lpa,
        Nf=Nf_ini
    )
    solNf.u[end][1]
end

getlambda(1.779445)
lambdaα=getmassmin(getlambda, 2.7844252,2.78443, N_kernels = 4, N_iters = 25)[2]
getlambda(lambdaα)
@eval @everywhere lambdaα = $lambdaα
#
# @everywhere config_self.atol=1e-8
# @everywhere config_self.rtol=1e-8
@everywhere getTcfun(Temper) = tchanelZeroSelfSolve(
    Temper,
    lpakrang,
    flowslefZeroNf,
    u0 = [mass2_UV, lambdaα * lambda_UV],
    Nf = Nf_ini,
    config = config_self,
)


# rmprocs(1,waitfor=0)
# nprocs()
# interrupt(4)
@everywhere config_self.dtmax=0.4
@everywhere config_self.atol=1e-5
@everywhere config_self.rtol=1e-5
Tcminmax
Tcminmax=getmassmin(getTcfun, 144.575, 144.578176, N_kernels = 4, N_iters = 10)
Tcminmax=getmassmin(getTcfun, Tcminmax[1], Tcminmax[2], N_kernels = 4, N_iters = 1)



solself=tchanelZeroSelfSolve2(
    Tcminmax[2],
    lpakrang,
    flowslefZeroNf,
    u0 = [mass2_UV, lambdaα * lambda_UV],
    Nf = Nf_ini,
    config = config_self,
)
solselft=solself.t
solselfu1=hcat(solself.u...)[1,:]
solselfu2=hcat(solself.u...)[2,:]

-solself[end][1]

@eval @everywhere solselft = $solselft
@eval @everywhere solselfu1 = $solselfu1
@eval @everywhere solselfu2 = $solselfu2

@everywhere tempsolmassfun=Spline1D(reverse(solselft),-reverse(solselfu1))
@everywhere tempsollambdafun=Spline1D(reverse(solselft),reverse(solselfu2))

@everywhere solmassfun = CubicSplineInterpolation(0.0:0.5:800.0,tempsolmassfun.(0.0:0.5:800.0))
@everywhere sollambdafun = CubicSplineInterpolation(0.0:0.5:800.0,tempsollambdafun.(0.0:0.5:800.0))

writedlm(
    joinpath(path, "N=1/realtime_self_N=1_k.dat"),
    solselft,
)

writedlm(
    joinpath(path, "N=1/realtime_self_N=1_m2k.dat"),
    -solselfu1,
)
writedlm(
    joinpath(path, "N=1/realtime_self_N=1_lambdak.dat"),
    solselfu2,
)
writedlm(
    joinpath(path, "N=1/realtime_self_N=1_α.dat"),
    lambdaα,
)


massselfT=collect(Tcminmax[2]:1.0:400.0)
massslpaT=collect(145:1.0:400.0)
massselfdata= @showprogress pmap(getTcfun,massselfT,batch_size=1)

masslpadata=map(massslpaT) do temp
    tchanelLPASolve(
        temp,
        lpakrang,
        massflow2,
        lambdaflow2;
        u0 = [-116585.9, lambdaα*6 * 8.0],
        config = config_lpa,
        Nf=2.0
    )[end][1]
end

plot(massselfT,sqrt.(massselfdata))
plot!(massslpaT,sqrt.(masslpadata))

writedlm(
    joinpath(path, "N=1/massselfT_N=1.dat"),
    massselfT,
)

writedlm(
    joinpath(path, "N=1/massself_N=1.dat"),
    sqrt.(massselfdata),
)

writedlm(
    joinpath(path, "N=2/masslpaT_N=2.dat"),
    massselfT,
)

writedlm(
    joinpath(path, "N=2/masslpa_N=2.dat"),
    sqrt.(massselfdata),
)
