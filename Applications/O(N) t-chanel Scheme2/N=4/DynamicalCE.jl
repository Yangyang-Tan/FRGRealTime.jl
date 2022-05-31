@eval @everywhere Tcminmax = $Tcminmax
@everywhere config_spec=deepcopy(config_self)

@everywhere config_spec.atol=1e-9
@everywhere config_spec.rtol=1e-9


@everywhere println(config_spec.atol)

# @everywhere gettwopointminfun(dTc) = tchanelSolveTwoPointParallel4(
#     0.0,
#     Tcminmax[2] - dTc,
#     lpakrang,
#     tempsolmassfun,
#     tempsollambdafun,
#     Imlambdaflow3,
#     Relambdaflow3,
#     8 * 2000,
#     ϵ5,
#     u0 = [-1e-10, lambdaα * lambda_UV],
#     Nf = Nf_ini,
#     config = config_spec,
# )[2]


solselftinput=reverse(deepcopy(solselft))
itp = interpolate(solselftinput, BSpline(Linear()))
@eval @everywhere itp = $itp
@eval @everywhere solselftinput = $solselftinput

@everywhere gettwopointminfun(dTc) = tchanelSolveTwoPointParallel5(
    0.0,
    Tcminmax[2] - dTc,
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    itp(1:0.2:length(solselftinput)),
    ϵ5,
    u0 = [-1e-15, lambdaα * lambda_UV],
    Nf = Nf_ini,
    config = config_spec,
)[2]

twopointReminmax =
    getmassmin(gettwopointminfun, -0.001, 0.001, N_kernels = 4, N_iters = 13)


tchanelSolveTwoPointParallel5(
    0.0001,
    Tcminmax[2] - twopointReminmax[1],
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    itp(1:0.2:length(solselftinput)),
    ϵ5,
    u0 = [-1e-12, lambdaα * lambda_UV],
    Nf = Nf_ini,
    config = config_spec,
)

@eval @everywhere twopointReminmax = $twopointReminmax

@time tsolTwoPointparallelTc = pmap(p0grid_Tc) do p0
    tchanelSolveTwoPointParallel5(
        p0,
        Tcminmax[2] - twopointReminmax[1],
        lpakrang,
        tempsolmassfun,
        tempsollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        itp(1:0.2:length(solselftinput)),
        ϵ5,
        u0 = [-1e-15, lambdaα * lambda_UV],
        Nf = Nf_ini,
        config = config_spec,
    )
end

p0grid_T=SharedArray(collect(0.1:0.5:400.0))
@time tsolTwoPointparallelT = @showprogress pmap(p0grid_T) do p0
    tchanelSolveTwoPointParallel5(
        p0,
        250.0,
        lpakrang,
        tempsolmassfun,
        tempsollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        itp(1:0.2:length(solselftinput)),
        ϵ5,
        u0 = [-1e-15, lambdaα * lambda_UV],
        Nf = Nf_ini,
        config = config_spec,
    )
end


@everywhere GC.gc()
@everywhere GC.gc(true)
redata_twopoint=(hcat(tsolTwoPointparallelT...)'[:,2])
imdata_twopoint=(hcat(tsolTwoPointparallelT...)'[:,1])


sepcdata=2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2)



dycefun=Spline1D(log.(p0grid_Tc),log.(sepcdata))


plot(x->derivative(dycefun,x),-9.2,-7.1)



plot(p0grid_T,log.(sepcdata))

plot(log.(p0grid_T),log.(imdata_twopoint))



using DelimitedFiles

writedlm(joinpath(path,"N=4/redataN=4.dat"),redata_twopoint)
writedlm(joinpath(path,"N=4/imdataN=4.dat"),imdata_twopoint)
writedlm(joinpath(path,"N=4/p0dataN=4.dat"),p0grid_Tc)

writedlm(joinpath(path,"N=4/redataN=4T=250.dat"),redata_twopoint)
writedlm(joinpath(path,"N=4/imdataN=4T=250.dat"),imdata_twopoint)
writedlm(joinpath(path,"N=4/p0dataN=4T=250.dat"),p0grid_T)
