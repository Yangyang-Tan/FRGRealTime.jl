@everywhere include(joinpath(path,"ini.jl"))
# @everywhere Tself=Tself+0.101
@everywhere include(joinpath(path,"Zero.jl"))



@time solself = solve(
    probself,
    Tsit5(),
    dense = false,
    save_on = false,
    save_start = false,
    save_timeseries=false,
    save_end = true,
    reltol = 1e-8,
    abstol = 1e-8,
    progress = true,
    dtmax=50.0,
    dtmin=1e-16,
)
# solself=1




@everywhere function getmass(Temp)
    probself2 = ODEProblem(
        flowslefZero,
        collect(vcat([116585.9, fill(6 * 8.0, length(q0grid))]...)),
        lpakrang.kspan,
        # (800.0,400.0),
        Temp,
    )
    solself3 = solve(
        probself2,
        Tsit5(),
        dense = false,
        save_on = false,
        save_start = false,
        save_timeseries=false,
        save_end = true,
        reltol = 1e-8,
        abstol = 1e-8,
        progress = true,
        dtmax=50.0,
        dtmin=1e-14,
    )
    -solself3.u[end][1]
end



massTc=pmap(getmass,Tcgrid_Tc)


plot(log.(Tcgrid_Tc.-Tc_Tc),seriestype=:scatter)

plot(log.((Tcgrid_Tc.-Tc_Tc)./Tcmin_Tc)[2:end],log.(1 ./sqrt.(massTc))[2:end],seriestype=:scatter)


massTcfun=Spline1D(log.((Tcgrid_Tc.-Tc_Tc)./Tcmin_Tc)[2:end],log.(1 ./sqrt.(massTc))[2:end])



plot(x->derivative(massTcfun,x),-7,-4.2)

-solself.u[end][1]

solselft=solself.t
solselfu1=hcat(solself.u...)[1,:]
solselfu2=hcat(solself.u...)[2,:]
-solself[end][1]
a
@eval @everywhere solselft = $solselft
@eval @everywhere solselfu1 = $solselfu1
@eval @everywhere solselfu2 = $solselfu2

@everywhere include(joinpath(path,"solver.jl"))

@everywhere tempsolmassfun=Spline1D(reverse(solselft),-reverse(solselfu1))
@everywhere tempsollambdafun=Spline1D(reverse(solselft),reverse(solselfu2))

@everywhere solmassfun = CubicSplineInterpolation(0.0:0.05:800.0,tempsolmassfun.(0.0:0.05:800.0))
@everywhere sollambdafun = CubicSplineInterpolation(0.0:0.05:800.0,tempsollambdafun.(0.0:0.05:800.0))
solmassfun(0.0)
p0grid_T=SharedArray(0.0:0.5:400.0|>collect)

@time tsolTwoPointparallelTc = pmap(p0grid_Tc) do p0
    tchanelSolveTwoPointParallel4(
        p0,
        Tself+0.007276130,
        lpakrang,
        tempsolmassfun,
        tempsollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        40*2000,
        ϵ5,
        u0 = [-1e-10, 6 * 8.0],
        config = config_spec,
    )
end


redata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,2])
imdata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,1])
sepcdata=2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2)
dycefun=Spline1D(log.(p0grid_Tc),log.(sepcdata))


plot(x->derivative(dycefun,x),-6.9,-1)

plot(log.(p0grid_Tc),log.(sepcdata))

tchanelSolveTwoPointParallel4(
    0.002,
    Tself+0.007276130,
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    40*2000,
    ϵ5,
    u0 = [-1e-10, 6 * 8.0],
    config = config_spec,
)



tchanelSolveTwoPointParallel4(
    0.00,
    Tself+0.0106555,
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    20*2000,
    ϵ5,
    u0 = [-1e-10, 6 * 8.0],
    config = config_spec,
)


plot(x->sqrt(tempsolmassfun(x)+x^2),0.01,0.1)


@everywhere config_spec.atol=1e-8
@everywhere config_spec.rtol=1e-8

redata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,2])
imdata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,1])
plot!(
    p0grid_T,
    2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2),yaxis=:log
)



using DelimitedFiles

writedlm("/home/tyy/nfs_share/FRGRealTime.jl/Applications/O(N) t-chanel Scheme2/Data/DynamicalCE/p0gridN=4.dat",p0grid_Tc)


writedlm("/home/tyy/nfs_share/FRGRealTime.jl/Applications/O(N) t-chanel Scheme2/Data/DynamicalCE/redataN=4.dat",redata_twopoint)
writedlm("/home/tyy/nfs_share/FRGRealTime.jl/Applications/O(N) t-chanel Scheme2/Data/DynamicalCE/imdataN=4.dat",imdata_twopoint)

plot(log.(p0grid_Tc),log.(imdata_twopoint))
plot(log.(p0grid_Tc[20:end]),log.(redata_twopoint[20:end]))
redata_twopoint
