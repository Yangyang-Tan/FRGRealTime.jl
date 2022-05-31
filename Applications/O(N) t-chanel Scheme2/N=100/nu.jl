@everywhere include(joinpath(path,"ini.jl"))
# @everywhere Tself=Tself+0.101
@everywhere include(joinpath(path,"Zero.jl"))

@time solself = solve(
    probselfNf,
    Tsit5(),
    dense = false,
    save_on = false,
    save_start = false,
    save_timeseries=false,
    save_end = true,
    reltol = 1e-5,
    abstol = 1e-5,
    progress = true,
    dtmax=50.0,
    dtmin=1e-16,
)

-solself.u[end][1]


@everywhere function getmass(Temp)
    probself2 = ODEProblem(
        flowslefZeroNf,
        collect(vcat([116585.9*sqrt(16.5/6), fill(6 * 8.0*sqrt(12/16.5), length(q0grid))]...)),
        (800.0,0.1),
        # (800.0,400.0),
        [Temp,10.0],
    )
    solself3 = solve(
        probself2,
        Tsit5(),
        dense = false,
        save_on = false,
        save_start = false,
        save_timeseries=false,
        save_end = true,
        reltol = 1e-7,
        abstol = 1e-7,
        progress = true,
        dtmax=50.0,
        dtmin=1e-14,
    )
    -solself3.u[end][1]
end



getmass(140.0)

mytempmin=140.0
mytempmax=145.0

mytempminmax=getmassmin(getmass,mytempmin, mytempmax, N_kernels=3, N_iters=9)

mytempminmax=getmassmin(getmass,mytempminmax[1], mytempminmax[2], N_kernels=3, N_iters=2)




Tc_Tc=mytempminmax[2]
Tcmin_Tc=Tc_Tc+0.01
Tcmax_Tc=Tcmin_Tc+2
N_Tcgrid_Tc=30
Tcstep_Tc = (log(Tcmax_Tc-Tc_Tc) - log(Tcmin_Tc-Tc_Tc)) / (N_Tcgrid_Tc-1)
Tcgrid_Tc = SharedArray(Tc_Tc .+exp.(collect(log(Tcmin_Tc-Tc_Tc):Tcstep_Tc:log(Tcmax_Tc-Tc_Tc))))


massTc=pmap(getmass,Tcgrid_Tc)
plot(log.((Tcgrid_Tc.-Tc_Tc)./Tcmin_Tc),log.(1 ./sqrt.(massTc)),seriestype=:scatter)
massTcfun=Spline1D(log.((Tcgrid_Tc.-Tc_Tc)./Tc_Tc),log.(1 ./sqrt.(massTc)))

plot(x->derivative(massTcfun,x),-9.572,-4.56)
plot(solself.t,-hcat(solself.u...)[1,:])
plot(solself.t,hcat(solself.u...)[2,:])
