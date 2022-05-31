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
    reltol = 1e-8,
    abstol = 1e-8,
    progress = true,
    dtmax=50.0,
    dtmin=1e-16,
)
-solself.u[end][1]





@everywhere function getmass(Temp)
    probself2 = ODEProblem(
        flowslefZeroNf,
        collect(vcat([0.3509389062993*116585.9, fill(6 * 8.0, length(q0grid))]...)),
        (800.0,0.1),
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

Tcgrid_Tc



Tself+0.001939
91.371010001
Tc_Tc=
Tcmin_Tc=Tc_Tc+0.00005
Tcmax_Tc=Tcmin_Tc+100
N_Tcgrid_Tc=120
Tcstep_Tc = (log(Tcmax_Tc-Tc_Tc) - log(Tcmin_Tc-Tc_Tc)) / (N_Tcgrid_Tc-1)
Tcgrid_Tc = SharedArray(Tc_Tc .+exp.(collect(log(Tcmin_Tc-Tc_Tc):Tcstep_Tc:log(Tcmax_Tc-Tc_Tc))))




massTc=pmap(getmass,Tcgrid_Tc)
plot(log.((Tcgrid_Tc.-Tc_Tc)./Tcmin_Tc)[2:end],log.(1 ./sqrt.(massTc))[2:end],seriestype=:scatter)
massTcfun=Spline1D(log.((Tcgrid_Tc.-Tc_Tc)./Tc_Tc),log.(1 ./sqrt.(massTc)))

plot(x->derivative(massTcfun,x),-16.2,-5.2)
plot(solself.t,-hcat(solself.u...)[1,:])
plot(solself.t,hcat(solself.u...)[2,:])

writedlm(joinpath(path,"N=0/massN=0.data"),massTc)
writedlm(joinpath(path,"N=0/TcgridN=0.data"),Tcgrid_Tc)
