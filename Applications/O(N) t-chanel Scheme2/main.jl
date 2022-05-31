using Distributed
addprocs(["xeon6146b","xeon6146c","xeon6146d"],tunnel=true)
nprocs()
@everywhere println(Threads.nthreads())
@everywhere using FRGRealTime
@everywhere import FRGRealTime.ThresholdFunctions as TF
@everywhere using DifferentialEquations
@everywhere using HCubature
@everywhere using Dierckx, ThreadsX, Referenceables, LoopVectorization, Interpolations
using LaTeXStrings
using Plots
@everywhere using BenchmarkTools, SharedArrays, Interpolations, Tullio
@everywhere path="/home/tyy/nfs_share/FRGRealTime.jl/Applications/O(N) t-chanel Scheme2/"
@everywhere include(joinpath(path,"struct.jl"))
@everywhere include(joinpath(path,"flow.jl"))

@everywhere include(joinpath(path,"solver.jl"))
@everywhere config_lpa = ODEConfig(dtmax = 100.0)
@everywhere config_spec = ODEConfig(dtmax = 50.0, progress = false, dense = false, save_on = false, save_start = false, save_end = true)
@everywhere include(joinpath(path,"ini.jl"))
@everywhere Tself=Tself+0.05
@everywhere include(joinpath(path,"Zero.jl"))

@time @everywhere sol1 = tchanelLPASolve(
    Tlpa,
    lpakrang,
    massflow2,
    lambdaflow2;
    u0 = [-116585.9, 6 * 8.0],
    config = config_lpa,
)

@time @everywhere solself = solve(
    probself,
    Vern9(),
    # dense = false,
    # save_on = false,
    # save_start = false,
    # save_end = true,
    reltol = 1e-8,
    abstol = 1e-8,
    dtmax=0.5,
    dtmin=1e-14,
)

solself.t
Tself
-solself.u[end][1]

@time @everywhere tsolall = tchanelZeroSolve(
    Tself,
    Imlambdaflow3,
    Relambdaflow3,
    lpakrang,
    0.01,
    u0 = [-0.0000001, 6 * 8.0, 116585.9],
    config = config_spec,
)
solself.t
@everywhere tempsolmassfun=Spline1D(reverse(solself.t),-reverse(hcat(solself.u...)[1,:]))
@everywhere tempsollambdafun=Spline1D(reverse(solself.t),reverse(hcat(solself.u...)[2,:]))

@everywhere solmassfun = CubicSplineInterpolation(0.0:0.2:800.0,tempsolmassfun.(0.0:0.2:800.0))
@everywhere sollambdafun = CubicSplineInterpolation(0.0:0.2:800.0,tempsollambdafun.(0.0:0.2:800.0))


@benchmark tempsolmassfun.(1:0.02:800)
@benchmark solmassfun.(1:0.02:800)


lines(tsolall.t,-hcat(tsolall.u...)[3,:])
lines!(sol1.t,hcat(sol1.u...)[1,:])

lines(tsolall.t,hcat(tsolall.u...)[2,:])
lines!(sol1.t,hcat(sol1.u...)[2,:])

writedlm(
    joinpath(path, "Data/FourPoint/kselfT=145.dat"),
    tsolall.t,
)

writedlm(
    joinpath(path, "Data/FourPoint/relambdaselfT=145.dat"),
    hcat(tsolall.u...)[2,:],
)


writedlm(
    joinpath(path, "Data/FourPoint/imlambdaselfT=145.dat"),
    hcat(tsolall.u...)[1,:],
)
writedlm(
    joinpath(path, "Data/FourPoint/m2selfT=145.dat"),
    -hcat(tsolall.u...)[3,:],
)

@everywhere solmassfun=Spline1D(reverse(sol1.t),reverse(hcat(sol1.u...)[1,:]))
@everywhere sollambdafun=Spline1D(reverse(sol1.t),reverse(hcat(sol1.u...)[2,:]))

testa=rand(1000)

@benchmark solmassfun.(testa)
myfun1=x->sol1(x)[1]
@benchmark myfun1.(testa)

plot(k -> sol1(k)[1], 0.001, 1.0)


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





@time tsolTwoPointparallel2 = tchanelSolveTwoPointParallel2(
    100.0,
    150.0,
    1.0,
    800.0,
    x -> sol1(x)[1],
    x -> sol1(x)[2],
    Imlambdaflow3,
    Relambdaflow3,
    1600,
    0.2,
    u0 = [-0.00001, 6 * 8.0],
)

using SharedArrays
v1=SharedArray(0.00000:0.01:0.05 |>collect)
config_spec
@time tsolTwoPointparallelTc = pmap(p0grid_Tc[375:400]) do p0
    tchanelSolveTwoPointParallel3(
        p0,
        Tself,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        6400,
        0.1,
        u0 = [-1e-8, 6 * 8.0],
        config = config_spec,
    )
end



@time tsolTwoPointparallelTc = pmap(p0grid_Tc) do p0
    tchanelSolveTwoPointParallel4(
        p0,
        Tself,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        32*8000,
        ϵ5,
        u0 = [-1e-8, 6 * 8.0],
        config = config_spec,
    )
end
config_spec

p0grid_Tc[330]
solmassfun
tchanelSolveTwoPointParallel4(
    p0grid_Tc[330],
    Tself,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    3200,
    0.1,
    u0 = [-1e-8, 6 * 8.0],
    config = config_spec,
)

typeof(solmassfun)



lines(p0grid_Tc,imdata_twopoint)


imdata_twopoint[250:255]
lines(p0grid_Tc,redata_twopoint)

redata_twopoint=hcat(tsolTwoPointparallelTc...)'[:,2]
imdata_twopoint=hcat(tsolTwoPointparallelTc...)'[:,1]


lines(
    p0grid_Tc,
    2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2),
    axis = (xscale = log,yscale = log,),
)


current_figure()



config_spec
Tlpa
@time tsolTwoPointparallel1 = pmap(v1) do p0
    tchanelSolveTwoPointParallel3(
        p0,
        Tlpa,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3_ac,
        Relambdaflow3_ac,
        800,
        0.2,
        u0 = [-0.00001, 6 * 8.0],
        config = config_spec,
    )
end



tsolTwoPointparallel3=tchanelSolveTwoPointParallel4(
    0.0,
    Tself,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    8000,
    0.01,
    u0 = [-0.0000001, 6 * 8.0],
    config = config_spec,
)




selfkrang=RGscale(800.0, 37.999)


lines(tsolTwoPointparallel3.t,-vcat(tsolTwoPointparallel3.u...))
lines!(sol1.t,hcat(sol1.u...)[1,:])

current_figure()

sqrt(800.0^2+solmassfun(800.0))
lines(v1,hcat(tsolTwoPointparallel2...)'[:,1],axis = (yscale = log,))
lines!(v1,-hcat(tsolTwoPointparallel1...)'[:,1],axis = (yscale = log,))

lines(v1,hcat(tsolTwoPointparallel2...)'[:,2])
lines!(v1,-hcat(tsolTwoPointparallel1...)'[:,2])

lines(v1,map(x->-2x[1]/(x[1]^2+x[2]^2),tsolTwoPointparallel2),axis = (yscale = log,))
lines!(v1,map(x->-2x[1]/(x[1]^2+x[2]^2),tsolTwoPointparallel1),axis = (yscale = log,))
current_figure()

config_spec

tsolTwoPointparallel2 = [[0.0, 0.0] for i = 1:400]


p0range = 1.0:1.0:400.0


@time Threads.@threads for i = 1:64
    tsolTwoPointparallel2[i] = tchanelSolveTwoPointParallel3(
        p0range[i],
        150.0, lpakrang,
        x -> sol1(x)[1],
        x -> sol1(x)[2],
        Imlambdaflow3,
        Relambdaflow3,
        800,
        0.2,
        u0 = [-0.00001, 6 * 8.0],)
end

@everywhere config_spec.dtmax=50.0
@everywhere config_spec.progress=false
@everywhere config_spec.atol=1e-6
@everywhere config_spec.rtol=1e-6
@time fourpointdatamk = tchanelSolveFourPointParallel(
    100.5,
    150.0,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3_ac,
    Relambdaflow3_ac,
    800,
    0.1,
    u0 = [-0.00001, 6 * 8.0],
    config = config_spec,
)

@time fourpointdatamk2 = tchanelSolveFourPointParallel(
    5.0,
    Tlpa,
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    6400,
    ϵ5,
    u0 = [-1e-15, 6 * 8.0],
    config = config_spec,
)


temf(p0,k0)=hcubature(
    k ->
        -48/3 * (
            TF.ReFb2(
                p0 + Eb(k0, tempsolmassfun(k0)),
                k[1],
                tempsolmassfun(k[1]),
                145.0,
                ϵ5,
            ) + TF.ReFb2(
                p0 - Eb(k0, tempsolmassfun(k0)),
                k[1],
                tempsolmassfun(k[1]),
                145.0,
                ϵ5,
            )
        ),
    [k0],
    [800.0],
    rtol=1e-6,
    atol=1e-6,
)[1]


plot(p0->temf(p0,400.0),10.0,400.0)


plot(fourpointdatamk2[1][10:2600], -fourpointdatamk2[2][10:2600], yaxis = :log)
plot!(fourpointdatamk2[1][10:3000], fourpointdatamk2[3][10:3000])


plot(fourpointdatamk2[1][10:1000], fourpointdatamk2[3][10:1000])


hcat(fourpointdatamk2...)


writedlm(
    joinpath(path, "Data/FourPoint/fourpointredataT=145q0k=100.dat"),
    hcat(fourpointdatamk2...),
)



tempsolmassfun(0.1)|>sqrt

config_spec

klistfourpoint=SharedArray(0.00001:0.5:800.00001 |>collect)

@time fourpointdatak0p0 = pmap(klistfourpoint) do k0
    tchanelSolveFourPointParallel(
        k0,
        Tlpa,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        800,
        0.001,
        u0 = [-0.000001, 6 * 8.0],
        config = config_spec,
    )
end



testsol1=tchanelSolveFourPointParallelSimple(
    0.0,
    Tlpa,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    0.01,
    u0 = [-0.0000001, 6 * 8.0],
    config = config_spec,
)

lines!(fourpointdatamk[1],-fourpointdatamk[2],axis = (yscale = log,))
lines(fourpointdatamk2[1],-fourpointdatamk2[2],axis = (yscale = log,))
current_figure()
lines(fourpointdatamk[1],-fourpointdatamk[3])
lines(fourpointdatamk2[1],-fourpointdatamk2[3])

lines(klistfourpoint,[fourpointdatak0p0[x][3][400] for x in 1:length(klistfourpoint)])

fourpointdatak0p0

[fourpointdatak0p0[x][3][400] for x in 1:length(klistfourpoint)]

redata_fourpoint = hcat([fourpointdatak0p0[x][3] for x = 1:length(klistfourpoint)]...)

imdata_fourpoint = hcat([fourpointdatak0p0[x][2] for x = 1:length(klistfourpoint)]...)
writedlm(
    joinpath(path, "Data/FourPoint/fourpointredataT=145q0=0.dat"),
    redata_fourpoint,
)
writedlm(
    joinpath(path, "Data/FourPoint/fourpointimdataT=145q0=0.dat"),
    imdata_fourpoint,
)
writedlm(
    joinpath(path, "Data/FourPoint/fourpointp0T=150.dat"),
    fourpointdatak0p0[1][1],
)
writedlm(
    joinpath(path, "Data/FourPoint/fourpointkT=150.dat"),
    klistfourpoint,
)

writedlm(
    joinpath(path, "Data/FourPoint/klpaT=145.dat"),
    sol1.t,
)

writedlm(
    joinpath(path, "Data/FourPoint/m2lpaT=145.dat"),
    hcat(sol1.u...)[1,:],
)

writedlm(
    joinpath(path, "Data/FourPoint/lambdalpaT=145.dat"),
    hcat(sol1.u...)[2,:],
)

lines(testsol1.t,hcat(testsol1.u...)[2,:])

lines!(sol1.t,hcat(sol1.u...)[2,:])
current_figure()


writedlm(joinpath(path,"Data/k=100mk.dat"),hcat(fourpointdatamk...))


config_spec.maxiters=10^10


plot(fourpointdata[1],log.(10,abs.(-fourpointdata[2])),label =L"\mathrm{Im}\lambda",title="k=200MeV",xlabel=L"p_0")

testfun1=Spline1D(fourpointdata[1],log.(10,abs.(-fourpointdata[2])))
q0=Epi(5.0,sol1(5.0)[1])

plot!([q0],[testfun1(q0)],seriestype = :scatter,label = L"q_0")

plot!([3q0],[testfun1(3q0)],seriestype = :scatter,label = L"3q_0")



writedlm("/home/wjfu1/tyy/O(4)Data/vertex/full/k=200T=150.dat",hcat(fourpointdata...))

writedlm("/home/wjfu1/tyy/O(4)Data/vertex/full/k=200.cofig",config_spec)

config_spec

config_spec.maxiters=10^10

# yaxis = :log,
plot!(fourpointdata[1],-fourpointdata[2],label =L"\mathrm{Im}\lambda",title="k=200MeV",xlabel=L"p_0",yaxis = :log)

plot(p0->TF.ImFb2(p0 - 10.0, 10.0, sol1(10.0)[1], 150.0, 0.1),1.0,100.0)


testfun1=Spline1D(fourpointdata[1],fourpointdata[2])
q0=Epi(200.0,sol1(200.0)[1])

plot!([q0],[-testfun1(q0)],yaxis = :log,seriestype = :scatter,label = L"q_0")

plot!([3q0],[-testfun1(3q0)],yaxis = :log,seriestype = :scatter,label = L"3q_0")


tsolTwoPointparallel2

plot(hcat(tsolTwoPointparallel2...)'[:, 1])


imvec = hcat((tsolTwoPointparallel2)...)'[:, 1]
revec = hcat((tsolTwoPointparallel2)...)'[:, 2]

lambda_omega = range(0.0, stop = 400.0, length = 100)


plot(p0range, revec)

plot!(lambda_omega, abs.(revec), yaxis = :log)


plot(p0range, -imvec ./ (imvec .^ 2 + revec .^ 2), yaxis = :log)
tsolparallel1[3](40.0)

plot(k -> sol1(k)[1], 1.0, 800.0)


plot(k -> -tsol1(k)[1], 1, 1000.0)


plot(k -> tsol1(k)[2], 1.0, 800.0, yaxis = :log)
