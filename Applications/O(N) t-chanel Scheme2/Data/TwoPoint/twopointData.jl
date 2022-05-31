@everywhere include(joinpath(path,"ini.jl"))
@everywhere Tself=Tself+0.101
@everywhere include(joinpath(path,"Zero.jl"))


testv1=collect(1:100.)
testv2=similar(testv1)

@tullio testv2[i]:=testv1[i]



@time solself = solve(
    probself,
    Tsit5(),
    dense = false,
    save_on = false,
    save_start = false,
    save_timeseries=false,
    save_end = true,
    reltol = 1e-6,
    abstol = 1e-6,
    progress = true,
    dtmax=50.0,
    dtmin=1e-14,
)
# solself=1
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
solmassfun(0.0)
p0grid_T=SharedArray(0.0:0.5:400.0|>collect)
@time tsolTwoPointparallelTc = pmap(p0grid_T) do p0
    tchanelSolveTwoPointParallel4(
        p0,
        Tself,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        20*8000,
        ϵ5,
        u0 = [-1e-8, 6 * 8.0],
        config = config_spec,
    )
end

@everywhere config_spec.atol=1e-6
@everywhere config_spec.rtol=1e-6

redata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,2])
imdata_twopoint=(hcat(tsolTwoPointparallelTc...)'[:,1])
plot!(
    p0grid_T,
    2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2),yaxis=:log
)


current_figure()


Tself

writedlm(
    joinpath(path, "Data/TwoPoint/twopointredataT=160.dat"),
    redata_twopoint,
)
writedlm(
    joinpath(path, "Data/TwoPoint/twopointimdataT=160.dat"),
    imdata_twopoint,
)

writedlm(
    joinpath(path, "Data/TwoPoint/twopointp0_T.dat"),
    p0grid_T,
)
Tlpa-Tself-8.07
