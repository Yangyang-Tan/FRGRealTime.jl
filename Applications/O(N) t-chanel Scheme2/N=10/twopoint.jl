solmassfun(0.0)
p0grid_T=SharedArray(0.0:0.5:400.0|>collect)
@time tsolTwoPointparallelT = pmap(p0grid_T) do p0
    tchanelSolveTwoPointParallel4(
        p0,
        Tlpa,
        lpakrang,
        solmassfun,
        sollambdafun,
        Imlambdaflow3,
        Relambdaflow3,
        5*8000,
        ϵ5,
        u0 = [-1e-8, lambdaα*6 * 8.0],
        Nf=10.0,
        config = config_spec,
    )
end
@everywhere config_spec.atol=1e-6
@everywhere config_spec.rtol=1e-6

redata_twopoint=(hcat(tsolTwoPointparallelT...)'[:,2])
imdata_twopoint=(hcat(tsolTwoPointparallelT...)'[:,1])
plot(
    p0grid_T,
    2 * imdata_twopoint ./ (imdata_twopoint .^ 2 + redata_twopoint .^ 2),yaxis=:log
)


current_figure()


Tself

writedlm(
    joinpath(path, "N=10/twopointredataT=145_N=10.dat"),
    redata_twopoint,
)
writedlm(
    joinpath(path, "N=10/twopointimdataT=145_N=10.dat"),
    imdata_twopoint,
)

writedlm(
    joinpath(path, "N=10/twopointp0_T_N=10.dat"),
    p0grid_T,
)

tchanelSolveTwoPointParallel4(
    0.0,
    Tlpa,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    5*8000,
    ϵ5,
    u0 = [-1e-8, lambdaα*6 * 8.0],
    Nf=2.0,
    config = config_spec,
)
