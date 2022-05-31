solmassfun(0.0)
p0grid_T=SharedArray(0.0:0.5:400.0|>collect)
@time tsolTwoPointparallelT = pmap(p0grid_T) do p0
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
config_spec.atol=1e-6
config_spec.rtol=1e-6

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
