solselfGPU_log = tchanelSolveFourPointGPU(
    # 9.4685f0,
    130.835,
    kmin = 1.0,
    ϵ = 5.0,
    pgridmax = 200.0,
    dp0 = 0.125 / 0.125,
    dq0 = 0.125 / 0.125,
    Nf = 4.0,
    config = config_spec,
)
heatmap(solselfGPU_log.sol[:,2:end,2])
heatmap(solselfGPU_log.sol[:,2:end,1])

solselfGPU_log = tchanelSolveFourPointGPU2(
    # 9.4685f0,
    145.835,
    kmin = 200.0,
    ϵ = 5.0,
    pgridmax = 100.0,
    dp0 = 0.125,
    dq0 = 0.125,
    Nf = 4.0,
    config = config_spec,
)
1
