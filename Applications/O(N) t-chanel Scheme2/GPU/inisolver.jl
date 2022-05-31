iniu0=tchanelSolveFourPointGPU_ini(
    145.0,
    maxiters=20,
    dp0 = 0.5,
    dq0 = 0.5,
    ϵ = 100.0,
    Nf = 4.0,
    pgridmax = 800.0,
)
heatmap(Array(iniu0[:,2:end,1]))
heatmap(Array(iniu0[:,2:end,1]))

Array(iniu0[:,2:end,1])
isnan.(Array(iniu0[:,2:end,1]))|>sum
TF.ReFb2(0.0, 800.0,-0.18216546875 * lpakrang.UV^2 , 145.0, 1.0)
