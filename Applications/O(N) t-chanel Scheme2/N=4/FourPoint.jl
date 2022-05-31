@everywhere config_spec=deepcopy(config_self)
@everywhere config_spec.atol=1e-8
@everywhere config_spec.rtol=1e-8



solselftinput=reverse(deepcopy(solselft))
itp = interpolate(solselftinput, BSpline(Linear()))
@eval @everywhere itp = $itp
@eval @everywhere solselftinput = $solselftinput

@time fourpointdatamk = tchanelSolveFourPointParallel(
    15.0,
    Tlpa,
    lpakrang,
    tempsolmassfun,
    tempsollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    6400,
    ϵ5,
    u0 = [-1e-15, lambdaα * lambda_UV],
    Nf=-4,
    config = config_spec,
)
plot!(fourpointdatamk[1][20:end], -fourpointdatamk[2][20:end], yaxis = :log)
plot!(fourpointdatamk[1], fourpointdatamk[3])
plot(solselft,solselfu1)
