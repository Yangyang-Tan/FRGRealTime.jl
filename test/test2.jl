FRGRealTime.ppfunps_delta(10.0, 10.0, 20.0, m2fun(20.0), 145.0)

plot(
    k -> FRGRealTime.dkF1TildeintqsAll_delta2(
        50.0,
        1.0,
        10.0,
        k,
        m2fun(k),
        145.0,
    ),
    10,
    50,
)

dkVReintqs_delta2(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
plot(
    k -> FRGRealTime.dkVReintqs_delta2(
        50.0,
        1.0,
        5.0,
        10.0,
        k,
        m2fun(k),
        145.0,
        4.0,
        lamfun(k),
    ),10,50
)

plot(
    k -> FRGRealTime.dkF1Tildeps_delta2(50.0, 10.0, k, m2fun(k), 145.0),
    10,
    50,
)


FRGRealTime.fastpropReintqs_delta2(
    50.0,
    1.0,
    145.0,
    1.0,
    800.0,
    4.0,
    m2fun,
    lamfun;maxevals=2000000)
