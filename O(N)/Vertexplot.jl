plot(
    p0 ->
        FRGRealTime.VImintqs(
            p0,
            0.0001,
            80.0,
            145.0,
            4.0,
            1.0,
            800.0,
            m2fun,
            lamfun,
        ) + FRGRealTime.VImintqs_delta1(
            p0,
            0.0001,
            80.0,
            145.0,
            4.0,
            1.0,
            800.0,
            m2fun,
            lamfun,
        ),
    1.0,
    800.0,
    label = "I2, p0-q0",
    legend = :topleft,
)


quadgk(k->FRGRealTime.dkV4piResimple(
    0.01,
    0.0001,
    Epi(100.0, m2fun(100.0)),
    100.0,
    0.0,
    k,
    m2fun(k),
    145.0,
    lamfun(k),
    4.0,
),800.0,100.0)[1]


plot(qs->quadgk(k->FRGRealTime.dkV4piResimple(
    0.000001,
    0.0001,
    0.0001,
    qs,
    0.0,
    k,
    m2fun(k),
    145.0,
    lamfun(k),
    4.0,
),800.0,100.0)[1]-64,1.0,100.0)


plot(
    k -> FRGRealTime.dkF1Tilde3(k,m2fun(k), 145.0),
    2.0,
    100.0,
)


plot(
    k -> central_fdm(5, 1)(
        qsmax ->
            FRGRealTime.PvdkF1Tildeps(0.001, qsmax, k, m2fun(k), 145.0),
        2.0,
    ),
    0.01,
    10.0,
)




v1 = collect(1.0:0.2:800.0)


v2 = exp.(collect(log(1e-8):(log(1e-4)-log(1e-8))/1500:log(1e-4)))

outv1 = map(v1) do p0
    FRGRealTime.VImintqs(
        p0,
        0.0001,
        k_temp,
        145.0,
        4.0,
        1.0,
        800.0,
        m2fun,
        lamfun,
    ) + FRGRealTime.VImintqs_delta1(
        p0,
        0.0001,
        k_temp,
        145.0,
        4.0,
        1.0,
        800.0,
        m2fun,
        lamfun,
    )
end

FRGRealTime.loopfunppfix.(v2, 1e-6, 2.0, 2.0, 10.0)

plot(
    p0 -> FRGRealTime.loopfunppfix(p0, 1e-6, 2.0, 2.0, 10.0),
    1e-8,
    1e-4,
    xaxis = :log,
)

writedlm("$dir/k=50.dat", 3 * outv1)

plot!(v1, outv1 / 100^3)
dir = "/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/Loop"
k_temp = 50
writedlm("$dir/looppp_m=2.dat", FRGRealTime.loopfunppfix.(v2, 1e-6, 2.0, 2.0, 10.0))

k_temp = [50, 80, 100, 200, 300, 400, 500, 600, 700]
writedlm("$dir/klist.dat", k_temp)
writedlm("$dir/q0list.dat", Epi.(k_temp, m2fun.(k_temp)))

Epi.(k_temp, m2fun.(k_temp))



Epi(200, m2fun(200))
