using FRGRealTime
using HCubature
using Plots


δk=0.02
ppfunps_delta(p0, k, kprim, m, T)

plot(p0 -> 2*ppfunps_delta(p0, 10.0, 10.0, -2.0, 145.0), 19.7, 19.9)
plot(
    p0 -> FRGRealTime.delta2_intcosthqs(p0, 0.001, 10.0, 10.0, -2.0, 145.0),
    19.7, 19.9
)





plot(
    p0 -> FRGRealTime.F4tildekSmooth_intcosthqs(p0,0.001,10.0,-2.0,145.0,400.0,0.1,0.00001),
    19.0, 20.9
)


function testf(p0, k)
    -hquadrature(
        x -> FRGRealTime.delta2Smooth_intcosthqs(p0, 0.001, k, x, -2.0, 145.0,0.00,0.01),
        k,
        400.0,initdiv=500,rtol=1e-8,atol=1e-8,
    )[1]
end




plot!(p0->testf(p0,10.0),19.0, 20.9)



plot(
    p0 -> FRGRealTime.delta2_intcosthqs(p0, 0.001, 10.0, 10.0, -2.0, 145.0),
    1.0, 30.9
)


plot(
    p0 -> FRGRealTime.PvdkF1Tildeps(p0, 0.001, 10.0, 10.0, -2.0, 145.0),
    30.0,
    40.0,
)



plot(
    p0 -> FRGRealTime.PvdkF1Tildeps(p0, 10.0, 10.0, -2.0, 145.0),
    30.0,
    40.0,
)



plot!(
    p0 -> FRGRealTime.PvdkF1Tildeps2(p0, 10.0, 10.0, -2.0, 145.0),
    30.0,
    40.0,
)



plot(
    p0 ->
        FRGRealTime.flowpp_intcostheqs(p0, 0.0000000, 10.0, 10.0, -50.0, 145.0),
    2 * Epi(10.0, -50),
    20.0,
)

plot!(
    p0 -> 2 * FRGRealTime.ppfunps(p0, 10.0, 10.0, -50.0, 145.0),
    2 * Epi(10.0, -50),
    20.0,
)


testf()=2 * Epi(k, m)

function flowppintp0qs(p0, qsmax, k, m, T)
    quadgk(
        ps -> flowppintp0(p0, ps, k, m, T),
        0.0,
        qsmax,
        atol = 1e-4,
        rtol = 1e-4,
    )
end


plot(p0->FRGRealTime.flowpm(p0, p0, 10.0, -50.0, 100.0),0.0000001,0.1)


plot(k->FRGRealTime.loopfunpm(0.1, 0.1, k, -10.0, 100.0),4.0,20.0)


plot(p0->FRGRealTime.loopfunpm(p0, p0+0.1, 4.0, -10.0, 100.0),0.0,0.00000001)

using Plots,FiniteDifferences

plot(p0->central_fdm(5,1)(x->FRGRealTime.deltasumkfix(
    x,
    130.001,
    10.0,
    145.0,
    400.0,
    1.0,
    800.0,
    m2fun,
    lamfun,
),p0),1.0,800.0)


plot(p0->FRGRealTime.deltasumkfix(
    p0,
    0.001,
    80.0,
    145.0,
    400.0,
    1.0,
    800.0,
    m2fun,
    lamfun,
),250.0,260.0)

using HCubature


hquadrature(
    x -> FRGRealTime.deltafun(x - 10.0,0.1),
    3.0,
    20.0,
    rtol = 1e-12,
    atol = 1e-12,
    initdiv = 20,
)[1]


hquadrature(k->hquadrature(
    x -> FRGRealTime.deltadxfun(x - 10.0,0.01),
    k,
    20.0)[1],1.0,20.0)[1]



plot(
    p0 -> FRGRealTime.deltasumkfix(
        p0,
        0.00001,
        1.0,
        145.0,
        4.0,
        1.0,
        150.0,
        m2fun,
        lamfun,
    ),38.0,42.0
)


plot(
    p0 -> FRGRealTime.VImintqs(
        p0,
        0.0001,
        50.0,
        100.0,
        4.0,
        1.0,
        200.0,
        m2funnew,
        lamfunnew,
    )+FRGRealTime.VImintqs_delta1(
        p0,
        0.0001,
        50.0,
        100.0,
        4.0,
        1.0,
        200.0,
        m2fun,
        lamfun,
    ),52.0,60.0
)


plot(p0->central_fdm(5,1)(x->FRGRealTime.VImintqs(
    p0,
    0.0001,
    x,
    50.0,
    100.0,
    4.0,
    1.0,
    200.0,
    m2fun,
    lamfun,
),Epi(50.0,100.0)),2.0,20.0)


plot!(
    p0 -> central_fdm(3,1)(x->FRGRealTime.VImintqsSimple(
        p0,
        0.0001,
        x,
        50.0,
        100.0,
        4.0,
        10.0^2,
        -10.0,
        200.0,
    ),Epi(50.0,100.0)),
    2.0,
    20.0,
)

using Plots
plot(
    p0 -> FRGRealTime.VImintqs(
        p0,
        0.0001,
        100.0,
        40.0,
        4.0,
        1.0,
        1600.0,
        m2fun,
        lamfun,
    ),
    1.0,
    800.0,label="I2, p0-q0",legend=:topleft,
)

sqrt(50^2+m2fun(50))
xaxis!("p0")

using QuadGK

function testf2(f;initdiv,kwarg...)
    quadgk(f,0.0,10.0;order=initdiv,kwarg...)
end
