function VImintqs_q0(p0, ps,q0, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -hquadrature(
        kprim -> FRGRealTime.dkVImintqs(
            p0,
            ps,
            q0,
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim),
        ),
        k,
        UVScale,
        rtol = 1e-2,
        atol = 1e-2,
    )[1]
end

using HCubature,FiniteDifferences
plot(
    q0 -> VImintqs_q0(
        0.01,
        2.0,
        q0,
        10.0,
        145.0,
        4.0,
        1.0,
        800.0,
        m2fun,
        lamdafun,
    ),-50.0,50.0,title="V(q0)",xaxis="q0"
)

plot(
    q0 -> central_fdm(5, 1)(x->VImintqs_q0(
        10.0,
        2.0,
        x,
        10.0,
        145.0,
        4.0,
        1.0,
        800.0,
        m2fun,
        lamdafun,
    ),q0),-50.0,50.0,title="V'(q0)",xaxis="q0"
)
