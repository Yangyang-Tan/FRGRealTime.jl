function Coeffgamm2_dVdq0(k, T, mfun)
    (
        k * (coth(Epi(k, mfun(k)) / (2 * T)) / Epi(k, mfun(k))^2)
    ) / (16 * pi^2)
end

function VImintqs_delta1_dq0(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -(2 + Npi) *
    π *
    3 *
    (
        deltasumkAll_dp0(p0 + Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun) -
        deltasumkAll_dp0(p0 - Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun)
    )
end



function VImintqs_delta1_dVdq0(
    p0,
    ps,
    k,
    T,
    Npi,
    IRScale,
    UVScale,
    mfun,
    lamfun,
)
    return -(2 + Npi) *
           π *
           3 *
           central_fdm(5, 1)(
               z -> (
                   deltasumkAllSmooth(
                       p0 + z,
                       ps,
                       k,
                       T,
                       Npi,
                       IRScale,
                       UVScale,
                       mfun,
                       lamfun,
                   ) + deltasumkAllSmooth(
                       p0 - z,
                       ps,
                       k,
                       T,
                       Npi,
                       IRScale,
                       UVScale,
                       mfun,
                       lamfun,
                   )
               ),
               Epi(k, mfun(k)),
           )
end



function fastpropImintqs_dVdq0(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return -hcubature(
        x ->
            (x[1] - UVScale) *
            2 *
            central_fdm(5, 1)(
                z -> dkVImintqs(
                    p0,
                    ps,
                    z,
                    x[1],
                    x[2] * (x[1] - UVScale) + UVScale,
                    mfun(x[2] * (x[1] - UVScale) + UVScale),
                    T,
                    Npi,
                    lamfun(x[2] * (x[1] - UVScale) + UVScale),
                ),
                Epi(x[1], mfun(x[1])),
            ) *
            Coeffgamm2_dVdq0(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end


function propImintqs_delta1_dVdq0(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return -hquadrature(
        k ->
            2 *
            VImintqs_delta1_dVdq0(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun) *
            Coeffgamm2_dVdq0(k, T, mfun),
        IRScale,
        UVScale,
        rtol = 1e-6,
        atol = 1e-6,
        initdiv=20,
        maxevals = 8000,
    )[1]
end


function fastpropImintqs_All_dVdq0(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return fastpropImintqs_dVdq0(
        p0,
        ps,
        T,
        IRScale,
        UVScale,
        Npi,
        mfun,
        lamfun;
        kwarg...,
    ) + propImintqs_delta1_dVdq0(
        p0,
        ps,
        T,
        IRScale,
        UVScale,
        Npi,
        mfun,
        lamfun,
    )
end







@doc raw"""
    propImintqs_delta1_dq0_delta(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)

Compute the delta function part appears
"""
function propImintqs_delta1_dq0_delta1(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun,
)
    deltaf1(x) = Epi(x, mfun(x)) - p0
    if deltaf1(IRScale) * deltaf1(UVScale) >= 0
        return 0.0
    else
        k01 = find_zero(deltaf1, (IRScale, UVScale))
        return (2 + Npi) *sqrt(k01^2+mfun(k01))/(2*k01+derivative(mfun,k01))*
               π *
               3 *
               2 *deltasumkfix_dp0_deltafun(p0+Epi(k01,mfun(k01)),
                   ps,
                   k01,
                   k01,
                   T,
                   Npi,
                   IRScale,
                   UVScale,
                   mfun,
                   lamfun,
               ) * Coeffgamm2_dVdq0(k01, T, mfun)
    end
end

function propImintqs_delta1_dq0_delta2(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun,
)
    deltaf2(x) = 3 * Epi(x, mfun(x)) - p0
    if deltaf2(IRScale) * deltaf2(UVScale) >= 0
        return 0.0
    else
        k02 = find_zero(deltaf2, (IRScale, UVScale))
        return -(2 + Npi) *sqrt(k02^2+mfun(k02))/(2*k02+derivative(mfun,k02))*
               π *
               3 *
               2 *deltasumkfix_dp0_deltafun(p0-Epi(k02,mfun(k02)),
                   ps,
                   k02,
                   k02,
                   T,
                   Npi,
                   IRScale,
                   UVScale,
                   mfun,
                   lamfun,
               ) * Coeffgamm2_dVdq0(k02, T, mfun)
    end
end

function propImintqs_delta1_dq0(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return -hquadrature(
        k ->
            2 *
            VImintqs_delta1_dq0(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun) *
            Coeffgamm2_dVdq0(k, T, mfun),
        IRScale,
        UVScale,
        rtol = 1e-6,
        atol = 1e-6,
        maxevals = 2000,
    )[1]
end

function fastpropImintqs_All_dq0(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return fastpropImintqs_dVdq0(
               p0,
               ps,
               T,
               IRScale,
               UVScale,
               Npi,
               mfun,
               lamfun;
               kwarg...,
           ) +
           propImintqs_delta1_dq0(
               p0,
               ps,
               T,
               IRScale,
               UVScale,
               Npi,
               mfun,
               lamfun,
           ) +
           propImintqs_delta1_dq0_delta1(
               p0,
               ps,
               T,
               IRScale,
               UVScale,
               Npi,
               mfun,
               lamfun,
           )+propImintqs_delta1_dq0_delta2(
               p0,
               ps,
               T,
               IRScale,
               UVScale,
               Npi,
               mfun,
               lamfun,
           )
end
