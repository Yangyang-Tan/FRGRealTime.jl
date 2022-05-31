# Zero Momentum
function lam4piflow(k, m2, T, Npi, lamda4pi)
    return lamda4pi^2 * (8 + Npi) * (dkF1TildeAll(k, m2, T) 
    + dkF2TildeAll(k, m2, T)
    )
end
function dkVintqs_Zero(q0, k, kprim, m, T, lam4pik, Npi; kwargs...)
    return lam4pik^2 *(2 + Npi) *(
        k^3 *
        (2 + Npi) *
        (dkF1TildeAll(kprim, m, T) 
        + dkF2TildeAll(kprim, m, T)
        ) / 3 +
        6 * dkF1TildeintqsAll(q0, k, kprim, m, T; kwargs...) 
        +
        6 * dkF2TildeintqsAll(q0, k, kprim, m, T; kwargs...)
    )
end

function Vintqs_Zero(k, T, Npi, hfun, UVScale; kwargs...)
    -2 * hquadrature(
        kprim -> dkVintqs_Zero(
            Epi(k, -hfun(k)[2]),
            k,
            kprim,
            -hfun(kprim)[2],
            T,
            hfun(kprim)[1],
            Npi;
            kwargs...,
        ),
        k,
        UVScale;
        kwargs...,
        # initdiv=200,
    )[1] + (2 * k^3 * hfun(UVScale)[1] * (2 + Npi)) / 3
end

function Vintqs_LPA(k, Npi, lamdapik)
    (2 * k^3 * lamdapik * (2 + Npi)) / 3
end

function propReZeroflow(k, T, Npi, hfun, UVScale; kwargs...)
    2 *
    Coeffgamm2Simple(k, T, -hfun(k)[2]) *
    Vintqs_Zero(k, T, Npi, hfun, UVScale; kwargs...)
end

function propReLPAflow(k, T, Npi, hfun)
    2 * Coeffgamm2Simple(k, T, -hfun(k)[2]) * Vintqs_LPA(k, Npi, hfun(k)[1])
end


function dkV4piResimple(p0, ps, q0, qs, costhe, k, m, T, lam4pik, Npi)
    return lam4pik^2 *
           (2 + Npi) *
           (
               2 * (2 + Npi) * (dkF1TildeAll(k, m, T) + dkF2TildeAll(k, m, T)) +
               6 * dkF1TildeAll(
                   p0 - q0,
                   max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                   k,
                   m,
                   T,
               ) +
               6 * dkF1TildeAll(
                   p0 + q0,
                   max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                   k,
                   m,
                   T,
               ) +
               6 * (
                   dkF2TildeAll(
                       p0 - q0,
                       max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                       k,
                       m,
                       T,
                   ) + dkF2TildeAll(
                       p0 + q0,
                       max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                       k,
                       m,
                       T,
                   )
               )
           )
end




@doc raw"""
    dkVReintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)

compute $\int_0^{qsmax}dq_s qs^2\int_{-1}^{1}d\cos\theta \tilde{\partial_k}\mathrm{Re}V(q_0)$.

`dkVImintqs` only contains $V(q_0)$, for $-q_0$, we have $\int d\cos\theta V(q_0)=\int d\cos\theta V(-q_0)$,
so we need an extra $2$ at somewhere.

`dkVReintqs` contains type-1 and type-2 delta function

# Arguments
- `qsmax`: we integrate $q_s$ from $0$ to $k$, `qsmax` will set to `k` when we do the integration $dk'$, it should be distinguished from $k'$
- `m`: mass square, it will be $m(k')$ when we do the integration $dk'$.
- `lam4pik`: $\lambda_{4\pi}$, it will be $\lambda_{4\pi}(k')$ when we do the integration $dk'$ .
"""
function dkVReintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik; kwarg...)
    return lam4pik^2 * (2 + Npi) * (3 * (
            dkF1TildeintqsAll(p0 - q0, ps, qsmax, k, m, T; kwarg...) +
            dkF1TildeintqsAll(p0 + q0, ps, qsmax, k, m, T; kwarg...) 
            +
            dkF2TildeintqsAll(p0 - q0, ps, qsmax, k, m, T; kwarg...) +
            dkF2TildeintqsAll(p0 + q0, ps, qsmax, k, m, T; kwarg...)
        ) +
        (Npi + 2) * 2 / 3 *
        qsmax^3 *
        (dkF1TildeAll(k, m, T) 
        + dkF2TildeAll(k, m, T)
        )
    )
end


function dkVReintqs_Compensate(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    return lam4pik^2 *
           (2 + Npi) *
           (
               3 * (
                   dkF1TildeintqsAll_Compensate(p0 - q0, ps, qsmax, k, m, T) +
                   dkF1TildeintqsAll_Compensate(p0 + q0, ps, qsmax, k, m, T)
               )
           )
end


function dkVReintqs_delta1(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    return lam4pik^2 *
           (2 + Npi) *
           (
               3 * (
                   dkF1TildeintqsAll_delta1(p0 - q0, ps, qsmax, k, m, T) +
                   dkF1TildeintqsAll_delta1(p0 + q0, ps, qsmax, k, m, T)
               )
           )
end

function dkVReintqs_delta2(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    return lam4pik^2 *
           (2 + Npi) *
           (
               3 * (
                   dkF1TildeintqsAll_delta2(p0 - q0, ps, qsmax, k, m, T) +
                   dkF1TildeintqsAll_delta2(p0 + q0, ps, qsmax, k, m, T)
               )
           )
end


@doc raw"""
    VReintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)

compute $\int_0^{k}dq_s qs^2\int_{-1}^{1}d\cos\theta \mathrm{Re}V(q_0,k)$.
In our code, we perform integration over `kprim`, `q0` & `qs` does not involved,
so `qs=k`, `q0=Epi(k, mfun(k))`.


`VReintqs` contains type-1 and type-2 delta function.

# Arguments
- `mfun::Function`: $m^2(k)$, input from zero momentum result
- `lampifun::Function`: $\lambda_{4\pi}(k)$, input from zero momentum result.
"""
function VReintqs(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun; kwarg...)
    -hquadrature(
        kprim -> dkVReintqs(
            p0,
            ps,
            Epi(k, mfun(k)),
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim);
            kwarg...,
        ),
        k,
        UVScale;
        kwarg...,
    )[1] + 2 * (k^3 * lamfun(UVScale) * (2 + Npi)) / 3
end




function propReintqs(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun; kwarg...)
    -hquadrature(
        k ->
            2 *
            VReintqs(
                p0,
                ps,
                k,
                T,
                Npi,
                IRScale,
                UVScale,
                mfun,
                lamfun;
                kwarg...,
            ) *
            Coeffgamm2(k, T, mfun),
        IRScale,
        UVScale,
        rtol=1e-8,
        atol=1e-8,
    )[1] + p0^2 - ps^2 - mfun(UVScale)
end

function fastpropReintqs(
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
            dkVReintqs(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale);
                kwarg...,
            ) *
            Coeffgamm2(x[1], T, mfun) +
            2 *
            Coeffgamm2(x[1], T, mfun) *
            (2 * x[1]^3 * lamfun(UVScale) * (2 + Npi)) / 3,
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1] + p0^2 - ps^2 - mfun(UVScale)
end



function fastpropReintqs_Compensate(
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
            dkVReintqs_Compensate(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end




function fastpropReintqs_delta1(
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
            dkVReintqs_delta1(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end


function fastpropReintqs_delta2(
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
            dkVReintqs_delta2(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end



function fastpropReintqs_All(
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
    return fastpropReintqs(
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
           fastpropReintqs_Compensate(
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
           fastpropReintqs_delta1(
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
           fastpropReintqs_delta2(
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
end



function fastpropReintqs2(
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
            dkVReintqs(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale);
                kwarg...,
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1] - hquadrature(
        k ->
            2 *
            Coeffgamm2(k, T, mfun) *
            (2 * (k^3 * lamfun(UVScale) * (2 + Npi)) / 3),
        IRScale,
        UVScale;
        kwarg...,
    )[1] + p0^2 - ps^2 - mfun(UVScale)
end
