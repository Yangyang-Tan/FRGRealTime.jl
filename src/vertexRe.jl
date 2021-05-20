#Zero Momentum
function lam4piflow(k, m2, T, Npi, lamda4pi)
    lamda4pi^2 * (8 + Npi) * (dkF1TildeAll(k, m2, T) + dkF2TildeAll(k, m2, T))
end
function dkVintqs_Zero(q0, k, kprim, m, T, lam4pik, Npi)
    lam4pik^2 *
    (2 + Npi) *
    (
        k^3 *
        (2 + Npi) *
        (dkF1TildeAll(kprim, m, T) + dkF2TildeAll(kprim, m, T)) / 3 +
        6 * PvdkF1Tildeps(q0, k, kprim, m, T) +
        6 * PvdkF2Tildeps(q0, k, kprim, m, T)
    )
end



function Vintqs_Zero(k, T, Npi, lamdapiΛ, hfun, UVScale)
    -2 * hquadrature(
        kprim -> dkVintqs_Zero(
            Epi(k, -hfun(k)[2]),
            k,
            kprim,
            -hfun(kprim)[2],
            T,
            hfun(kprim)[1],
            Npi,
        ),
        k,
        UVScale,
        rtol = error2,
        atol = error2,
        # initdiv=200,
    )[1] + (2 * k^3 * lamdapiΛ * (2 + Npi)) / 3
end


function Vintqs_LPA(k, Npi, lamdapik)
    (2 * k^3 * lamdapik * (2 + Npi)) / 3
end

function propReZeroflow(k, T, Npi, lamdapiΛ, hfun, UVScale)
    2 *
    Coeffgamm2Simple(k, T, Npi, -hfun(k)[2]) *
    Vintqs_Zero(k, T, Npi, lamdapiΛ, hfun, UVScale)
end
propReLPAflow(k, T, Npi, hfun) =
        2*
        Gamma2coeff(k, -hfun(k)[2],T) *
        Vintqs_LPA(k, Npi, hfun(k)[1])



################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################
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
    lam4pik^2 *
    (2 + Npi) *
    (
        3 * (
            dkF1TildeintqsAll(p0 - q0, ps, qsmax, k, m, T; kwarg...) +
            dkF1TildeintqsAll(p0 + q0, ps, qsmax, k, m, T; kwarg...) +
            dkF2TildeintqsAll(p0 - q0, ps, qsmax, k, m, T; kwarg...) +
            dkF2TildeintqsAll(p0 + q0, ps, qsmax, k, m, T; kwarg...)
        ) +
        (Npi + 2) * 2 / 3 *
        qsmax^3 *
        (dkF1TildeAll(k, m, T) + dkF2TildeAll(k, m, T))
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
        UVScale,
        rtol = 1e-3,
        atol = 1e-3,
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
            Coeffgamm2(k, T, Npi, mfun),
        IRScale,
        UVScale,
        rtol = 1e-3,
        atol = 1e-3,
    )[1] - p0^2 + ps^2 - mfun(UVScale)
end
