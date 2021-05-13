################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################


@doc raw"""
    dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)

compute $\int_0^{qsmax}dq_s qs^2\int_{-1}^{1}d\cos\theta \tilde{\partial_k}\mathrm{Im}V(q_0)$.

`dkV4piImintqs` only contains $V(q_0)$, for $-q_0$, we have $\int d\cos\theta V(q_0)=\int d\cos\theta V(-q_0)$,
so we need an extra $2$ at somewhere.

# Arguments
- `qsmax`: we integrate $q_s$ from $0$ to $k$, `qsmax` will set to `k` when we do the integration $dk'$, it should be distinguished from $k'$
- `m`: mass square, it will be $m(k')$ when we do the integration $dk'$.
- `lam4pik`: $\lambda_{4\pi}$, it will be $\lambda_{4\pi}(k')$ when we do the integration $dk'$ .
"""
function dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    π *
    3 *
    (
        # dkF1Allintqs(p0 - q0, ps, qsmax, k, m, T) +
        # dkF1Allintqs(p0 + q0, ps, qsmax, k, m, T) +
        dkF2Allintqs(p0 - q0, ps, qsmax, k, m, T) +
        dkF2Allintqs(p0 + q0, ps, qsmax, k, m, T)
    )
end


@doc raw"""
    VImintqs(p0, ps, k, T, Npi,UVScale,mfun::Function,lamfun::Function)

compute $\int_0^{k}dq_s qs^2\int_{-1}^{1}d\cos\theta \mathrm{Im}V(q_0,k)$.
In our code, we perform integration over `kprim`, `q0` & `qs` does not involved,
so `qs=k`, `q0=Epi(k, mfun(k))`.

# Arguments
- `mfun::Function`: $m^2(k)$, input from zero momentum result
- `lampifun::Function`: $\lambda_{4\pi}(k)$, input from zero momentum result.
"""
function VImintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -hquadrature(
        kprim -> dkVImintqs(
            p0,
            ps,
            Epi(k, mfun(k)),
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim),
        ),
        k,
        UVScale,
        rtol = 1e-5,
        atol = 1e-5,
        maxevals = 1000,
    )[1]
    # +
    # (2 + Npi) *
    # π *
    # 3 *
    # (
    #     deltasumkAll(p0 + Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun) +
    #     deltasumkAll(p0 - Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun)
    # )
end



function Coeffgamm2(k, T, Npi, mfun)
    (
        k * (
            -(coth(Epi(k, mfun(k)) / (2 * T)) / Epi(k, mfun(k))^3) -
            csch(Epi(k, mfun(k)) / (2 * T))^2 / (2 * T * Epi(k, mfun(k))^2)
        )
    ) / (16 * pi^2)
end


function propImsimpleintqs(p0, ps, T,IRScale,UVScale, Npi, mfun, lamfun)
    -hquadrature(
        k ->
            2 *
            VImintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun) *
            Coeffgamm2(k, T, Npi, mfun),
        IRScale,
        UVScale,
        rtol = 1e-5,
        atol = 1e-5,
        maxevals = 1000,
    )[1]
end
