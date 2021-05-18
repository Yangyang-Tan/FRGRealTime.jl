################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################






@doc raw"""
    dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)

compute $\int_0^{qsmax}dq_s qs^2\int_{-1}^{1}d\cos\theta \tilde{\partial_k}\mathrm{Im}V(q_0)$.

`dkVImintqs` only contains $V(q_0)$, for $-q_0$, we have $\int d\cos\theta V(q_0)=\int d\cos\theta V(-q_0)$,
so we need an extra $2$ at somewhere.

`dkVImintqs` doesn't have any delta function contribution, we include the type-1 delta function in `VImintqs_delta1` separately.

# Arguments
- `qsmax`: we integrate $q_s$ from $0$ to $k$, `qsmax` will set to `k` when we do the integration $dk'$, it should be distinguished from $k'$
- `m`: mass square, it will be $m(k')$ when we do the integration $dk'$.
- `lam4pik`: $\lambda_{4\pi}$, it will be $\lambda_{4\pi}(k')$ when we do the integration $dk'$ .
"""
function dkVReintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    (
        3 * (
            dkF1Allintqs(p0 - q0, ps, qsmax, k, m, T) +
            dkF1Allintqs(p0 + q0, ps, qsmax, k, m, T) +
            dkF2Allintqs(p0 - q0, ps, qsmax, k, m, T) +
            dkF2Allintqs(p0 + q0, ps, qsmax, k, m, T)
        ) +
        (Npi + 2) *2/3 *qsmax^3*(
            dkF1All(1e-8 - 1e-14, 1e-8, k, m, T) +
            dkF2All(1e-8 - 1e-14, 1e-8, k, m, T)
        )
    )
end





@doc raw"""
    VImintqs(p0, ps, k, T, Npi,UVScale,mfun::Function,lamfun::Function)

compute $\int_0^{k}dq_s qs^2\int_{-1}^{1}d\cos\theta \mathrm{Im}V(q_0,k)$.
In our code, we perform integration over `kprim`, `q0` & `qs` does not involved,
so `qs=k`, `q0=Epi(k, mfun(k))`.


`VImintqs` don't have any delta function contribution, we include the type-1 delta function in `VImintqs_delta1` separately.

# Arguments
- `mfun::Function`: $m^2(k)$, input from zero momentum result
- `lampifun::Function`: $\lambda_{4\pi}(k)$, input from zero momentum result.
"""
function VReintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
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
        rtol = 1e-4,
        atol = 1e-4,
    )[1]
end




function propReintqs(p0, ps, T,IRScale,UVScale, Npi, mfun, lamfun)
    -hquadrature(
        k ->
            2 *
            VImintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun) *
            Coeffgamm2(k, T, Npi, mfun),
        IRScale,
        UVScale,
        rtol = 1e-4,
        atol = 1e-4,
    )[1]
end
