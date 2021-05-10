################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################


@doc raw"""
    dkV4piImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)

compute $\int_0^{qsmax}dq_s qs^2\int_{-1}^{1}d\cos\theta \tilde{\partial_k}\mathrm{Im}V(q_0)$.

`dkV4piImintqs` only contains $V(q_0)$, $\int d\cos\theta V(q_0)=V(q_0)$, $\int d\cos\theta V(-q_0)$,
so we need an extra $2$ at somewhere.

# Arguments
- `qsmax`: we integrate $q_s$ from $0$ to $k$, `qsmax` will set to `k` when we do the integration $dk'$, it should be distinguished from $k'$ 
- `m`: mass square, it will be $m(k')$ when we do the integration $dk'$.
- `lam4pik`: $\lambda_{4\pi}$, it will be $\lambda_{4\pi}(k')$ when we do the integration $dk'$ .
"""
dkV4piImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik) =
    lam4pik^2 *
    (2 + Npi) *
    π *
    3 *
    (
        dkF1Allintqs(p0 - q0, ps, qsmax, k, m, T) +
        dkF1Allintqs(p0 + q0, ps, qsmax, k, m, T) +
        dkF2Allintqs(p0 - q0, ps, qsmax, k, m, T) +
        dkF2Allintqs(p0 + q0, ps, qsmax, k, m, T)
    )


V4piImintqs(p0, ps, k, T, Npi) =
    -quadgk(
        x ->
            dkV4piImintqs(p0, ps, Epi(k, msgfun2(k)), k, x, msgfun2(x), T, Npi),
        k,
        Λ,
        rtol = 1e-8,
        atol = 1e-8,
        order = 100,maxevals=8000
    )[1] +deltasumkAll(p0 + Epi(k, msgfun2(k)), ps, k, T, Npi) + deltasumkAll(p0 - Epi(k, msgfun2(k)), ps, k, T, Npi)



propImsimpleintqs(p0, ps, T, Npi) =
    -hcubature(
        k ->
            (
                (V4piImintqs(p0, ps, k[1], T, Npi)) *
                k[1] *
                (
                    -(
                        coth(Epi(k[1], msgfun2(k[1])) / (2 * T)) /
                        Epi(k[1], msgfun2(k[1]))^3
                    ) -
                    csch(Epi(k[1], msgfun2(k[1])) / (2 * T))^2 /
                    (2 * T * Epi(k[1], msgfun2(k[1]))^2)
                )
            ) / (16 * pi^2),
        [kmin],
        [Λ],
        rtol = 1e-8,
        atol = 1e-8,maxevals=8000
    )[1]
