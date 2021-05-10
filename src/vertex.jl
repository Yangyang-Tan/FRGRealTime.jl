################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################
dkV4piImintqs(p0, ps, q0, qsmax, k, m, T, Npi) =
    lampifun(k)^2 *
    (2 + Npi) *
    pi *
    (
        6 * dkF1Allintqs(p0 - q0, ps, qsmax, k, m, T)
        +
        6 * dkF1Allintqs(p0 + q0, ps, qsmax, k, m, T)
        +
        6 * (
            dkF2Allintqs(p0 - q0, ps, qsmax, k, m, T) +
            dkF2Allintqs(p0 + q0, ps, qsmax, k, m, T)
        )
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
