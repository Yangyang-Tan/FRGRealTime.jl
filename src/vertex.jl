
# lam4piflow(k, m, T, Npi, lamda4pi) =
#     lamda4pi^2 * (8 + Npi) * (dkF1TildeAll(k, m, T) + dkF2TildeAll(k, m, T))




lam4piflowIm(p0, ps, q0, qs, costhe, k, m, T, lam4pik, Npi) =
    2 *
    lam4pik^2 *
    pi *
    (
        dkF1All(
            p0 - q0,
            max(0.1, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        dkF1All(
            p0 + q0,
            max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        ((4 + Npi) * (dkF1All(0, 0.1, k, m, T) + dkF2All(0, 0.1, k, m, T))) /
        2 +
        dkF2All(
            p0 - q0,
            max(0.1, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        dkF2All(
            p0 + q0,
            max(0.1, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        )
    )



############################################
#  The original computation of Im parts  #
############################################


#contains two part V(q0)+V(-q0)
dkV4piImsimple(p0, ps, q0, qs, costhe, k, m, T, lamda4pik, Npi) =
    lamda4pik^2 *
    (2 + Npi) *
    pi *
    (
        3 * dkF1All(
            p0 - q0,
            max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * dkF1All(
            p0 - q0,
            max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * dkF1All(
            p0 + q0,
            max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * dkF1All(
            p0 + q0,
            max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * (
            dkF2All(
                p0 - q0,
                max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            dkF2All(
                p0 - q0,
                max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            dkF2All(
                p0 + q0,
                max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            dkF2All(
                p0 + q0,
                max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            )
        )
    )


#Im part the initial condition is zero
V4piImsimple(p0, ps, qs, costhe, k, T, Npi) = quadgk(
        x -> dkV4piImsimple(
            p0,
            ps,
            Epi(k, msgfun2(k)),
            qs,
            costhe,
            x,
            msgfun2(x),
            T,
            lampifun(x),
            Npi,
        ),
        Λ,
        k,
        rtol = 1e-5,
        atol = 1e-5,
        order = 100,
        maxevals = 1024,
    )[1]


dkpropIm(p0, ps, k, Npi) =
    2 * hcubature(
        x -> x[1]^2 * V4piImsimple(p0, ps, x[1], x[2], k, Tc, Npi),
        [0.0, 0.0],
        [k, 1.0],
        atol = 1e-2,
        rtol = 1e-2,
        maxevals = 1000,
    )[1]


parallelpropIm(p0, ps, Npi) = -parallelintegro(
    k ->
        (
            (dkpropIm(p0, ps, k, Npi)) *
            k *
            (
                -(coth(Epi(k, msgfun2(k)) / (2 * Tc)) / Epi(k, msgfun2(k))^3) -
                csch(Epi(k, msgfun2(k)) / (2 * Tc))^2 /
                (2 * Tc * Epi(k, msgfun2(k))^2)
            )
        ) / (16 * pi^2),
    kmin,
    Λ,
    20,
)
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




################################################################################
#  λ not running version, Im part                                              #
################################################################################



dkV4piIm_static(p0, ps, q0, qs, costhe, k, m, T, lamda4pik, Npi) =
    lamda4pik^2 *
    (2 + Npi) *
    pi *
    (
        3 * F1All(
            p0 - q0,
            max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * F1All(
            p0 - q0,
            max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * F1All(
            p0 + q0,
            max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * F1All(
            p0 + q0,
            max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
            k,
            m,
            T,
        ) +
        3 * (
            F2All(
                p0 - q0,
                max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            F2All(
                p0 - q0,
                max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            F2All(
                p0 + q0,
                max(0.01, sqrt(ps^2 - 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            ) +
            F2All(
                p0 + q0,
                max(0.01, sqrt(ps^2 + 2 * costhe * ps * qs + qs^2)),
                k,
                m,
                T,
            )
        )
    )


V4piIm_static(p0, ps, qs, costhe, k, T, Npi) =
    dkV4piIm_static(
        p0,
        ps,
        Epi(k, msgfun2(k)),
        qs,
        costhe,
        k,
        msgfun2(k),
        T,
        lampifun(k),
        Npi,
    ) - dkV4piIm_static(
        p0,
        ps,
        Epi(k, msgfun2(k)),
        qs,
        costhe,
        Λ,
        msgfun2(k),
        T,
        lampifun(k),
        Npi,
    )




dkpropIm_static(p0, ps, k, Npi) =
        2 * hcubature(
            x -> x[1]^2 * V4piIm_static(p0, ps, x[1], x[2], k, Tc, Npi),
            [0.0, 0.0],
            [k, 1.0],
            atol = 1e-2,
            rtol = 1e-2,
            maxevals = 2000,
        )[1]

parallelpropIm_static(p0, ps, Npi) = -parallelintegro(
            k ->
                (
                    (dkpropIm_static(p0, ps, k, Npi)) *
                    k *
                    (
                        -(coth(Epi(k, msgfun2(k)) / (2 * Tc)) / Epi(k, msgfun2(k))^3) -
                        csch(Epi(k, msgfun2(k)) / (2 * Tc))^2 /
                        (2 * Tc * Epi(k, msgfun2(k))^2)
                    )
                ) / (16 * pi^2),
            kmin,
            Λ,
            20,
        )



################################################################################
#  other methods                                                               #
################################################################################

dkV4piIm(p0, ps, q0, qs, k, m, T, Npi) =
    lampifun(k)^2 *
    (2 + Npi) *
    pi *
    (
        6 * dkF1All(p0 - q0, ps, qs, k, m, T) +
        6 * dkF1All(p0 + q0, ps, qs, k, m, T) +
        6 *
        (dkF2All(p0 - q0, ps, qs, k, m, T) + dkF2All(p0 + q0, ps, qs, k, m, T))
    )


V4piIm(p0, ps, k, T, Npi) =
        -hcubature(
            x ->
                x[2]^2 * dkV4piIm(
                    p0,
                    ps,
                    Epi(k, msgfun2(k)),
                    x[2],
                    x[1],
                    msgfun2(x[1]),
                    T,
                    Npi,
                ),
            [k, 0.0],
            [Λ, k],
            rtol = 1e-3,
            atol = 1e-3,
        )[1] +
        deltasumkAll(p0 + Epi(k, msgfun2(k)), ps, k, T, Npi) +
        deltasumkAll(p0 - Epi(k, msgfun2(k)), ps, k, T, Npi)










propImsimple(p0, ps, T, Npi) =
    -hcubature(
        k ->
            (
                (V4piIm(p0, ps, k[1], T, Npi)) *
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
        [0.1],
        [Λ],
        rtol = 1e-3,
        atol = 1e-3,
    )[1]





propImsimple2(p0, ps, T, Npi) = integro(
    k ->
        (
            (V4piIm(p0, ps, k, T, Npi)) *
            k *
            (
                -(coth(Epi(k, msgfun2(k)) / (2 * T)) / Epi(k, msgfun2(k))^3) -
                csch(Epi(k, msgfun2(k)) / (2 * T))^2 /
                (2 * T * Epi(k, msgfun2(k))^2)
            )
        ) / (16 * pi^2),
    Λ,
    1.0,
    10,
)
