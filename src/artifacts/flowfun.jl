function peak(p0, ps, m2, T,δ)
    ((sqrt(-4 * m2 + p0^2) - ps) * coth(p0 / (4 * T))) / (8 * p0 * δ * pi^2)
end

function ppfun(p0, ps, k, Ek, T)
    (
        k * (
            (
                2 *
                (
                    p0^2 * (sqrt(k^2 - 2 * Ek * p0 + p0^2) - 2 * ps) +
                    2 * Ek * p0 * ps +
                    ps *
                    (2 * Ek^2 - 2 * k^2 + sqrt(k^2 - 2 * Ek * p0 + p0^2) * ps)
                ) *
                coth(Ek / (2 * T))
            ) / sqrt(k^2 + p0 * (-2 * Ek + p0)) -
            (
                2 *
                (
                    p0^2 * (sqrt(k^2 - 2 * Ek * p0 + p0^2) - 2 * ps) +
                    2 * Ek * p0 * ps +
                    ps *
                    (2 * Ek^2 - 2 * k^2 + sqrt(k^2 - 2 * Ek * p0 + p0^2) * ps)
                ) *
                coth((Ek - p0) / (2 * T))
            ) / sqrt(k^2 + p0 * (-2 * Ek + p0)) -
            (
                16 *
                Ek *
                exp((2 * Ek + p0) / T) *
                (
                    2 * Ek * T * cosh(p0 / (2 * T)) -
                    2 * Ek * T * cosh((-2 * Ek + p0) / (2 * T)) +
                    (
                        -2 * Ek * p0 + p0^2 -
                        2 * sqrt(k^2 - 2 * Ek * p0 + p0^2) * ps + ps^2
                    ) * sinh((2 * Ek - p0) / (2 * T))
                ) *
                sinh(p0 / (2 * T))
            ) / (T * (-1 + exp(Ek / T))^2 * (exp(Ek / T) - exp(p0 / T))^2)
        )
    ) / (64 * Ek^3 * pi^2 * ps)
end
function pmfun(p0, ps, k, Ek, T)
    (
        4 *
        Ek^2 *
        k *
        ((1 - exp(Ek / T))^(-1) + (-1 + exp((Ek + p0) / T))^(-1)) +
        (
            k *
            csch(Ek / (2 * T))^2 *
            csch((Ek + p0) / (2 * T))^2 *
            sinh(p0 / (2 * T)) *
            (
                (
                    -(p0^2 * (sqrt(k^2 + 2 * Ek * p0 + p0^2) - 2 * ps)) +
                    2 * Ek * p0 * ps -
                    ps *
                    (2 * Ek^2 - 2 * k^2 + sqrt(k^2 + 2 * Ek * p0 + p0^2) * ps)
                ) *
                T *
                cosh(p0 / (2 * T)) +
                (
                    p0^2 * (sqrt(k^2 + 2 * Ek * p0 + p0^2) - 2 * ps) -
                    2 * Ek * p0 * ps +
                    ps *
                    (2 * Ek^2 - 2 * k^2 + sqrt(k^2 + 2 * Ek * p0 + p0^2) * ps)
                ) *
                T *
                cosh((2 * Ek + p0) / (2 * T)) +
                Ek *
                (
                    2 * Ek * p0 * sqrt(k^2 + p0 * (2 * Ek + p0)) +
                    p0^2 * sqrt(k^2 + p0 * (2 * Ek + p0)) - 4 * Ek * p0 * ps -
                    2 * (k^2 + p0^2) * ps +
                    sqrt(k^2 + p0 * (2 * Ek + p0)) * ps^2
                ) *
                sinh((2 * Ek + p0) / (2 * T))
            )
        ) / (2 * sqrt(k^2 + p0 * (2 * Ek + p0)) * T)
    ) / (32 * Ek^3 * pi^2 * ps)
end

function ppfunps(p0, psmax, k, m, T)
    (
        k *
        csch(sqrt(k^2 + m) / (2 * T))^2 *
        csch((sqrt(k^2 + m) - p0) / (2 * T))^2 *
        (
            -(
                (
                    -8 * k^3 * (m + (sqrt(k^2 + m) - p0) * p0) +
                    p0^4 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) -
                    6 *
                    p0^2 *
                    sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) *
                    psmax^2 - 8 * sqrt(k^2 + m) * p0 * psmax^3 +
                    8 * p0^2 * psmax^3 -
                    3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^4 +
                    4 *
                    k^2 *
                    sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) *
                    (2 * m + 2 * sqrt(k^2 + m) * p0 - 3 * p0^2 + 3 * psmax^2) +
                    4 *
                    m *
                    (
                        -2 * p0^2 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) +
                        2 *
                        p0 *
                        sqrt((k^2 + m) * (k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)) +
                        3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2 -
                        2 * psmax^3
                    )
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                -8 * k^3 * (m + (sqrt(k^2 + m) - p0) * p0) +
                p0^4 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) -
                6 * p0^2 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2 -
                8 * sqrt(k^2 + m) * p0 * psmax^3 + 8 * p0^2 * psmax^3 -
                3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^4 +
                4 *
                k^2 *
                sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) *
                (2 * m + 2 * sqrt(k^2 + m) * p0 - 3 * p0^2 + 3 * psmax^2) +
                4 *
                m *
                (
                    -2 * p0^2 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) +
                    2 *
                    p0 *
                    sqrt((k^2 + m) * (k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)) +
                    3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2 -
                    2 * psmax^3
                )
            ) *
            T *
            cosh((-2 * sqrt(k^2 + m) + p0) / (2 * T)) +
            (
                8 * k^5 * (sqrt(k^2 + m) - 2 * p0) +
                8 * k^3 * p0 * (-2 * m + sqrt(k^2 + m) * p0) +
                8 *
                k^4 *
                (-sqrt(k^2 + m) + p0) *
                sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) +
                4 *
                m *
                p0 *
                (
                    p0 *
                    sqrt((k^2 + m) * (k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)) -
                    p0^2 * sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) +
                    (
                        3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) -
                        4 * psmax
                    ) * psmax^2
                ) +
                4 *
                k^2 *
                (
                    -(p0^3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)) +
                    2 * m * p0 * sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) +
                    p0 *
                    (
                        3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) -
                        4 * psmax
                    ) *
                    psmax^2 +
                    2 * sqrt(k^2 + m) * psmax^3
                ) +
                sqrt(k^2 + m) * (
                    p0^4 * sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) -
                    3 * sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) * psmax^4 +
                    2 *
                    p0^2 *
                    psmax^2 *
                    (
                        -3 * sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2) +
                        4 * psmax
                    )
                )
            ) * sinh((2 * sqrt(k^2 + m) - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (
        768 *
        (k^2 + m)^(3 / 2) *
        sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0)) *
        pi^2 *
        T
    )
end

function pmfunps(p0, psmax, k, m, T)
    -(
        k *
        csch(sqrt(k^2 + m) / (2 * T)) *
        csch((sqrt(k^2 + m) + p0) / (2 * T)) *
        sinh(p0 / (2 * T)) *
        (
            2 *
            (
                p0^4 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) +
                8 * k^3 * (-m + p0 * (sqrt(k^2 + m) + p0)) -
                6 * p0^2 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2 +
                8 * sqrt(k^2 + m) * p0 * psmax^3 +
                8 * p0^2 * psmax^3 -
                3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^4 +
                4 *
                k^2 *
                (
                    2 * m * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) -
                    3 * p0^2 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) -
                    2 *
                    p0 *
                    sqrt((k^2 + m) * (k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2)) +
                    3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2
                ) -
                4 *
                m *
                (
                    2 * p0^2 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) +
                    2 *
                    p0 *
                    sqrt((k^2 + m) * (k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2)) -
                    3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) * psmax^2 +
                    2 * psmax^3
                )
            ) *
            T +
            (
                8 * k^5 * (sqrt(k^2 + m) + 2 * p0) +
                8 * k^3 * p0 * (2 * m + sqrt(k^2 + m) * p0) +
                p0^4 * sqrt((k^2 + m) * (k^2 + p0 * (2 * sqrt(k^2 + m) + p0))) -
                8 *
                k^4 *
                (
                    p0 * sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0)) +
                    sqrt((k^2 + m) * (k^2 + p0 * (2 * sqrt(k^2 + m) + p0)))
                ) -
                6 *
                p0^2 *
                sqrt((k^2 + m) * (k^2 + p0 * (2 * sqrt(k^2 + m) + p0))) *
                psmax^2 + 8 * sqrt(k^2 + m) * p0^2 * psmax^3 -
                3 *
                sqrt((k^2 + m) * (k^2 + p0 * (2 * sqrt(k^2 + m) + p0))) *
                psmax^4 +
                4 *
                m *
                p0 *
                (
                    p0^2 * sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0)) +
                    p0 *
                    sqrt((k^2 + m) * (k^2 + p0 * (2 * sqrt(k^2 + m) + p0))) +
                    psmax^2 * (
                        -3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) +
                        4 * psmax
                    )
                ) +
                4 *
                k^2 *
                (
                    p0^3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) -
                    2 * m * p0 * sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0)) +
                    2 * sqrt(k^2 + m) * psmax^3 +
                    p0 *
                    psmax^2 *
                    (
                        -3 * sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2) +
                        4 * psmax
                    )
                )
            ) *
            csch(sqrt(k^2 + m) / (2 * T)) *
            csch((sqrt(k^2 + m) + p0) / (2 * T)) *
            sinh((2 * sqrt(k^2 + m) + p0) / (2 * T))
        )
    ) / (
        768 *
        (k^2 + m)^(3 / 2) *
        sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0)) *
        pi^2 *
        T
    )
end

function ppfuncosthe(p0, ps, qs, k, Ek, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs)) -
            min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)
        ) *
        (
            T *
            cosh(p0 / (2 * T)) *
            (
                3 * (2 * Ek^2 - p0^2) * sqrt(k^2 - 2 * Ek * p0 + p0^2) -
                sqrt(k^2 + p0 * (-2 * Ek + p0)) *
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs))^2 -
                3 *
                (Ek^2 - k^2 + Ek * p0 - p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs) -
                sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)^2 -
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs)) * (
                    3 * Ek^2 - 3 * k^2 + 3 * Ek * p0 - 3 * p0^2 +
                    sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                    min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)
                )
            ) -
            T *
            cosh((-2 * Ek + p0) / (2 * T)) *
            (
                3 * (2 * Ek^2 - p0^2) * sqrt(k^2 - 2 * Ek * p0 + p0^2) -
                sqrt(k^2 + p0 * (-2 * Ek + p0)) *
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs))^2 -
                3 *
                (Ek^2 - k^2 + Ek * p0 - p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs) -
                sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)^2 -
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs)) * (
                    3 * Ek^2 - 3 * k^2 + 3 * Ek * p0 - 3 * p0^2 +
                    sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                    min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)
                )
            ) +
            Ek *
            (
                3 * p0 * (-2 * Ek + p0) * sqrt(k^2 - 2 * Ek * p0 + p0^2) +
                sqrt(k^2 + p0 * (-2 * Ek + p0)) *
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs))^2 -
                3 *
                (k^2 - 2 * Ek * p0 + p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs) +
                sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)^2 +
                max(-k + sqrt(k^2 + p0 * (-2 * Ek + p0)), abs(ps - qs)) * (
                    -3 * k^2 + 6 * Ek * p0 - 3 * p0^2 +
                    sqrt(k^2 - 2 * Ek * p0 + p0^2) *
                    min(k + sqrt(k^2 + p0 * (-2 * Ek + p0)), ps + qs)
                )
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (192 * Ek^3 * sqrt(k^2 + p0 * (-2 * Ek + p0)) * pi^2 * ps * qs * T)
end
function ppfuncostheps1(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (Ek2 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                -(
                    (Ek2 + k - ps) * (
                        5 * Ek * p0 * (Ek2 - 3 * k - ps) +
                        Ek2 * (Ek2 - 4 * k - ps) * (Ek2 + k - ps) +
                        5 * Ek^2 * (-Ek2 + 3 * k + ps)
                    )
                ) +
                2 *
                (
                    5 * Ek * p0 * (Ek2 - 3 * k - ps) +
                    Ek2 * (Ek2 - 4 * k - ps) * (Ek2 + k - ps) +
                    5 * Ek^2 * (-Ek2 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (-Ek + p0) + Ek2 * (7 * Ek2 - 8 * k - 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                -5 *
                Ek^2 *
                (
                    Ek2^2 - 3 * k^2 - 2 * Ek2 * (k + ps + qsmax) +
                    2 * k * (ps + 3 * qsmax) +
                    (ps - qsmax) * (ps + 3 * qsmax)
                ) +
                5 *
                Ek *
                p0 *
                (
                    Ek2^2 - 3 * k^2 - 2 * Ek2 * (k + ps + qsmax) +
                    2 * k * (ps + 3 * qsmax) +
                    (ps - qsmax) * (ps + 3 * qsmax)
                ) +
                Ek2 *
                (Ek2 + k - ps + qsmax) *
                (
                    Ek2^2 - 4 * k^2 - Ek2 * (3 * k + 2 * ps + 3 * qsmax) +
                    3 * k * (ps + 4 * qsmax) +
                    (ps - qsmax) * (ps + 4 * qsmax)
                )
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek2 *
            (Ek2 + k - ps + qsmax) *
            (
                Ek2^2 - 4 * k^2 - Ek2 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end
function ppfuncostheps2(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (Ek2 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                -(
                    (Ek2 + k - ps) * (
                        5 * Ek * p0 * (Ek2 - 3 * k - ps) +
                        Ek2 * (Ek2 - 4 * k - ps) * (Ek2 + k - ps) +
                        5 * Ek^2 * (-Ek2 + 3 * k + ps)
                    )
                ) +
                2 *
                (
                    5 * Ek * p0 * (Ek2 - 3 * k - ps) +
                    Ek2 * (Ek2 - 4 * k - ps) * (Ek2 + k - ps) +
                    5 * Ek^2 * (-Ek2 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (-Ek + p0) + Ek2 * (7 * Ek2 - 8 * k - 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                -5 *
                Ek^2 *
                (
                    Ek2^2 - 3 * k^2 - 2 * Ek2 * (k + ps + qsmax) +
                    2 * k * (ps + 3 * qsmax) +
                    (ps - qsmax) * (ps + 3 * qsmax)
                ) +
                5 *
                Ek *
                p0 *
                (
                    Ek2^2 - 3 * k^2 - 2 * Ek2 * (k + ps + qsmax) +
                    2 * k * (ps + 3 * qsmax) +
                    (ps - qsmax) * (ps + 3 * qsmax)
                ) +
                Ek2 *
                (Ek2 + k - ps + qsmax) *
                (
                    Ek2^2 - 4 * k^2 - Ek2 * (3 * k + 2 * ps + 3 * qsmax) +
                    3 * k * (ps + 4 * qsmax) +
                    (ps - qsmax) * (ps + 4 * qsmax)
                )
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek2 *
            (Ek2 + k - ps + qsmax) *
            (
                Ek2^2 - 4 * k^2 - Ek2 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end

function ppfuncostheps3(p0, ps, qsmax, k, Ek, Ek2, T)
    -(
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                10 * Ek^2 * (Ek2 - ps) + 10 * Ek * p0 * (-Ek2 + ps) -
                Ek2 * (-5 * k^2 + 5 * (Ek2 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((-2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek2 *
            (-5 * k^2 + 5 * (Ek2 - ps)^2 + qsmax^2) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (480 * Ek^3 * Ek2 * pi^2 * ps * T)
end


function ppfuncostheps4(p0, ps, qsmax, k, Ek, Ek2, T)
    -(
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                10 * Ek^2 * (Ek2 - ps) + 10 * Ek * p0 * (-Ek2 + ps) -
                Ek2 * (-5 * k^2 + 5 * (Ek2 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((-2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek2 *
            (-5 * k^2 + 5 * (Ek2 - ps)^2 + qsmax^2) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (480 * Ek^3 * Ek2 * pi^2 * ps * T)
end


function ppfuncostheps5(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (-Ek2 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek2 *
            (Ek2 - k - ps - qsmax) *
            (
                Ek2^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek2 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end

function ppfuncostheps6(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (-Ek2 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek2 *
            (Ek2 - k - ps - qsmax) *
            (
                Ek2^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek2 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end

function ppfuncostheps7(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek2 *
            (
                5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek2 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (1920 * Ek^3 * Ek2 * pi^2 * T)
end

function ppfuncostheps8(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek2 *
            (
                5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek2 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (1920 * Ek^3 * Ek2 * pi^2 * T)
end

function ppfuncostheps9(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (-Ek2 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek2 *
            (Ek2 - k - ps - qsmax) *
            (
                Ek2^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek2 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end

function ppfuncostheps10(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                20 *
                Ek^2 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                20 *
                Ek *
                p0 *
                (
                    Ek2^3 +
                    Ek2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek2 * (
                    5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek2 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek2 *
            (
                5 * Ek2^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek2^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek2 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (1920 * Ek^3 * Ek2 * pi^2 * T)
end

function ppfuncostheps11(p0, ps, qsmax, k, Ek, Ek2, T)
    (
        k *
        (-Ek2 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek - p0) / (2 * T))^2 *
        (
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) -
            (
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) * (Ek2 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek2 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek2 + 3 * k - ps) +
                    Ek2 * (Ek2 - k - ps) * (Ek2 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek - p0) + Ek2 * (-7 * Ek2 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek2 * qsmax^3
            ) *
            T *
            cosh((-2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek2 *
            (Ek2 - k - ps - qsmax) *
            (
                Ek2^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek2 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek - p0) / (2 * T))
        ) *
        sinh(p0 / (2 * T))
    ) / (3840 * Ek^3 * Ek2 * pi^2 * ps * T)
end
function pmfuncosthe1(p0, ps, qs, k, Ek, T)
    (
        k *
        csch(Ek / (2 * T)) *
        csch((Ek + p0) / (2 * T)) *
        (
            max(2 * k, abs(ps - qs)) -
            min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs)
        ) *
        sinh(p0 / (2 * T)) *
        (
            12 * Ek^2 - 6 * p0^2 +
            (
                csch(Ek / (2 * T)) *
                csch((Ek + p0) / (2 * T)) *
                (
                    -3 *
                    Ek *
                    p0 *
                    (2 * Ek + p0) *
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    sinh((2 * Ek + p0) / (2 * T)) -
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    max(2 * k, abs(ps - qs))^2 *
                    (
                        -(T * cosh(p0 / (2 * T))) +
                        T * cosh((2 * Ek + p0) / (2 * T)) +
                        Ek * sinh((2 * Ek + p0) / (2 * T))
                    ) -
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs)^2 *
                    (
                        -(T * cosh(p0 / (2 * T))) +
                        T * cosh((2 * Ek + p0) / (2 * T)) +
                        Ek * sinh((2 * Ek + p0) / (2 * T))
                    ) +
                    3 *
                    min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs) *
                    (
                        (Ek^2 - k^2 - Ek * p0 - p0^2) *
                        T *
                        (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
                        Ek *
                        (k^2 + 2 * Ek * p0 + p0^2) *
                        sinh((2 * Ek + p0) / (2 * T))
                    ) +
                    max(2 * k, abs(ps - qs)) * (
                        3 *
                        (Ek^2 - k^2 - Ek * p0 - p0^2) *
                        T *
                        (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
                        3 *
                        Ek *
                        (k^2 + 2 * Ek * p0 + p0^2) *
                        sinh((2 * Ek + p0) / (2 * T)) -
                        sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                        min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs) *
                        (
                            -(T * cosh(p0 / (2 * T))) +
                            T * cosh((2 * Ek + p0) / (2 * T)) +
                            Ek * sinh((2 * Ek + p0) / (2 * T))
                        )
                    )
                )
            ) / (sqrt(k^2 + p0 * (2 * Ek + p0)) * T)
        )
    ) / (192 * Ek^3 * pi^2 * ps * qs)
end


function pmfuncosthe2(p0, ps, qs, k, Ek, T)
    (
        k *
        csch(Ek / (2 * T)) *
        csch((Ek + p0) / (2 * T)) *
        (
            max(-k + sqrt(k^2 + 2 * Ek * p0 + p0^2), abs(ps - qs)) -
            min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs)
        ) *
        sinh(p0 / (2 * T)) *
        (
            12 * Ek^2 - 6 * p0^2 +
            (
                csch(Ek / (2 * T)) *
                csch((Ek + p0) / (2 * T)) *
                (
                    -3 *
                    Ek *
                    p0 *
                    (2 * Ek + p0) *
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    sinh((2 * Ek + p0) / (2 * T)) -
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    max(-k + sqrt(k^2 + 2 * Ek * p0 + p0^2), abs(ps - qs))^2 *
                    (
                        -(T * cosh(p0 / (2 * T))) +
                        T * cosh((2 * Ek + p0) / (2 * T)) +
                        Ek * sinh((2 * Ek + p0) / (2 * T))
                    ) -
                    sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                    min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs)^2 *
                    (
                        -(T * cosh(p0 / (2 * T))) +
                        T * cosh((2 * Ek + p0) / (2 * T)) +
                        Ek * sinh((2 * Ek + p0) / (2 * T))
                    ) +
                    3 *
                    min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs) *
                    (
                        (Ek^2 - k^2 - Ek * p0 - p0^2) *
                        T *
                        (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
                        Ek *
                        (k^2 + 2 * Ek * p0 + p0^2) *
                        sinh((2 * Ek + p0) / (2 * T))
                    ) +
                    max(-k + sqrt(k^2 + 2 * Ek * p0 + p0^2), abs(ps - qs)) * (
                        3 *
                        (Ek^2 - k^2 - Ek * p0 - p0^2) *
                        T *
                        (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
                        3 *
                        Ek *
                        (k^2 + 2 * Ek * p0 + p0^2) *
                        sinh((2 * Ek + p0) / (2 * T)) -
                        sqrt(k^2 + 2 * Ek * p0 + p0^2) *
                        min(k + sqrt(k^2 + 2 * Ek * p0 + p0^2), ps + qs) *
                        (
                            -(T * cosh(p0 / (2 * T))) +
                            T * cosh((2 * Ek + p0) / (2 * T)) +
                            Ek * sinh((2 * Ek + p0) / (2 * T))
                        )
                    )
                )
            ) / (sqrt(k^2 + p0 * (2 * Ek + p0)) * T)
        )
    ) / (192 * Ek^3 * pi^2 * ps * qs)
end

function pmfuncostheps1(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (Ek1 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                (Ek1 + k - ps) * (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) -
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 - 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                (-Ek1 - k + ps) * (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) +
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax -
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 + 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (Ek1 + k - ps + qsmax) *
            (
                Ek1^2 - 4 * k^2 - Ek1 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end

function pmfuncostheps2(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-2 * k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (2 * k - ps) * (
                        5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                        5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                        Ek1 * (
                            10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                            5 * Ek1 * (6 * k + ps)
                        )
                    ) +
                    2 *
                    (
                        5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                        5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                        Ek1 * (
                            10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                            5 * Ek1 * (6 * k + ps)
                        )
                    ) *
                    qsmax +
                    (
                        15 * Ek * (Ek + p0) +
                        Ek1 * (-15 * Ek1 + 16 * k + 7 * ps)
                    ) * qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (2 * k - ps) * (
                    5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                    5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                    Ek1 * (
                        10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                        5 * Ek1 * (6 * k + ps)
                    )
                ) +
                2 *
                (
                    5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                    5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                    Ek1 * (
                        10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                        5 * Ek1 * (6 * k + ps)
                    )
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-15 * Ek1 + 16 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (
                (2 * k - ps) * (
                    10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                    5 * Ek1 * (6 * k + ps)
                ) +
                2 *
                (
                    10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                    5 * Ek1 * (6 * k + ps)
                ) *
                qsmax +
                (-15 * Ek1 + 16 * k + 7 * ps) * qsmax^2 +
                4 * qsmax^3
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end

function pmfuncostheps3(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                10 * Ek^2 * (Ek1 - ps) + 10 * Ek * p0 * (Ek1 - ps) -
                Ek1 * (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek1 *
            (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (480 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps4(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (Ek1 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                (Ek1 + k - ps) * (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) -
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 - 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                -(
                    (Ek1 + k - ps) * (
                        Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                        5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                        5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                    )
                ) +
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax -
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 + 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (Ek1 + k - ps + qsmax) *
            (
                Ek1^2 - 4 * k^2 - Ek1 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps5(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-2 * k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (2 * k - ps) * (
                        5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                        5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                        Ek1 * (
                            10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                            5 * Ek1 * (6 * k + ps)
                        )
                    ) +
                    2 *
                    (
                        5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                        5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                        Ek1 * (
                            10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                            5 * Ek1 * (6 * k + ps)
                        )
                    ) *
                    qsmax +
                    (
                        15 * Ek * (Ek + p0) +
                        Ek1 * (-15 * Ek1 + 16 * k + 7 * ps)
                    ) * qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (2 * k - ps) * (
                    5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                    5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                    Ek1 * (
                        10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                        5 * Ek1 * (6 * k + ps)
                    )
                ) +
                2 *
                (
                    5 * Ek^2 * (-4 * Ek1 + 6 * k + ps) +
                    5 * Ek * p0 * (-4 * Ek1 + 6 * k + ps) +
                    Ek1 * (
                        10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                        5 * Ek1 * (6 * k + ps)
                    )
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-15 * Ek1 + 16 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (
                (2 * k - ps) * (
                    10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                    5 * Ek1 * (6 * k + ps)
                ) +
                2 *
                (
                    10 * Ek1^2 + 14 * k^2 + 6 * k * ps + ps^2 -
                    5 * Ek1 * (6 * k + ps)
                ) *
                qsmax +
                (-15 * Ek1 + 16 * k + 7 * ps) * qsmax^2 +
                4 * qsmax^3
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps6(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (Ek1 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                (Ek1 + k - ps) * (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) -
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 - 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                -(
                    (Ek1 + k - ps) * (
                        Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                        5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                        5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                    )
                ) +
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax -
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 + 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (Ek1 + k - ps + qsmax) *
            (
                Ek1^2 - 4 * k^2 - Ek1 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end

function pmfuncostheps7(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                10 * Ek^2 * (Ek1 - ps) + 10 * Ek * p0 * (Ek1 - ps) -
                Ek1 * (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek1 *
            (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (480 * Ek^3 * Ek1 * pi^2 * ps * T)
end

function pmfuncostheps8(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (Ek1 + k - ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                (Ek1 + k - ps) * (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) -
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 - 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                -(
                    (Ek1 + k - ps) * (
                        Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                        5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                        5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                    )
                ) +
                2 *
                (
                    Ek1 * (Ek1 - 4 * k - ps) * (Ek1 + k - ps) +
                    5 * Ek^2 * (-Ek1 + 3 * k + ps) +
                    5 * Ek * p0 * (-Ek1 + 3 * k + ps)
                ) *
                qsmax -
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 + 8 * k + 7 * ps)) *
                qsmax^2 + 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (Ek1 + k - ps + qsmax) *
            (
                Ek1^2 - 4 * k^2 - Ek1 * (3 * k + 2 * ps + 3 * qsmax) +
                3 * k * (ps + 4 * qsmax) +
                (ps - qsmax) * (ps + 4 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps9(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                10 * Ek^2 * (Ek1 - ps) + 10 * Ek * p0 * (Ek1 - ps) -
                Ek1 * (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek1 *
            (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (480 * Ek^3 * Ek1 * pi^2 * ps * T)
end

function pmfuncostheps10(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-Ek1 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) * (Ek1 - k - ps) +
                    2 *
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) *
                    qsmax +
                    (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                    qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (Ek1 - k - ps - qsmax) *
            (
                Ek1^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek1 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps11(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    20 *
                    Ek^2 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) +
                    20 *
                    Ek *
                    p0 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) -
                    Ek1 * (
                        5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                        15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                        10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        40 * Ek1 * (k^3 + qsmax^3)
                    )
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                20 *
                Ek^2 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) +
                20 *
                Ek *
                p0 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek1 * (
                    5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek1 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (
                5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek1 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (1920 * Ek^3 * Ek1 * pi^2 * T)
end


function pmfuncostheps12(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-Ek1 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) * (Ek1 - k - ps) +
                    2 *
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) *
                    qsmax +
                    (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                    qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (Ek1 - k - ps - qsmax) *
            (
                Ek1^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek1 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps13(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                10 * Ek^2 * (Ek1 - ps) + 10 * Ek * p0 * (Ek1 - ps) -
                Ek1 * (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2)
            ) *
            T *
            (cosh(p0 / (2 * T)) - cosh((2 * Ek + p0) / (2 * T))) +
            Ek *
            Ek1 *
            (-5 * k^2 + 5 * (Ek1 - ps)^2 + qsmax^2) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (480 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps14(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-Ek1 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) * (Ek1 - k - ps) +
                    2 *
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) *
                    qsmax +
                    (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                    qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (Ek1 - k - ps - qsmax) *
            (
                Ek1^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek1 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps15(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    20 *
                    Ek^2 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) +
                    20 *
                    Ek *
                    p0 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) -
                    Ek1 * (
                        5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                        15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                        10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        40 * Ek1 * (k^3 + qsmax^3)
                    )
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                20 *
                Ek^2 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) +
                20 *
                Ek *
                p0 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek1 * (
                    5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek1 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (
                5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek1 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (1920 * Ek^3 * Ek1 * pi^2 * T)
end


function pmfuncostheps16(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-Ek1 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) * (Ek1 - k - ps) +
                    2 *
                    (
                        -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                        5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                        Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                    ) *
                    qsmax +
                    (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                    qsmax^2 +
                    4 * Ek1 * qsmax^3
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (Ek1 - k - ps - qsmax) *
            (
                Ek1^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek1 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end


function pmfuncostheps17(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            -(
                (
                    20 *
                    Ek^2 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) +
                    20 *
                    Ek *
                    p0 *
                    (
                        Ek1^3 +
                        Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        2 * (k^3 + qsmax^3)
                    ) -
                    Ek1 * (
                        5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                        15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                        10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                        40 * Ek1 * (k^3 + qsmax^3)
                    )
                ) *
                T *
                cosh(p0 / (2 * T))
            ) +
            (
                20 *
                Ek^2 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) +
                20 *
                Ek *
                p0 *
                (
                    Ek1^3 +
                    Ek1 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    2 * (k^3 + qsmax^3)
                ) -
                Ek1 * (
                    5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                    15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                    10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                    40 * Ek1 * (k^3 + qsmax^3)
                )
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) -
            Ek *
            Ek1 *
            (
                5 * Ek1^4 - 15 * k^4 + ps^4 - 10 * ps^2 * qsmax^2 -
                15 * qsmax^4 - 10 * k^2 * (ps^2 - 3 * qsmax^2) +
                10 * Ek1^2 * (-3 * k^2 + ps^2 - 3 * qsmax^2) +
                40 * Ek1 * (k^3 + qsmax^3)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (1920 * Ek^3 * Ek1 * pi^2 * T)
end

function pmfuncostheps18(p0, ps, qsmax, k, Ek, Ek1, T)
    (
        k *
        (-Ek1 + k + ps + qsmax)^2 *
        csch(Ek / (2 * T))^2 *
        csch((Ek + p0) / (2 * T))^2 *
        sinh(p0 / (2 * T)) *
        (
            (
                (
                    5 * Ek^2 * (Ek1 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) -
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    5 * Ek^2 * (Ek1 + 3 * k - ps) +
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) -
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (-15 * Ek * (Ek + p0) + Ek1 * (7 * Ek1 + 8 * k - 7 * ps)) *
                qsmax^2 - 4 * Ek1 * qsmax^3
            ) *
            T *
            cosh(p0 / (2 * T)) +
            (
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) * (Ek1 - k - ps) +
                2 *
                (
                    -5 * Ek^2 * (Ek1 + 3 * k - ps) -
                    5 * Ek * p0 * (Ek1 + 3 * k - ps) +
                    Ek1 * (Ek1 - k - ps) * (Ek1 + 4 * k - ps)
                ) *
                qsmax +
                (15 * Ek * (Ek + p0) + Ek1 * (-7 * Ek1 - 8 * k + 7 * ps)) *
                qsmax^2 +
                4 * Ek1 * qsmax^3
            ) *
            T *
            cosh((2 * Ek + p0) / (2 * T)) +
            Ek *
            Ek1 *
            (Ek1 - k - ps - qsmax) *
            (
                Ek1^2 - 4 * k^2 - 3 * k * (ps - 4 * qsmax) +
                (ps - 4 * qsmax) * (ps + qsmax) +
                Ek1 * (3 * k - 2 * ps + 3 * qsmax)
            ) *
            sinh((2 * Ek + p0) / (2 * T))
        )
    ) / (3840 * Ek^3 * Ek1 * pi^2 * ps * T)
end
