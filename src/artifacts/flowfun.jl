function peak(p0, ps, m2, T, δ)
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


function delta2funcosthqs1(p0, ps, qsmax, k, Ek, T)
    -1 / 57600 * (
        k *
        csch(Ek / (2 * T))^2 *
        (
            Ek * (
                (4 * Ek^2 - p0^2) * (
                    (2 * k - ps)^3 * (248 * k^2 + 702 * k * ps + 137 * ps^2) +
                    300 * ps * (-2 * k + ps)^2 * (4 * k + ps) * qsmax -
                    300 * (4 * k^3 - 6 * k^2 * ps + ps^3) * qsmax^2 +
                    200 * (-2 * k^2 + ps^2) * qsmax^3 - 75 * ps * qsmax^4 +
                    12 * qsmax^5
                ) +
                60 * (
                    4 *
                    (
                        8 * k^5 * p0^2 +
                        40 * Ek^2 * k^3 * ps^2 +
                        5 * k^2 * p0^2 * ps^3 +
                        Ek^2 * ps^5
                    ) *
                    log((-2 * k + ps) / qsmax) +
                    (
                        p0^2 * ps^2 * (40 * k^3 + ps^3) +
                        16 * Ek^2 * (8 * k^5 + 5 * k^2 * ps^3)
                    ) * log(-(qsmax / (2 * k - ps)))
                )
            ) +
            T *
            (
                (2 * k - ps + qsmax) * (
                    -1600 *
                    Ek^4 *
                    (
                        10 * k^2 + 17 * k * ps - 11 * ps^2 - 5 * k * qsmax +
                        7 * ps * qsmax - 2 * qsmax^2
                    ) +
                    p0^2 * (
                        -(
                            (-2 * k + ps)^2 *
                            (248 * k^2 + 702 * k * ps + 137 * ps^2)
                        ) +
                        (2 * k - ps) *
                        (248 * k^2 - 498 * k * ps - 163 * ps^2) *
                        qsmax +
                        (352 * k^2 - 102 * k * ps - 137 * ps^2) * qsmax^2 +
                        3 * (8 * k + 21 * ps) * qsmax^3 - 12 * qsmax^4
                    ) +
                    4 *
                    Ek^2 *
                    (
                        (2 * k - ps) * (
                            1488 * k^3 +
                            500 * k * p0^2 +
                            3468 * k^2 * ps +
                            1100 * p0^2 * ps - 1284 * k * ps^2 - 411 * ps^3
                        ) -
                        (
                            1488 * k^3 + 500 * k * p0^2 - 3732 * k^2 * ps -
                            700 * p0^2 * ps +
                            516 * k * ps^2 +
                            489 * ps^3
                        ) * qsmax +
                        (-1056 * k^2 - 200 * p0^2 + 306 * k * ps + 411 * ps^2) *
                        qsmax^2 - 9 * (8 * k + 21 * ps) * qsmax^3 +
                        36 * qsmax^4
                    )
                ) +
                240 *
                (
                    40 * Ek^4 * (4 * k^3 + ps^3) +
                    3 * Ek^2 * ps^2 * (40 * k^3 + 10 * k * p0^2 + ps^3) +
                    k^2 * p0^2 * (8 * k^3 + 5 * ps^3)
                ) *
                log((-2 * k + ps) / qsmax) +
                60 *
                (
                    480 * Ek^4 * k * ps^2 +
                    p0^2 * ps^2 * (40 * k^3 + ps^3) +
                    8 *
                    Ek^2 *
                    (48 * k^5 + 20 * k^3 * p0^2 + 5 * (6 * k^2 + p0^2) * ps^3)
                ) *
                log(-(qsmax / (2 * k - ps)))
            ) *
            sinh(Ek / T)
        )
    ) / (Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T)
end

function delta2funcosthqs2(p0, ps, qsmax, k, Ek, T)
    (
        k * (
            -(
                Ek *
                (4 * Ek^2 - p0^2) *
                (
                    (2 * k - ps)^3 * (248 * k^2 + 702 * k * ps + 137 * ps^2) +
                    300 * ps * (-2 * k + ps)^2 * (4 * k + ps) * qsmax -
                    300 * (4 * k^3 - 6 * k^2 * ps + ps^3) * qsmax^2 +
                    200 * (-2 * k^2 + ps^2) * qsmax^3 - 75 * ps * qsmax^4 +
                    12 * qsmax^5 +
                    60 *
                    (-2 * k + ps)^3 *
                    (4 * k^2 + 6 * k * ps + ps^2) *
                    log((2 * k - ps) / qsmax)
                )
            ) +
            T *
            (
                (2 * k - ps + qsmax) * (
                    1600 *
                    Ek^4 *
                    (
                        10 * k^2 + 17 * k * ps - 11 * ps^2 - 5 * k * qsmax +
                        7 * ps * qsmax - 2 * qsmax^2
                    ) +
                    4 *
                    Ek^2 *
                    (
                        -(
                            (2 * k - ps) * (
                                1488 * k^3 +
                                500 * k * p0^2 +
                                3468 * k^2 * ps +
                                1100 * p0^2 * ps - 1284 * k * ps^2 - 411 * ps^3
                            )
                        ) +
                        (
                            1488 * k^3 + 500 * k * p0^2 - 3732 * k^2 * ps -
                            700 * p0^2 * ps +
                            516 * k * ps^2 +
                            489 * ps^3
                        ) * qsmax +
                        (1056 * k^2 + 200 * p0^2 - 306 * k * ps - 411 * ps^2) *
                        qsmax^2 +
                        9 * (8 * k + 21 * ps) * qsmax^3 - 36 * qsmax^4
                    ) +
                    p0^2 * (
                        (-2 * k + ps)^2 *
                        (248 * k^2 + 702 * k * ps + 137 * ps^2) -
                        (2 * k - ps) *
                        (248 * k^2 - 498 * k * ps - 163 * ps^2) *
                        qsmax +
                        (-352 * k^2 + 102 * k * ps + 137 * ps^2) * qsmax^2 -
                        3 * (8 * k + 21 * ps) * qsmax^3 + 12 * qsmax^4
                    )
                ) +
                60 *
                (-2 * k + ps)^2 *
                (
                    -160 * Ek^4 * (k + ps) +
                    4 *
                    Ek^2 *
                    (
                        24 * k^3 +
                        10 * k * p0^2 +
                        24 * k^2 * ps +
                        10 * p0^2 * ps - 12 * k * ps^2 - 3 * ps^3
                    ) +
                    p0^2 * (-8 * k^3 - 8 * k^2 * ps + 4 * k * ps^2 + ps^3)
                ) *
                log((2 * k - ps) / qsmax)
            ) *
            sinh(Ek / T)
        )
    ) /
    (28800 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T * (-1 + cosh(Ek / T)))
end

function delta2funcosthqs3(p0, ps, qsmax, k, Ek, T)
    (
        k *
        qsmax *
        csch(Ek / (2 * T))^2 *
        (
            -(
                Ek *
                (4 * Ek^2 - p0^2) *
                (
                    75 * ps * (-2 * k + ps)^2 * (4 * k + ps) +
                    50 * (-2 * k^2 + ps^2) * qsmax^2 +
                    3 * qsmax^4
                )
            ) +
            (
                -75 *
                (2 * k - ps) *
                ps *
                (
                    -96 * Ek^4 +
                    p0^2 * (-8 * k^2 + 2 * k * ps + ps^2) +
                    12 * Ek^2 * (2 * p0^2 + (2 * k - ps) * (4 * k + ps))
                ) +
                50 *
                (
                    -16 * Ek^4 +
                    4 * Ek^2 * (6 * k^2 + p0^2 - 3 * ps^2) +
                    p0^2 * (-2 * k^2 + ps^2)
                ) *
                qsmax^2 +
                3 * (-12 * Ek^2 + p0^2) * qsmax^4
            ) *
            T *
            sinh(Ek / T)
        )
    ) / (7200 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T)
end


function delta2funcosthqs4(p0, ps, qsmax, k, Ek, T)
    (
        k *
        csch(Ek / (2 * T))^2 *
        (
            Ek *
            (4 * Ek^2 - p0^2) *
            (
                -137 * ps^4 - 4800 * k^3 * qsmax - 300 * ps^2 * qsmax^2 -
                75 * qsmax^4 +
                200 * k^2 * (11 * ps^2 + 9 * qsmax^2) +
                60 * ps^2 * (-20 * k^2 + ps^2) * log(ps / qsmax)
            ) +
            T *
            (
                4 *
                Ek^2 *
                (
                    1100 * (6 * k^2 + p0^2) * ps^2 - 411 * ps^4 -
                    3600 * k * (4 * k^2 + p0^2) * qsmax +
                    900 * (6 * k^2 + p0^2 - ps^2) * qsmax^2 - 225 * qsmax^4
                ) - 1600 * Ek^4 * (11 * ps^2 + 9 * qsmax * (-4 * k + qsmax)) +
                p0^2 * (
                    137 * ps^4 +
                    4800 * k^3 * qsmax +
                    300 * ps^2 * qsmax^2 +
                    75 * qsmax^4 - 200 * k^2 * (11 * ps^2 + 9 * qsmax^2)
                ) +
                60 *
                ps^2 *
                (
                    20 * (8 * Ek^4 + k^2 * p0^2 - 2 * Ek^2 * (6 * k^2 + p0^2)) +
                    (12 * Ek^2 - p0^2) * ps^2
                ) *
                log(ps / qsmax)
            ) *
            sinh(Ek / T)
        )
    ) / (28800 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * T)
end
