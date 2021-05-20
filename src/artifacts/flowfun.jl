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


function ppfunps_delta(p0, qsmax, k, m, T)
    (
        qsmax^3 *
        csch(Epi(k, m) / (2 * T))^2 *
        (
            k *
            (32 * k^3 - 18 * k^2 * qsmax + qsmax^3) *
            Epi(k, m) *
            (p0^2 - 4 * Epi(k, m)^2) +
            k *
            T *
            (
                p0^2 * (32 * k^3 - 18 * k^2 * qsmax + qsmax^3) -
                12 *
                (
                    8 * k * (4 * k^2 + p0^2) - 3 * (6 * k^2 + p0^2) * qsmax +
                    qsmax^3
                ) *
                Epi(k, m)^2 + 48 * (8 * k - 3 * qsmax) * Epi(k, m)^4
            ) *
            sinh(Epi(k, m) / T)
        )
    ) / (1152 * pi^2 * T * Epi(k, m)^3 * (p0^2 - 4 * Epi(k, m)^2)^2)
end

function deltasumps(p0, psmax, k, m2, T, δk = 0.02)
    if p0 + δk > 2 * Epi(k, m2) > p0 - δk
        return (
            k *
            (sqrt(k^2 + m2) + p0) *
            (8 * k - 3 * psmax) *
            psmax^3 *
            coth(sqrt(k^2 + m2) / (2 * T))
        ) / (24 * p0^2 * (2 * k^2 + 2 * m2 + sqrt(k^2 + m2) * p0) * pi^2)
    else
        return (
            k * psmax^3 * (-8 * k + 3 * psmax) * coth(sqrt(k^2 + m2) / (2 * T))
        ) / (24 * sqrt(k^2 + m2) * (4 * (k^2 + m2) - p0^2) * pi^2)
    end
end

function pmfun_zerofix(p0, p, k, m, T, a)
    (
        k *
        (2 * k - p) *
        csch(Epi(k, m) / (2 * T))^2 *
        (
            2 * (2 * k - p) * (4 * k + p) * T +
            (2 * k - p) * (4 * k + p) * coth(Epi(k, m) / (2 * T)) * Epi(k, m) -
            24 * T * Epi(k, m)^2
        ) *
        (-1 + tanh(10^a * p0))
    ) / (768 * pi^2 * T^2 * Epi(k, m)^4)
end

function pmfunps_zerofix(p0, psmax, k, m, T, a)
    (
        k *
        psmax^3 *
        (1 + coth(Epi(k, m) / (2 * T))) *
        (
            2 * (32 * k^3 - 18 * k^2 * psmax + psmax^3) * T +
            (32 * k^3 - 18 * k^2 * psmax + psmax^3) *
            coth(Epi(k, m) / (2 * T)) *
            Epi(k, m) +
            12 * (-8 * k + 3 * psmax) * T * Epi(k, m)^2
        ) *
        (-1 + tanh(10^a * p0))
    ) / (2304 * pi^2 * T^2 * Epi(k, m)^4 * (-1 + exp(Epi(k, m) / T)))
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
    (
        k *
        (2 * k - ps + qsmax)^3 *
        csch(Ek / (2 * T))^2 *
        (
            Ek *
            (4 * Ek^2 - p0^2) *
            (2 * k - ps + qsmax) *
            (
                36 * k^3 + 2 * k^2 * (ps - 36 * qsmax) -
                8 * k * (ps - qsmax) * (ps + 6 * qsmax) -
                (ps - qsmax)^2 * (ps + 6 * qsmax)
            ) +
            (
                -(
                    p0^2 *
                    (2 * k - ps + qsmax) *
                    (
                        36 * k^3 + 2 * k^2 * (ps - 36 * qsmax) -
                        8 * k * (ps - qsmax) * (ps + 6 * qsmax) -
                        (ps - qsmax)^2 * (ps + 6 * qsmax)
                    )
                ) -
                336 *
                Ek^4 *
                (
                    6 * k^2 - (ps - qsmax) * (ps + 4 * qsmax) -
                    k * (ps + 9 * qsmax)
                ) +
                12 *
                Ek^2 *
                (
                    72 * k^4 - 4 * k^3 * (8 * ps + 27 * qsmax) +
                    6 * k^2 * (7 * p0^2 - (ps - qsmax) * (3 * ps + 4 * qsmax)) +
                    (ps - qsmax) * (
                        -7 * p0^2 * (ps + 4 * qsmax) +
                        (ps - qsmax)^2 * (ps + 6 * qsmax)
                    ) +
                    k * (
                        6 * (ps - qsmax)^2 * (ps + 6 * qsmax) -
                        7 * p0^2 * (ps + 9 * qsmax)
                    )
                )
            ) *
            T *
            sinh(Ek / T)
        )
    ) / (40320 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T)
end

function delta2funcosthqs2(p0, ps, qsmax, k, Ek, T)
    (
        k *
        (2 * k - ps + qsmax)^3 *
        csch(Ek / (2 * T))^2 *
        (
            Ek *
            (4 * Ek^2 - p0^2) *
            (2 * k - ps + qsmax) *
            (
                36 * k^3 + 2 * k^2 * (ps - 36 * qsmax) -
                8 * k * (ps - qsmax) * (ps + 6 * qsmax) -
                (ps - qsmax)^2 * (ps + 6 * qsmax)
            ) +
            (
                -(
                    p0^2 *
                    (2 * k - ps + qsmax) *
                    (
                        36 * k^3 + 2 * k^2 * (ps - 36 * qsmax) -
                        8 * k * (ps - qsmax) * (ps + 6 * qsmax) -
                        (ps - qsmax)^2 * (ps + 6 * qsmax)
                    )
                ) -
                336 *
                Ek^4 *
                (
                    6 * k^2 - (ps - qsmax) * (ps + 4 * qsmax) -
                    k * (ps + 9 * qsmax)
                ) +
                12 *
                Ek^2 *
                (
                    72 * k^4 - 4 * k^3 * (8 * ps + 27 * qsmax) +
                    6 * k^2 * (7 * p0^2 - (ps - qsmax) * (3 * ps + 4 * qsmax)) +
                    (ps - qsmax) * (
                        -7 * p0^2 * (ps + 4 * qsmax) +
                        (ps - qsmax)^2 * (ps + 6 * qsmax)
                    ) +
                    k * (
                        6 * (ps - qsmax)^2 * (ps + 6 * qsmax) -
                        7 * p0^2 * (ps + 9 * qsmax)
                    )
                )
            ) *
            T *
            sinh(Ek / T)
        )
    ) / (40320 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T)
end

function delta2funcosthqs3(p0, ps, qsmax, k, Ek, T)
    (
        k *
        qsmax^3 *
        csch(Ek / (2 * T))^2 *
        (
            Ek *
            (-4 * Ek^2 + p0^2) *
            (
                35 * ps * (-2 * k + ps)^2 * (4 * k + ps) +
                42 * (-2 * k^2 + ps^2) * qsmax^2 +
                3 * qsmax^4
            ) +
            (
                12 *
                Ek^2 *
                (
                    35 *
                    (2 * k - ps) *
                    ps *
                    (8 * Ek^2 - 8 * k^2 - 2 * p0^2 + 2 * k * ps + ps^2) +
                    14 * (-4 * Ek^2 + 6 * k^2 + p0^2 - 3 * ps^2) * qsmax^2 -
                    3 * qsmax^4
                ) +
                p0^2 * (
                    35 * ps * (-2 * k + ps)^2 * (4 * k + ps) +
                    42 * (-2 * k^2 + ps^2) * qsmax^2 +
                    3 * qsmax^4
                )
            ) *
            T *
            sinh(Ek / T)
        )
    ) / (10080 * Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * ps * T)
end


function delta2funcosthqs4(p0, ps, qsmax, k, Ek, T)
    -1 / 20160 * (
        k *
        csch(Ek / (2 * T))^2 *
        (
            Ek *
            (4 * Ek^2 - p0^2) *
            (
                -ps^6 +
                21 * ps^4 * qsmax^2 +
                1120 * k^3 * qsmax^3 +
                105 * ps^2 * qsmax^4 +
                35 * qsmax^6 +
                42 * k^2 * (ps^4 - 10 * ps^2 * qsmax^2 - 15 * qsmax^4)
            ) -
            (
                336 *
                Ek^4 *
                (
                    ps^4 - 10 * ps^2 * qsmax^2 +
                    5 * (8 * k - 3 * qsmax) * qsmax^3
                ) +
                12 *
                Ek^2 *
                (
                    -7 * (6 * k^2 + p0^2) * ps^4 +
                    ps^6 +
                    7 * ps^2 * (60 * k^2 + 10 * p0^2 - 3 * ps^2) * qsmax^2 -
                    280 * k * (4 * k^2 + p0^2) * qsmax^3 +
                    105 * (6 * k^2 + p0^2 - ps^2) * qsmax^4 - 35 * qsmax^6
                ) +
                p0^2 * (
                    -ps^6 +
                    21 * ps^4 * qsmax^2 +
                    1120 * k^3 * qsmax^3 +
                    105 * ps^2 * qsmax^4 +
                    35 * qsmax^6 +
                    42 * k^2 * (ps^4 - 10 * ps^2 * qsmax^2 - 15 * qsmax^4)
                )
            ) *
            T *
            sinh(Ek / T)
        )
    ) / (Ek^3 * (-4 * Ek^2 + p0^2)^2 * pi^2 * T)
end
