@doc raw"""
    flowpp(p0, ps, k, m, T,δ=0.02)

compute $-\frac{1}{\pi}\tilde{\partial_k}\Im I_{1, k}(p)$

To be noticed that, `flowpp` doesn't constains $\tilde{\partial_k}\mathcal{F}_4$,
we will consider it seperately.

At $p_0=2E_{pi,k}$, `flowpp` has a $\delta$ function contribution, we use a rectangle
function with width $\delta$ to approximate.
"""
function flowpp(p0, ps, k, m, T,δ=0.02)
    if k > ps / 2
        if p0 >= Epi(k + ps, m) + Epi(k, m)
            return 0.0
        elseif 2 * Epi(k, m) < p0 < Epi(k + ps, m) + Epi(k, m)
            return ppfun(p0, ps, k, Epi(k, m), T)
        elseif 2 * Epi(k - δ, m) <= p0 <= 2 * Epi(k, m)
            return -peak(p0, ps, m, T)
        elseif p0 < 2 * Epi(k - δ, m)
            return 0.0
        end
    elseif k <= ps / 2
        if p0 >= Epi(k + ps, m) + Epi(k, m)
            return 0.0
        elseif Epi(ps, m) + Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
            return ppfun(p0, ps, k, Epi(k, m), T)
        elseif Epi(k, m) + Epi(ps - k, m) <= p0 < Epi(ps, m) + Epi(k, m)
            return ppfun(p0, ps, k, Epi(k, m), T)
        elseif p0 < Epi(k, m) + Epi(ps - k, m)
            return 0.0
        end
    end
end


integroconst(p0, ps) = p0 / (16 * pi^2 * ps)

deltasumpole(p0, ps, k, m2, T) =
    (k * (-2 * sqrt(k^2) + ps) * coth(sqrt(k^2 + m2) / (2 * T))) /
    (2 * sqrt(k^2 + m2) * (4 * (k^2 + m2) - p0^2) * pi^2)

# deltasum(k, m, T) =
#     -1 / 4 * (k^2 * coth(Epi(k, m) / (2 * T))) / (pi^2 * Epi(k, m)^3)

function deltasum(p0, ps, k, m2, T)
    if p0 + δk > 2 * Epi(k, m2) > p0 - δk
        return (k * (-2 * k + ps) * coth(Epi(k, m2) / (2 * T))) /
               (8 * Epi(k, m2)^2 * (2 * Epi(k, m2) + p0) * pi^2)
    else
        return (k * (-2 * k + ps) * coth(sqrt(k^2 + m2) / (2 * T))) /
               (2 * sqrt(k^2 + m2) * (4 * (k^2 + m2) - p0^2) * pi^2)
    end
end


function F1tilde_delta(p0, p, k, m, T)
    (
        k *
        (2 * k - p) *
        csch(sqrt(k^2 + m) / (2 * T))^2 *
        (
            sqrt(k^2 + m) *
            (2 * k - p) *
            (4 * k + p) *
            (-4 * (k^2 + m) + p0^2) +
            (
                (2 * k - p) * (4 * k + p) * p0^2 +
                12 * (k^2 + m) * (8 * m + 2 * k * p + p^2 - 2 * p0^2)
            ) *
            T *
            sinh(sqrt(k^2 + m) / T)
        )
    ) / (192 * (k^2 + m)^(3 / 2) * (-4 * (k^2 + m) + p0^2)^2 * pi^2 * T)
end


function deltasumps(p0, psmax, k, m2, T)
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


function pmfun_zerofix(p0, p, k, m, T,a)
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

function pmfunps_zerofix(p0, psmax, k, m, T,a)
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






function peak(p0, ps, m2, T)
    ((sqrt(-4 * m2 + p0^2) - ps) * coth(p0 / (4 * T))) / (8 * p0 * δ * pi^2)
end




function flowpm(p0, ps, k, m, T)
    if k > ps / 2
        if p0 > Epi(k + ps, m) - Epi(k, m)
            return 0.0
        elseif p0 <= Epi(k + ps, m) - Epi(k, m)
            return pmfun(p0, ps, k, Epi(k, m), T)
        end
    elseif k <= ps / 2
        if p0 > Epi(k + ps, m) - Epi(k, m)
            return 0.0
        elseif Epi(k - ps, m) - Epi(k, m) < p0 <= Epi(k + ps, m) - Epi(k, m)
            return pmfun(p0, ps, k, Epi(k, m), T)
        elseif p0 <= Epi(k - ps, m) - Epi(k, m)
            return 0.0
        end
    end
end


function flowpm(p0, k, m, T)
    if p0 > 0.0
        return 0.0
    elseif p0 <= 0
        return pmfun(p0, ps, k, Epi(k, m), T)
    end
end

function flowpm_ps(p0, ps, k, m, T)
    if 0.0 <= ps <= -k + sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2)
        return 0.0
    elseif ps > -k + sqrt(k^2 + 2 * sqrt(k^2 + m) * p0 + p0^2)
        return pmfun(p0, ps, k, Epi(k, m), T)
    end
end

function flowpm_intps(p0, psmax, k, m, T)
    if 0 <= p0 <= Epi(k + psmax, m) - Epi(k, m)
        return pmfunps(p0, k, m, T)
    elseif p0 >= Epi(k + psmax, m) - Epi(k, m)
        return 0.0
    end
end










function flowpp2(p0, k, m, T)
    if p0 >= 2 * Epi(k, m)
        return 0.0
    elseif 2 * Epi(k - δ, m) <= p0 < 2 * Epi(k, m)
        return -peak(p0, ps, m, T)
    elseif p0 < 2 * Epi(k - δ, m)
        return 0.0
    end
end





function flowpp_ps(p0, ps, k, m, T)
    if 0.0 <= ps <= -k + sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)
        return 0.0
    elseif ps > -k + sqrt(k^2 - 2 * sqrt(k^2 + m) * p0 + p0^2)
        if p0 <= 2 * Epi(k, m)
            return 0.0
        elseif p0 > 2 * Epi(k, m)
            return ppfun(p0, ps, k, Epi(k, m), T)
        end
    end
end

#integrate qs from 0 to qsmax
function flowpp_intps(p0, psmax, k, m, T)
    if p0 <= 2 * Epi(k, m)
        return 0.0
    elseif p0 > 2 * Epi(k, m)
        if p0 >= sqrt(k^2 + m) + sqrt(m + (k + psmax)^2)
            return 0.0
        elseif p0 < sqrt(k^2 + m) + sqrt(m + (k + psmax)^2)
            return ppfunps(p0, psmax, k, m, T)
        end
    end
end

function flowpp_intps2(p0, psmax, k, m, T)
    if p0 <= 2 * Epi(k, m)
        return 0.0
    elseif Epi(k, m) + Epi(k + psmax, m) > p0 > 2 * Epi(k, m)
        return ppfunps(p0, psmax, k, m, T)
    elseif p0 >= Epi(k, m) + Epi(k + psmax, m)
        return 0.0
    end
end
