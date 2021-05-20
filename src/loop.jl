"""
    deltafun(x,dϵ=0.02)
compute ``\\delta`` function with the approximation ``\\frac{\\epsilon}{\\pi(\\epsilon^2+x^2)}``
and ``\\epsilon`` is set to be ``0.02`` by default
"""
deltafun(x, dϵ = 0.02) = dϵ / (pi * (dϵ^2 + x^2))


"""
    Epi(k,m)

compute ``\\sqrt{(x^2+m)}``
"""
function Epi(q, m2)
    sqrt(q^2 + m2)
end


"""
    star1fun(qp, qm, ps, m, T)

compute ``-\\frac{1}{\\pi}\\mathcal{F}_{1}\\left(q_{+}, q_{-}, p, \\bar{m}_{\\pi, k}^{2}\\right)``
"""
function star1fun(qp, qm, ps, m, T)
    (2 * π)^-2 *
    (4 * ps)^-1 *
    (
        Epi(qp, m) - Epi(qm, m) +
        2 * T * log((1 - exp(-Epi(qp, m) / T)) / (1 - exp(-Epi(qm, m) / T)))
    )
end

@doc raw"""
    star1funpm(qp, qm, ps, m, T)

compute $-\frac{1}{\pi}\mathcal{F}_{1}^{\prime}\left(q_{+}, q_{-}, p, \bar{m}_{\pi, k}^{2}\right)$

Compared with $\mathcal{F}$, the $\mathcal{F'}$ has subtracted vacuum contributions. it will be used in $\mathrm{Im} I_2$
"""
function star1funpm(qp, qm, ps, m, T)
    (2 * π)^-2 *
    (4 * ps)^-1 *
    (2 * T * log((1 - exp(-Epi(qp, m) / T)) / (1 - exp(-Epi(qm, m) / T))))
end


@doc raw"""
    star2fun(qp, qm, ps, k, m, T)

compute $-\frac{1}{\pi}\mathcal{F}_{2}\left(q_{+}, q_{-}, k, p, \bar{m}_{\pi, k}^{2}\right)$
"""
function star2fun(qp, qm, ps, k, m, T)
    (2 * π)^-2 *
    (4 * ps)^-1 *
    (
        Epi(k, m)^-1 * coth(Epi(k, m) / (2 * T)) * (1 / 2) * (k^2 - qm^2) + (
            Epi(qp, m) - Epi(k, m) +
            2 * T * log((1 - exp(-Epi(qp, m) / T)) / (1 - exp(-Epi(k, m) / T)))
        )
    )
end


@doc raw"""
    star2funpm(qp, qm, ps, k, m, T)

compute $-\frac{1}{\pi}\mathcal{F}^{\prime}_{2}\left(q_{+}, q_{-}, k, p, \bar{m}_{\pi, k}^{2}\right)$
"""
function star2funpm(qp, qm, ps, k, m, T)
    (
        (k - qm) * (k + qm) * (-1 + coth(Epi(k, m) / (2 * T))) +
        4 *
        T *
        Epi(k, m) *
        log((-1 + exp(-(Epi(qp, m) / T))) / (-1 + exp(-(Epi(k, m) / T))))
    ) / (32 * pi^2 * ps * Epi(k, m))
end


@doc raw"""
    star3fun(qp, ps, k, m, T)

compute $-\frac{1}{\pi}\mathcal{F}_{3}\left(q_{+},k, p, \bar{m}_{\pi, k}^{2}\right)$
"""
function star3fun(qp, ps, k, m, T)
    (2 * π)^-2 *
    (4 * Epi(k, m))^-1 *
    qp *
    coth(Epi(qp, m) / (2 * T)) *
    (1 - (qp^2 + ps^2 - k^2) / (2 * qp * ps))
end


@doc raw"""
    star3funpm(qp, ps, k, m, T)

compute $-\frac{1}{\pi}\mathcal{F}^{\prime}_{3}\left(q_{+}, k, p, \bar{m}_{\pi, k}^{2}\right)$
"""
function star3funpm(qp, ps, k, m, T)
    ((k + ps - qp) * (k - ps + qp) * (-1 + coth(Epi(qp, m) / (2 * T)))) /
    (32 * pi^2 * ps * Epi(k, m))
end


@doc raw"""
    star4fun(p0, p, k, m, T)

compute $-\frac{1}{\pi}\mathcal{F}_{4}\left(p_{0}, k, p, \bar{m}_{\pi, k}^{2}\right)$
"""
function star4fun(p0, p, k, m, T, dϵ = 0.02)
    if k > p / 2
        return ((-2 * k + p)^2 * (4 * k + p) * coth(Epi(k, m) / (2 * T))) /
               (384 * pi^2 * Epi(k, m)^2) * deltafun(p0 - 2 * Epi(k, m), dϵ)
    else
        return 0.0
    end
end



@doc raw"""
    loopfunppfix(p0, ps, k, m, T)

compute $-\frac{1}{\pi}\Im I_{1, k}(p)$
"""
function loopfunppfix(p0, ps, k, m, T, dϵ = 0.02)
    star4fun(p0, ps, k, m, T, dϵ) + loopfunpp(p0, ps, k, m, T)
end



@doc raw"""
    loopfunpm(p0, ps, k, m, T)

compute $-\frac{1}{\pi}\Im I_{2, k}(p)$
"""
function loopfunpm(p0, ps, k, m, T)
    if m >= 0
        if k > ps
            if p0 > ps
                return 0.0
            elseif Epi(k + ps, m) - Epi(k, m) < p0 <= ps
                qm =
                    -ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1funpm(qp, qm, ps, m, T)
            elseif p0 <= Epi(k + ps, m) - Epi(k, m)
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            end
        elseif ps / 2 < k <= ps
            if p0 > ps
                return 0.0
            elseif Epi(k + ps, m) - Epi(k, m) < p0 <= ps
                qm =
                    -ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1funpm(qp, qm, ps, m, T)
            elseif Epi(ps, m) - Epi(k, m) < p0 <= Epi(k + ps, m) - Epi(k, m)
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif p0 <= Epi(ps, m) - Epi(k, m)
                qm = ps - sqrt((p0 + Epi(k, m))^2 - m)
                qp = sqrt((p0 + Epi(k, m))^2 - m)
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            end
        elseif k <= ps / 2
            if p0 > ps
                return 0.0
            elseif Epi(k + ps, m) - Epi(k, m) < p0 <= ps
                qm =
                    -ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1funpm(qp, qm, ps, m, T)
            elseif Epi(ps, m) - Epi(k, m) < p0 <= Epi(k + ps, m) - Epi(k, m)
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif Epi(k - ps, m) - Epi(k, m) < p0 <= Epi(ps, m) - Epi(k, m)
                qm = ps - sqrt((p0 + Epi(k, m))^2 - m)
                qp = sqrt((p0 + Epi(k, m))^2 - m)
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif p0 <= Epi(k - ps, m) - Epi(k, m)
                qm =
                    ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1funpm(qp, qm, ps, m, T)
            end
        end
    elseif m < 0
        if k > ps
            if p0 > Epi(k + ps, m) - Epi(k, m)
                return 0.0
            elseif ps < p0 <= Epi(k + ps, m) - Epi(k, m)
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                qmp =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (-4 * m + p0^2 - ps^2)) /
                    (2 * (p0^2 - ps^2))
                qpp = qmp + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T) -
                       star1funpm(qpp, qmp, ps, m, T)
            elseif p0 <= ps
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            end
        elseif ps / 2 < k <= ps
            if p0 > Epi(k + ps, m) - Epi(k, m)
                return 0.0
            elseif ps < p0 <= Epi(k + ps, m) - Epi(k, m)
                qm3 = sqrt(-m + (sqrt(k^2 + m) + p0)^2) - ps
                qm3p =
                    (
                        -ps +
                        sqrt(p0^2 * (-p0^2 + ps^2) * (4 * m - p0^2 + ps^2)) /
                        (p0^2 - ps^2)
                    ) / 2
                qp3 = sqrt(-m + (sqrt(k^2 + m) + p0)^2)
                qp3p =
                    (
                        ps +
                        sqrt(p0^2 * (-p0^2 + ps^2) * (4 * m - p0^2 + ps^2)) /
                        (p0^2 - ps^2)
                    ) / 2
                return -star1funpm(qp3p, qm3p, ps, m, T) +
                       star2funpm(qp3, qm3, ps, k, m, T) -
                       star3funpm(qp3, ps, k, m, T)
            elseif Epi(ps, m) - Epi(k, m) < p0 <= ps
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif p0 <= Epi(ps, m) - Epi(k, m)
                qm = ps - sqrt((p0 + Epi(k, m))^2 - m)
                qp = sqrt((p0 + Epi(k, m))^2 - m)
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            end
        elseif k <= ps / 2
            if p0 > Epi(k + ps, m) - Epi(k, m)
                return 0.0
            elseif ps < p0 <= Epi(k + ps, m) - Epi(k, m)
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                qmp =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (-4 * m + p0^2 - ps^2)) /
                    (2 * (p0^2 - ps^2))
                qpp = qmp + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T) -
                       star1funpm(qpp, qmp, ps, m, T)
            elseif Epi(ps, m) - Epi(k, m) < p0 <= ps
                qm = -ps + sqrt((p0 + Epi(k, m))^2 - m)
                qp = qm + ps
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif Epi(k - ps, m) - Epi(k, m) < p0 <= Epi(ps, m) - Epi(k, m)
                qm = ps - sqrt((p0 + Epi(k, m))^2 - m)
                qp = sqrt((p0 + Epi(k, m))^2 - m)
                return star2funpm(qp, qm, ps, k, m, T) -
                       star3funpm(qp, ps, k, m, T)
            elseif p0 <= Epi(k - ps, m) - Epi(k, m)
                qm =
                    ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1funpm(qp, qm, ps, m, T)
            end
        end
    end
end



function compensate_delta(ps, k, m, T)
    ((2 * k - ps) * coth(Epi(k, m) / (2 * T))) / (16 * pi^2 * Epi(k, m))
end

@doc raw"""
    loopfunpp(p0, ps, k, m, T)

compute $-\frac{1}{\pi}\Im I_{2, k}(p)$, but without $\mathcal{F}_4$
contribution
"""
function loopfunpp(p0, ps, k, m, T)
    if m >= 0
        if k > ps
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif 2 * Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif p0 < 2 * Epi(k, m)
                return 0.0
            end
        elseif ps / 2 < k <= ps
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif Epi(ps, m) + Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif 2 * Epi(k, m) <= p0 < Epi(ps, m) + Epi(k, m)
                qm = ps - sqrt((p0 - Epi(k, m))^2 - m)
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif p0 < 2 * Epi(k, m)
                return 0.0
            end
        elseif k <= ps / 2
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif Epi(ps, m) + Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif Epi(k, m) + Epi(ps - k, m) <= p0 < Epi(ps, m) + Epi(k, m)
                qm = ps - sqrt((p0 - Epi(k, m))^2 - m)
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif 2 * Epi(ps / 2, m) <= p0 < Epi(k, m) + Epi(ps - k, m)
                qm =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp =
                    ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star1fun(qp, qm, ps, m, T)
            elseif p0 < 2 * Epi(ps / 2, m)
                return 0.0
            end
        end
    elseif m < 0
        if k > ps
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif 2 * Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif p0 < 2 * Epi(k, m)
                return 0.0
            end
        elseif ps / 2 < k <= ps
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif Epi(ps, m) + Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif 2 * Epi(k, m) <= p0 < Epi(ps, m) + Epi(k, m)
                qm = ps - sqrt((p0 - Epi(k, m))^2 - m)
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif p0 < 2 * Epi(k, m)
                return 0.0
            end
        elseif k <= ps / 2
            if p0 >= Epi(k + ps, m) + Epi(k, m)
                qm =
                    -ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qp = qm + ps
                return star1fun(qp, qm, ps, m, T)
            elseif Epi(ps, m) + Epi(k, m) <= p0 < Epi(k + ps, m) + Epi(k, m)
                qm = sqrt((p0 - Epi(k, m))^2 - m) - ps
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif 2 * Epi(ps / 2, m) <= p0 < Epi(ps, m) + Epi(k, m)
                qm = ps - sqrt((p0 - Epi(k, m))^2 - m)
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                return star2fun(qp, qm, ps, k, m, T) + star3fun(qp, ps, k, m, T)
            elseif Epi(k, m) + Epi(ps - k, m) <= p0 < 2 * Epi(ps / 2, m)
                qm = ps - sqrt((p0 - Epi(k, m))^2 - m)
                qp = sqrt((p0 - Epi(k, m))^2 - m)
                qpp =
                    ps / 2 -
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                qmp =
                    ps / 2 +
                    sqrt(p0^2 * (p0^2 - ps^2) * (p0^2 - ps^2 - 4 * m)) /
                    (2 * (p0^2 - ps^2))
                return star2fun(qp, qm, ps, k, m, T) +
                       star3fun(qp, ps, k, m, T) - star1fun(qpp, qmp, ps, m, T)
            elseif p0 < Epi(k, m) + Epi(ps - k, m)
                return 0.0
            end
        end
    end
end
