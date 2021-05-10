@doc raw"""
    flowpp(p0, ps, k, m, T,δ=0.02)

compute $-\frac{1}{\pi}\tilde{\partial_k}\Im I_{1, k}(p)$

To be noticed that, `flowpp` doesn't constains $\tilde{\partial_k}\mathcal{F}_4$,
we will consider it seperately.

At $p_0=2E_{\pi,k}$, `flowpp` has a $\delta$ function contribution, we use a rectangle
function with width $\delta$ to approximate.
"""
function flowpp(p0, ps, k, m, T,δ=0.02)
    if k > ps / 2
        if p0 >= Epi(k + ps, m) + Epi(k, m)
            return 0.0
        elseif 2 * Epi(k, m) < p0 < Epi(k + ps, m) + Epi(k, m)
            return ppfun(p0, ps, k, Epi(k, m), T)
        elseif 2 * Epi(k - δ, m) <= p0 <= 2 * Epi(k, m)
            return -peak(p0, ps, m, T,δ)
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

@doc raw"""
    flowpm(p0, ps, k, m, T)

compute $-\frac{1}{\pi}\tilde{\partial_k}\Im I_{2, k}(p)$
"""
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



@doc raw"""
    flowpp_intcostheqs(p0, ps, qsmax, k, m, T)

compute
```math
\begin{aligned}
&\int_0^{qsmax}d\!\!q_s q_s^2\int_{-1}^{1}d\!\!\cos\theta \tilde{\partial_k}F_1\left(\sqrt{p_s^2+q_s^2+2*p_s*q_s}\right)\\
&=\int_0^{qsmax}d\!\!q_s q_s^2\int_{-1}^{1}d\!\!\cos\theta \tilde{\partial_k}F_1\left(\sqrt{p_s^2+q_s^2-2*p_s*q_s}\right)
\end{aligned}
```

"""
function flowpp_intcostheqs(p0, ps, qsmax, k, m, T)
    if p0 <= 2 * Epi(k, m)
        return 0.0
    elseif p0 > 2 * Epi(k, m)
        Ek2 = sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0))
        Ek = sqrt(k^2 + m)
        if ps >= k + Ek2 + qsmax
            return 0.0
        elseif k + Ek2 <= ps < k + Ek2 + qsmax
            return ppfuncostheps1(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 + k - qsmax <= ps < k + Ek2
            return ppfuncostheps2(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 <= ps < Ek2 + k - qsmax
            return ppfuncostheps3(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k + qsmax <= ps < Ek2
            return ppfuncostheps4(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < Ek2 - k + qsmax && qsmax <= ps
            return ppfuncostheps5(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < qsmax && ps >= -Ek2 + k + qsmax
            return ppfuncostheps6(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < qsmax && ps < -Ek2 + k + qsmax
            return ppfuncostheps7(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps <= (Ek2 - k) / 2 &&
               Ek2 - k + ps <= qsmax
            return ppfuncostheps8(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps <= (Ek2 - k) / 2 &&
               Ek2 - k + ps > qsmax
            return ppfuncostheps9(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps > (Ek2 - k) / 2 &&
               Ek2 - k + ps <= qsmax
            return ppfuncostheps10(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps > (Ek2 - k) / 2 &&
               Ek2 - k + ps > qsmax
            return ppfuncostheps11(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif ps <= Ek2 - k - qsmax
            return 0.0
        end
    end
end

@doc raw"""
    flowpm_intcostheqs(p0, ps, qsmax, k, m, T)

compute
```math
\begin{aligned}
&\int_0^{qsmax}d\!\!q_s q_s^2\int_{-1}^{1}d\!\!\cos\theta \tilde{\partial_k}F_2\left(\sqrt{p_s^2+q_s^2+2*p_s*q_s}\right)\\
&=\int_0^{qsmax}d\!\!q_s q_s^2\int_{-1}^{1}d\!\!\cos\theta \tilde{\partial_k}F_2\left(\sqrt{p_s^2+q_s^2-2*p_s*q_s}\right)
\end{aligned}
```
"""
function flowpm_intcostheqs(p0, ps, qsmax, k, m, T)
    Ek1 = sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0))
    Ek = sqrt(k^2 + m)
    # println(ps-qsmax,",",Ek1 - k)
    if ps > Ek1 + k + qsmax
        return 0.0
    elseif Ek1 + k < ps < Ek1 + k + qsmax
        #check
        return pmfuncostheps6(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 + k - qsmax < ps < Ek1 + k
        return pmfuncostheps8(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 < ps < Ek1 + k - qsmax
        return pmfuncostheps7(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 - k + qsmax <= ps < Ek1
        return pmfuncostheps9(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 - k <= ps < Ek1 - k + qsmax
        if p0 >= -sqrt(k^2 + m) + sqrt(m + (k + qsmax)^2)
            return pmfuncostheps10(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif p0 < -sqrt(k^2 + m) + sqrt(m + (k + qsmax)^2)
            if ps <= qsmax / 2
                return pmfuncostheps11(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif qsmax >= ps > qsmax / 2 && qsmax - ps >= Ek1 - k
                return pmfuncostheps11(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif qsmax >= ps > qsmax / 2 && qsmax - ps < Ek1 - k
                return pmfuncostheps12(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif ps > qsmax
                return pmfuncostheps12(p0, ps, qsmax, k, Ek, Ek1, T)
                # elseif ps > qsmax &&ps-qsmax > Ek1 - k
                #     println("su2")
                #     return pmfuncostheps13(p0, ps, qsmax, k, Ek, Ek1, T)
            end
        end
    elseif ps < Ek1 - k
        if ps <= -k + Ek1 - qsmax
            return 0.0
        elseif ps > -k + Ek1 - qsmax && ps >= qsmax
            return pmfuncostheps14(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps <= (Ek1 - k) / 2 &&
               Ek1 - k + ps <= qsmax
            return pmfuncostheps15(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps <= (Ek1 - k) / 2 &&
               Ek1 - k + ps > qsmax
            return pmfuncostheps16(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps > (Ek1 - k) / 2 &&
               Ek1 - k + ps <= qsmax
            return pmfuncostheps17(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps > (Ek1 - k) / 2 &&
               Ek1 - k + ps > qsmax
            return pmfuncostheps18(p0, ps, qsmax, k, Ek, Ek1, T)
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
