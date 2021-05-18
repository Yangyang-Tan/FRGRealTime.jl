@doc raw"""
    flowpp(p0, ps, k, m, T,δ=0.02)

compute $-\frac{1}{\pi}\tilde{\partial_k}\Im I_{1, k}(p)$

To be noticed that, `flowpp` doesn't contains $\tilde{\partial_k}\mathcal{F}_4$,
we will consider it separately.

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
&\int_0^{qsmax}\!\!dq_s q_s^2\int_{-1}^{1}\!\!d\cos\theta \tilde{\partial_k}F_1\left(\sqrt{p_s^2+q_s^2+2*p_s*q_s\cos\theta}\right)\\
&=\int_0^{qsmax}\!\!dq_s q_s^2\int_{-1}^{1}\!\!d\cos\theta \tilde{\partial_k}F_1\left(\sqrt{p_s^2+q_s^2-2*p_s*q_s\cos\theta}\right)
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
&\int_0^{qsmax}\!\!dq_s q_s^2\int_{-1}^{1}\!\!d\cos\theta \tilde{\partial_k}F_2\left(\sqrt{p_s^2+q_s^2+2p_sq_s\cos\theta}\right)\\
&=\int_0^{qsmax}\!\!dq_s q_s^2\int_{-1}^{1}\!\!d\cos\theta \tilde{\partial_k}F_2\left(\sqrt{p_s^2+q_s^2-2*p_s*q_s\cos\theta}\right)
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


#integrate delta function in F1 , we integrate out qs, cos(θ) and k' this will be used in Im part calculation
function deltasumkfix(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun)
    #find the location of k0 where p0==2Epi(k,m)
    deltaf(x) = 2 * Epi(x, mfun(x)) - p0
    if deltaf(IRScale) * deltaf(UVScale) >= 0
        return 0.0
    else
        k0 = find_zero(deltaf, (IRScale, UVScale))
        # println(" k0=", k0)
        # k0 should lies between k~Λ, when k>k0 you will get 0
        if k > k0
            return 0.0
            #δ function only appears in p<2k so we have the following division
        elseif k <= k0
            if ps <= 2 * k0 - k
                if ps > k
                    return (
                        k^3 *
                        k0 *
                        (k^2 + 5 * ps * (-2 * k0 + ps)) *
                        coth(p0 / (4 * T)) *
                        lamfun(k0)^2
                    ) / (
                        30 *
                        p0 *
                        pi^2 *
                        ps *
                        abs(2 * k0 + derivative(mfun, k0))
                    )
                elseif ps <= k
                    return -1 / 120 * (
                        k0 *
                        (5 * k^3 * (-3 * k + 8 * k0) - 10 * k^2 * ps^2 + ps^4) *
                        coth(p0 / (4 * T)) *
                        lamfun(k0)^2
                    ) / (p0 * pi^2 * abs(2 * k0 + derivative(mfun, k0)))
                end
            elseif max(ps - k, 0.0) > 2 * k0
                return 0.0
            elseif 2*k < ps < 2 * k0 + k
                # println("locate end")
                return (k0*(k+2*k0-ps)^3*(6*k0^2+(k-ps)*(4*k+ps)-k0*(9*k+ps))*coth(p0/(4*T))*lamfun(k0)^2)/(240*p0*pi^2*ps*abs(2*k0+derivative(mfun,k0)))
            elseif 2 * k0 - k<ps < 2 * k0
                return (k0*(k+2*k0-ps)^3*(6*k0^2+(k-ps)*(4*k+ps)-k0*(9*k+ps))*coth(p0/(4*T))*lamfun(k0)^2)/(240*p0*pi^2*ps*abs(2*k0+derivative(mfun,k0)))
            end
        end
    end
end






# deltasum(k, m, T) =
#     -1 / 4 * (k^2 * coth(Epi(k, m) / (2 * T))) / (pi^2 * Epi(k, m)^3)


#integrate type1 delta function in F1 , we integrate out qs, cos(θ) and q0' this will be used in Re part calculation
function delta1_intcosthqs(p0, ps, qsmax, k, m, T)
    #δ function only appears in p<2k so we have the following division
    if ps >= 2 * k + qsmax
        return 0.0
        #for Principal Value integration we ommit k0-δ ~ k0+δ
    elseif 2 * k + qsmax > ps >= 2 * k
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return -1 / 120 * (
                k *
                (2 * k - ps + qsmax)^3 *
                (
                    6 * k^2 - (ps - qsmax) * (ps + 4 * qsmax) -
                    k * (ps + 9 * qsmax)
                ) *
                coth(Epi(k, m) / (2 * T)) *
                (p0 + Epi(k, m))
            ) / (p0^2 * pi^2 * ps * Epi(k, m) * (p0 + 2 * Epi(k, m)))
        else
            return (
                k *
                (2 * k - ps + qsmax)^3 *
                (
                    6 * k^2 - (ps - qsmax) * (ps + 4 * qsmax) -
                    k * (ps + 9 * qsmax)
                ) *
                coth(Epi(k, m) / (2 * T))
            ) / (120 * pi^2 * ps * Epi(k, m) * (-p0^2 + 4 * Epi(k, m)^2))
        end
    elseif 2 * k > ps >= 2 * k - qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return (
                k *
                qsmax^2 *
                (
                    10 * (-2 * k + ps)^2 * (k + ps) +
                    20 * (2 * k - ps) * ps * qsmax +
                    15 * (-k + ps) * qsmax^2 - 4 * qsmax^3
                ) *
                coth(Epi(k, m) / (2 * T)) *
                (p0 + Epi(k, m))
            ) / (120 * p0^2 * pi^2 * ps * Epi(k, m) * (p0 + 2 * Epi(k, m)))
        else
            return -1 / 120 * (
                k *
                qsmax^2 *
                (
                    -10 * (-2 * k + ps)^2 * (k + ps) +
                    20 * ps * (-2 * k + ps) * qsmax +
                    15 * (k - ps) * qsmax^2 +
                    4 * qsmax^3
                ) *
                coth(Epi(k, m) / (2 * T))
            ) / (pi^2 * ps * Epi(k, m) * (p0^2 - 4 * Epi(k, m)^2))
        end
    elseif 2 * k - qsmax > ps >= qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return (
                k *
                qsmax^3 *
                (10 * k * ps - 5 * ps^2 - qsmax^2) *
                coth(Epi(k, m) / (2 * T)) *
                (p0 + Epi(k, m))
            ) / (15 * p0^2 * pi^2 * ps * Epi(k, m) * (p0 + 2 * Epi(k, m)))
        else
            return (
                k *
                qsmax^3 *
                (5 * ps * (-2 * k + ps) + qsmax^2) *
                coth(Epi(k, m) / (2 * T))
            ) / (15 * pi^2 * ps * Epi(k, m) * (-p0^2 + 4 * Epi(k, m)^2))
        end
    elseif ps < qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return (
                k *
                (
                    ps^4 - 10 * ps^2 * qsmax^2 +
                    5 * (8 * k - 3 * qsmax) * qsmax^3
                ) *
                coth(Epi(k, m) / (2 * T)) *
                (p0 + Epi(k, m))
            ) / (60 * p0^2 * pi^2 * Epi(k, m) * (p0 + 2 * Epi(k, m)))
        else
            return (
                k *
                (
                    ps^4 - 10 * ps^2 * qsmax^2 +
                    5 * (8 * k - 3 * qsmax) * qsmax^3
                ) *
                coth(Epi(k, m) / (2 * T))
            ) / (60 * pi^2 * Epi(k, m) * (p0^2 - 4 * Epi(k, m)^2))
        end
    end
end



#integrate type2 delta function in F1 , we integrate out qs, cos(θ) and q0' this will be used in Re part calculation
function delta2_intcosthqs(p0, ps, qsmax, k, m, T,δk=0.02)
    #δ function only appears in p<2k so we have the following division
    Ek=Epi(k,m)
    if ps >= 2 * k + qsmax
        return 0.0
        #for Principal Value integration we ommit k0-δ ~ k0+δ
    elseif 2 * k + qsmax > ps >= 2 * k
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return 0.0
        else
            return delta2funcosthqs1(p0, ps, qsmax, k, Ek, T)
        end
    elseif 2 * k > ps >= 2 * k - qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return 0.0
        else
            return delta2funcosthqs2(p0, ps, qsmax, k, Ek, T)
        end
    elseif 2 * k - qsmax > ps >= qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return 0.0
        else
            return delta2funcosthqs3(p0, ps, qsmax, k, Ek, T)
        end
    elseif ps < qsmax
        if p0 + δk > 2 * Epi(k, m) > p0 - δk
            return 0.0
        else
            return delta2funcosthqs4(p0, ps, qsmax, k, Ek, T)
        end
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
