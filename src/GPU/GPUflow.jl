function flowpp_intcostheqs_gpu(p0, ps, qsmax, k, m, T)
    if p0 <= 2 * Epi(k, m)
        return 0.0f0
    elseif p0 > 2 * Epi(k, m)
        Ek2 = sqrt(k^2 + p0 * (-2 * sqrt(k^2 + m) + p0))
        Ek = sqrt(k^2 + m)
        if ps >= k + Ek2 + qsmax
            return 0.0f0
        elseif k + Ek2 <= ps < k + Ek2 + qsmax
            return ppfuncostheps1_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 + k - qsmax <= ps < k + Ek2
            return ppfuncostheps2_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 <= ps < Ek2 + k - qsmax
            return ppfuncostheps3_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k + qsmax <= ps < Ek2
            return ppfuncostheps4_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < Ek2 - k + qsmax && qsmax <= ps
            return ppfuncostheps5_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < qsmax && ps >= -Ek2 + k + qsmax
            return ppfuncostheps6_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k <= ps < qsmax && ps < -Ek2 + k + qsmax
            return ppfuncostheps7_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps <= (Ek2 - k) / 2 &&
               Ek2 - k + ps <= qsmax
            return ppfuncostheps8_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps <= (Ek2 - k) / 2 &&
               Ek2 - k + ps > qsmax
            return ppfuncostheps9_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps > (Ek2 - k) / 2 &&
               Ek2 - k + ps <= qsmax
            return ppfuncostheps10_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif Ek2 - k - qsmax < ps < Ek2 - k &&
               ps > (Ek2 - k) / 2 &&
               Ek2 - k + ps > qsmax
            return ppfuncostheps11_gpu(p0, ps, qsmax, k, Ek, Ek2, T)
        elseif ps <= Ek2 - k - qsmax
            return 0.0f0
        else
            return 0.0f0
        end
    else
        return 0.0f0
    end
end

function flowpm_intcostheqs_gpu(p0, ps, qsmax, k, m, T)
    Ek1 = sqrt(k^2 + p0 * (2 * sqrt(k^2 + m) + p0))
    Ek = sqrt(k^2 + m)
    # println(ps-qsmax,",",Ek1 - k)
    if ps > Ek1 + k + qsmax
        return 0.0f0
    elseif Ek1 + k < ps < Ek1 + k + qsmax
        #check
        return pmfuncostheps6_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 + k - qsmax < ps < Ek1 + k
        return pmfuncostheps8_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 < ps < Ek1 + k - qsmax
        return pmfuncostheps7_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 - k + qsmax <= ps < Ek1
        return pmfuncostheps9_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
    elseif Ek1 - k <= ps < Ek1 - k + qsmax
        if p0 >= -sqrt(k^2 + m) + sqrt(m + (k + qsmax)^2)
            return pmfuncostheps10_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif p0 < -sqrt(k^2 + m) + sqrt(m + (k + qsmax)^2)
            if ps <= qsmax / 2
                return pmfuncostheps11_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif qsmax >= ps > qsmax / 2 && qsmax - ps >= Ek1 - k
                return pmfuncostheps11_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif qsmax >= ps > qsmax / 2 && qsmax - ps < Ek1 - k
                return pmfuncostheps12_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
            elseif ps > qsmax
                return pmfuncostheps12_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
                # elseif ps > qsmax &&ps-qsmax > Ek1 - k
                #     println("su2")
                #     return pmfuncostheps13_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
            else
                return 0.0f0
            end
        else
            return 0.0f0
        end
    elseif ps < Ek1 - k
        if ps <= -k + Ek1 - qsmax
            return 0.0f0
        elseif ps > -k + Ek1 - qsmax && ps >= qsmax
            return pmfuncostheps14_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps <= (Ek1 - k) / 2 &&
               Ek1 - k + ps <= qsmax
            return pmfuncostheps15_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps <= (Ek1 - k) / 2 &&
               Ek1 - k + ps > qsmax
            return pmfuncostheps16_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps > (Ek1 - k) / 2 &&
               Ek1 - k + ps <= qsmax
            return pmfuncostheps17_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        elseif ps > -k + Ek1 - qsmax &&
               ps < qsmax &&
               ps > (Ek1 - k) / 2 &&
               Ek1 - k + ps > qsmax
            return pmfuncostheps18_gpu(p0, ps, qsmax, k, Ek, Ek1, T)
        else
            return 0.0f0
        end
    else
        return 0.0f0
    end
end
