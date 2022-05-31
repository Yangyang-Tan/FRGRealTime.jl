function deltasum(p0, ps, k, m2, T,δk=0.1)
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




function dkF1Tilde(p0, ps, k, m, T;kwarg...)
    if k > ps / 2
        hquadrature(
            x -> 2 * x * (-p0^2 + x^2)^-1 * ppfun(x, ps, k, Epi(k, m), T),
            2 * Epi(k, m),
            Epi(k + ps, m) + Epi(k, m);kwarg...
        )[1] +
        deltasum(p0, ps, k, m, T) +
        F1tilde_delta(p0, ps, k, m, T)
    elseif k <= ps / 2
        hquadrature(
            x -> 2 * x * (-p0^2 + x^2)^-1 * ppfun(x, ps, k, Epi(k, m), T),
            2 * Epi(k, m),
            Epi(k + ps, m) + Epi(k, m);kwarg...
        )[1]
    end
end



function dkF2Tilde(p0, ps, k, m, T;kwarg...)
    if k > ps / 2
        if Epi(k + ps, m) - Epi(k, m) < p0
            return hquadrature(
                x -> 2 * x * (-p0^2 + x^2)^-1 * pmfun(x, ps, k, Epi(k, m), T),
                0.0,
                Epi(k + ps, m) - Epi(k, m);kwarg...
            )[1]+pmfun_zerofix(p0, ps, k, m, T,9)
        else
            return quadgk_PV(
                x -> pmfun(x, ps, k, Epi(k, m), T),
                0.0,
                p0,
                Epi(ps + k, m) - Epi(k, m);kwarg...
            )[1]+pmfun_zerofix(p0, ps, k, m, T,9)
        end
    elseif k <= ps / 2
        if Epi(k + ps, m) - Epi(k, m) < p0 || p0 < Epi(k - ps, m) - Epi(k, m)
            return hquadrature(
                x -> 2 * x * (-p0^2 + x^2)^-1 * pmfun(x, ps, k, Epi(k, m), T),
                Epi(k - ps, m) - Epi(k, m),
                Epi(k + ps, m) - Epi(k, m);kwarg...
            )[1]
        else
            return quadgk_PV(
                x -> pmfun(x, ps, k, Epi(k, m), T),
                Epi(k - ps, m) - Epi(k, m),
                p0,
                Epi(ps + k, m) - Epi(k, m);kwarg...
            )[1]
        end
    end
end


dkF1TildeAll(q0, qs, k, m, T;kwarg...)=dkF1Tilde(abs(q0), qs, k, m, T;kwarg...)
dkF2TildeAll(q0, qs, k, m, T;kwarg...)=dkF2Tilde(abs(q0), qs, k, m, T;kwarg...)

function dklamdabar_Zero(q0, qs, k, m, T, lam4pik, Npi;kwarg...)
        lam4pik^2 *
        (
            (2 + Npi) * (dkF1TildeAll(k, m, T)+dkF2TildeAll(k, m, T)) +
             6 * dkF1TildeAll(q0, qs, k, m, T;kwarg...) +
             6 * dkF2TildeAll(q0, qs, k, m, T;kwarg...)
        )
end

function lamdabar_Zero(q0,qs,k, T, Npi, mfun,lamfun,UVScale;kwarg...)
    -hquadrature(
        kprim -> dklamdabar_Zero(q0, qs, kprim, mfun(kprim), T, lamfun(kprim), Npi;kwarg...),
        k,
        UVScale;kwarg...
        # initdiv=200,
    )[1] + lamfun(UVScale)
end
