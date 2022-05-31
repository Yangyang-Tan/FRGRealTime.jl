

#Zero Momentum

#without a maxiter

function dkF1Tilde3(k, m, T)
    -1 / 48 * (
        k^4 * (
            6 * T * coth(Epi(k, m) / (2 * T)) +
            csch(Epi(k, m) / (2 * T))^2 * Epi(k, m)
        )
    ) / (pi^2 * T * Epi(k, m)^5)
end

function dkF2Tilde3(k, m, T)
    -1 / 48 * (
        k^4 *
        csch(Epi(k, m) / (2 * T))^2 *
        (2 * T + coth(Epi(k, m) / (2 * T)) * Epi(k, m))
    ) / (pi^2 * T^2 * Epi(k, m)^4)
end




function PvdkF1Tildeps(p0, qsmax, k, m, T; kwargs...)
    if k > qsmax / 2
        return quadgk_PV(
            x -> ppfunps(x, qsmax, k, m, T),
            2 * Epi(k, m),
            p0,
            Epi(k + qsmax, m) + Epi(k, m);
            kwargs...,
        )
    elseif k <= qsmax / 2
        return quadgk_PV(
            x -> ppfunps(x, qsmax, k, m, T),
            Epi(k, m) + Epi(qsmax - k, m),
            p0,
            Epi(k + qsmax, m) + Epi(k, m);
            kwargs...,
        )
    end
end



function dkF1Tildeps_delta1(p0, qsmax, k, m, T)
    if k > qsmax / 2
        deltasumps(p0, qsmax, k, m, T)
    elseif k <= qsmax / 2
        0.0
    end
end

function dkF1Tildeps_delta2(p0, qsmax, k, m, T)
    if k > qsmax / 2
        ppfunps_delta(p0, qsmax, k, m, T)
    elseif k <= qsmax / 2
        0.0
    end
end




function PvdkF2Tildeps(p0, qsmax, k, m, T; kwargs...)
    return quadgk_PV(
        x -> pmfunps(x, qsmax, k, m, T),
        0.0,
        p0,
        Epi(k + qsmax, m) - Epi(k, m);
        kwargs...,
    )[1] + pmfunps_zerofix(p0, qsmax, k, m, T, 9)
end








#dkF1Tilde_intcostheqs has a delta function contribution in ps < 2 * k + qsmax brunch
#but we have seprated it

@doc raw"""
    PvdkF1Tildeps(p0, ps, qsmax, k, m, T; kwargs...)

Piecewise hints:

We only integrate the `p0'` domain which the  `flowpp_intcostheqs` is non zero.
"""
function PvdkF1Tildeps(p0, ps, qsmax, k, m, T; kwargs...)
    if ps >= 2 * k + qsmax
        return quadgk_PV2(
            x -> flowpp_intcostheqs(x, ps, qsmax, k, m, T),
            Epi(ps - k - qsmax, m) + Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) + Epi(k, m);
            kwargs...,
        )
    elseif ps < 2 * k + qsmax
        return quadgk_PV3(
            x -> flowpp_intcostheqs(x, ps, qsmax, k, m, T),
            2 * Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) + Epi(k, m);
            kwargs...,
        )
    end
end



function PvdkF1Tildeps_Compensate(p0, ps, qsmax, k, m, T)
    if ps >= 2 * k + qsmax
        return 0.0
    elseif ps < 2 * k + qsmax
        if Epi(k + ps + qsmax, m) + Epi(k, m) > p0 > 2 * Epi(k, m)
            return flowpp_intcostheqs(p0, ps, qsmax, k, m, T) * log(
                abs(
                    ((Epi(k + ps + qsmax, m) + Epi(k, m))^2 - p0^2) /
                    ((2 * Epi(k, m))^2 - p0^2),
                ),
            )
        else
            return flowpp_intcostheqs(nextfloat(2 * Epi(k, m)), ps, qsmax, k, m, T) * log(
                abs(
                    ((Epi(k + ps + qsmax, m) + Epi(k, m))^2 - p0^2) /
                    ((2 * Epi(k, m))^2 - p0^2),
                ),
            )
        end
    end
end


function dkF1Tildeps_delta1(p0, ps, qsmax, k, m, T)
    if ps >= 2 * k + qsmax
        0.0
    elseif ps < 2 * k + qsmax
        delta1_intcosthqs(p0, ps, qsmax, k, m, T,0.5)
    end
end


function dkF1Tildeps_delta2(p0, ps, qsmax, k, m, T)
    if ps >= 2 * k + qsmax
        0.0
    elseif ps < 2 * k + qsmax
        delta2Smooth_intcosthqs(p0, ps, qsmax, k, m, T)
    end
end


function PvdkF2Tildeps(p0, ps, qsmax, k, m, T; kwargs...)
    if ps > 2 * k + qsmax
        return quadgk_PV2(
            x -> flowpm_intcostheqs(x, ps, qsmax, k, m, T),
            Epi(ps - k - qsmax, m) - Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) - Epi(k, m);
            kwargs...,
        )
    elseif ps < 2 * k + qsmax
        return quadgk_PV2(
            x -> flowpm_intcostheqs(x, ps, qsmax, k, m, T),
            0.0,
            p0,
            Epi(k + ps + qsmax, m) - Epi(k, m);
            kwargs...,
        )
    end
end



@doc raw"""
    dkF1TildeintqsAll(p0, qsmax, k, m, T; kwargs...)

`costh` is not integrated we need an extra `2` at somewhere
"""
function dkF1TildeintqsAll(p0, qsmax, k, m, T; kwargs...)
    return PvdkF1Tildeps(abs(p0), qsmax, k, m, T; kwargs...) +
           dkF1Tildeps_delta1(abs(p0), qsmax, k, m, T) +
           dkF1Tildeps_delta2(abs(p0), qsmax, k, m, T)
end
@doc raw"""
    dkF2TildeintqsAll(p0, qsmax, k, m, T; kwargs...)

`costh` is not integrated we need an extra `2` at somewhere
"""
function dkF2TildeintqsAll(p0, qsmax, k, m, T; kwargs...)
    PvdkF2Tildeps(abs(p0), qsmax, k, m, T; kwargs...)
end

function dkF1TildeintqsAll(p0, ps, qsmax, k, m, T; kwargs...)
    return PvdkF1Tildeps(abs(p0), ps, qsmax, k, m, T; kwargs...)
    # +
           # dkF1Tildeps_delta1(abs(p0), ps, qsmax, k, m, T) +
           # dkF1Tildeps_delta2(abs(p0), ps, qsmax, k, m, T)    # +
end

function dkF1TildeintqsAll_delta1(p0, ps, qsmax, k, m, T)
    return dkF1Tildeps_delta1(abs(p0), ps, qsmax, k, m, T)
end
function dkF1TildeintqsAll_delta2(p0, ps, qsmax, k, m, T)
    return dkF1Tildeps_delta2(abs(p0), ps, qsmax, k, m, T)
end
function dkF1TildeintqsAll_Compensate(p0, ps, qsmax, k, m, T)
    return PvdkF1Tildeps_Compensate(abs(p0), ps, qsmax, k, m, T)
end

function dkF2TildeintqsAll(p0, ps, qsmax, k, m, T; kwargs...)
    PvdkF2Tildeps(abs(p0), ps, qsmax, k, m, T; kwargs...)
end



function dkF1TildeAll(k, m, T)
    dkF1Tilde3(k, m, T)
end
function dkF2TildeAll(k, m, T)
    dkF2Tilde3(k, m, T)
end
