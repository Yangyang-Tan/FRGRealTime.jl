



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






# function PvdkF1Tildeps(p0, k, kprim, m, T)
#     if kprim > k / 2
#         if Epi(kprim + k, m) + Epi(kprim, m) > p0 > 2 * Epi(kprim, m)
#             quadgk_PV(
#                 x -> ppfunps(x, k, kprim, m, T),
#                 2 * Epi(kprim, m),
#                 p0,
#                 Epi(kprim + k, m) + Epi(kprim, m),
#                 rtol = error2,
#                 atol = error2,
#             ) +
#             deltasumps(p0, k, kprim, m, T) +
#             ppfunps_delta(p0, k, kprim, m, T)
#         else
#             hquadrature(
#                 x -> 2 * x * (-p0^2 + x^2)^-1 * ppfunps(x, k, kprim, m, T),
#                 2 * Epi(kprim, m),
#                 Epi(kprim + k, m) + Epi(kprim, m),
#                 rtol = error2,
#                 atol = error2,
#             )[1] +
#             deltasumps(p0, k, kprim, m, T) +
#             ppfunps_delta(p0, k, kprim, m, T)
#         end
#     elseif kprim <= k / 2
#         if Epi(kprim + k, m) + Epi(kprim, m) >
#            p0 >
#            Epi(kprim, m) + Epi(k - kprim, m)
#             quadgk_PV(
#                 x -> ppfunps(x, k, kprim, m, T),
#                 Epi(kprim, m) + Epi(k - kprim, m),
#                 p0,
#                 Epi(kprim + k, m) + Epi(kprim, m),
#                 rtol = error2,
#                 atol = error2,
#             )
#         else
#             hquadrature(
#                 x -> 2 * x * (-p0^2 + x^2)^-1 * ppfunps(x, k, kprim, m, T),
#                 Epi(kprim, m) + Epi(k - kprim, m),
#                 Epi(kprim + k, m) + Epi(kprim, m),
#                 rtol = error2,
#                 atol = error2,
#             )[1]
#         end
#     end
# end



#dkF1Tilde_intcostheqs has a delta function contribution in ps < 2 * k + qsmax brunch
#but we have seprated it

@doc raw"""
    PvdkF1Tildeps(p0, ps, qsmax, k, m, T; kwargs...)

Piecewise hints:

We only integrate the `p0'` domain which the  `flowpp_intcostheqs` is non zero.
"""
function PvdkF1Tildeps(p0, ps, qsmax, k, m, T; kwargs...)
    if ps >= 2 * k + qsmax
        quadgk_PV(
            x -> flowpp_intcostheqs(x, ps, qsmax, k, m, T),
            Epi(ps - k - qsmax, m) + Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) + Epi(k, m);
            kwargs...,
        )
    elseif ps < 2 * k + qsmax
        quadgk_PV(
            x -> flowpp_intcostheqs(x, ps, qsmax, k, m, T),
            2 * Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) + Epi(k, m);
            kwargs...,
        )
        +
        delta1_intcosthqs(p0, ps, qsmax, k, m, T)
        # +
        # delta2_intcosthqs(p0, ps, qsmax, k, m, T)
    end
end


function PvdkF2Tildeps(p0, ps, qsmax, k, m, T; kwargs...)
    if ps > 2 * k + qsmax
        quadgk_PV(
            x -> flowpm_intcostheqs(x, ps, qsmax, k, m, T),
            Epi(ps - k - qsmax, m) - Epi(k, m),
            p0,
            Epi(k + ps + qsmax, m) - Epi(k, m);
            kwargs...,
        )
    elseif ps < 2 * k + qsmax
        quadgk_PV(
            x -> flowpm_intcostheqs(x, ps, qsmax, k, m, T),
            0.0,
            p0,
            Epi(k + ps + qsmax, m) - Epi(k, m);
            kwargs...,
        )
    end
end


dkF1TildeintqsAll(p0, ps, qsmax, k, m, T; kwargs...) =
    PvdkF1Tildeps(abs(p0), ps, qsmax, k, m, T; kwargs...)




dkF2TildeintqsAll(p0, ps, qsmax, k, m, T; kwargs...) =
    PvdkF2Tildeps(abs(p0), ps, qsmax, k, m, T; kwargs...)







# test_dkF1Tilde_intcostheqs(10.0, 20.0, 60.0, 100.0, msgfun2(100.0), Tc)
#
#
# deltasum_intcosthqs(p0, ps, qsmax, k, m, T)=0.0
# deltasum(p0, ps, k, m2, T)=0.0
#
# plot(
#     p0 -> test_dkF1Tilde_intcostheqs(p0, 100.0, 10.0, 100.0, msgfun2(100.0), Tc),
#     475,
#     550.0,
# )
#
# plot(
#     p0 -> dkF1Tilde_intcostheqs(p0, 100.0, 20.0, 300.0, msgfun2(300.0), Tc),
#     600,
#     800.0,
# )


#
# F1TildeAll(p0, ps, k, m, T, p0UV) = F1Tilde(abs(p0), ps, k, m, T, p0UV)
# F2TildeAll(p0, ps, k, m, T) = F2Tilde(abs(p0), ps, k, m, T)

# dkF1TildeAll(p0::Float64, ps::Float64, k::Float64, m::Float64, T::Float64) = dkF1Tilde3(abs(p0), ps, k, m, T)
# dkF2TildeAll(p0::Float64, ps::Float64, k::Float64, m::Float64, T::Float64) = dkF2Tilde3(abs(p0), ps, k, m, T)


# dkF1TildeAll(p0, ps, k, m, T) = dkF1Tilde3(abs(p0), ps, k, m, T)
# dkF2TildeAll(p0, ps, k, m, T) = dkF2Tilde3(abs(p0), ps, k, m, T)
#
#
dkF1TildeAll(k, m, T) = dkF1Tilde3(k, m, T)
dkF2TildeAll(k, m, T) = dkF2Tilde3(k, m, T)
#
#
#
# dkF1TildeAll(p0, ps, k, m, T,p0UV)=central_fdm(5, 1)(x -> F1TildeAll(p0, ps, x, m, T, p0UV),k)
