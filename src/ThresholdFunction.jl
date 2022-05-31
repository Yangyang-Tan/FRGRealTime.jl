#Compute Threshold Functions

module ThresholdFunctions

using ..MathFuns

@doc raw"""
    Fb1(k,m,T)
Compute the one bosonic propagator threshold function without external momentum dependence.
```math
\tilde{\partial_k}\int \frac{d^3q}{2\pi^3}\sum_{q_0}G_b(q_0)
```
"""
function Fb1(k, m2, T)
    1 / 12 *
    (k^4 * (2 * T * coth(Eb(k, m2) / (2 * T)) + Eb(k, m2) * csch(Eb(k, m2) / (2 * T))^2)) /
    (Eb(k, m2)^3 * pi^2 * T)
end

function Fb1(k::Float32, m2::Float32, T::Float32)
    1 / 12 *
    (k^4 * (2 * T * coth(Eb(k, m2) / (2 * T)) + Eb(k, m2) * csch(Eb(k, m2) / (2 * T))^2)) /
    (Eb(k, m2)^3 * pi^2.0f0 * T)
end



@doc raw"""
    Fb2(k,m,T)
Compute the two bosonic propagator threshold function without external momentum dependence.
```math
\tilde{\partial_k}\int \frac{d^3q}{2\pi^3}\sum_{q_0}G_b(q_0)G_b(q_0)
```
If we want to start from the momentum dependent case and take the zero momentum limit, we need to
take the $p_0\to 0$ first then take the $p_s\to 0$, they are noncommutative, or we need to take the $p_0$
goes faster to $0$ than $p_s$
"""
function Fb2(k, m2, T)
    (
        k^4 * (
            6 * T^2 * coth(Eb(k, m2) / (2 * T)) +
            csch(Eb(k, m2) / (2 * T))^2 *
            Eb(k, m2) *
            (3 * T + coth(Eb(k, m2) / (2 * T)) * Eb(k, m2))
        )
    ) / (48 * pi^2 * T^2 * Eb(k, m2)^5)
end

function Fb2(k::Float32, m2::Float32, T::Float32)
    (
        k^4 * (
            6 * T^2 * coth(Eb(k, m2) / (2 * T)) +
            csch(Eb(k, m2) / (2 * T))^2 *
            Eb(k, m2) *
            (3 * T + coth(Eb(k, m2) / (2 * T)) * Eb(k, m2))
        )
    ) / (48 * pi^2.0f0 * T^2 * Eb(k, m2)^5)
end


# function ReFb2positive(p0, k, m2, T, ϵ = 0.1, ϵ2 = 1e-14)
#     (
#         k^4 *
#         csch(Eb(k, m2) / (2 * T))^2 *
#         (
#             (
#                 ϵ2^2 *
#                 Eb(k, m2)^3 *
#                 (
#                     k^2 *
#                     coth(Eb(k, m2) / (2 * T)) *
#                     (k * ϵ2 + Eb(k, m2)^2) *
#                     (k^2 * ϵ2^2 - p0^2 * Eb(k, m2)^2) +
#                     T *
#                     Eb(k, m2) *
#                     (
#                         k^3 * ϵ2^2 * (2 * k + ϵ2) +
#                         p0^2 * (-2 * k^2 + Eb(k, m2)^2) * (k * ϵ2 + 2 * Eb(k, m2)^2)
#                     )
#                 )
#             ) / (T^2 * (k * ϵ2 + Eb(k, m2)^2)^2 * (k^2 * ϵ2^2 - p0^2 * Eb(k, m2)^2)^2) +
#             (2 * (p0 - 2 * Eb(k, m2))^2 * Eb(k, m2) * sinh(Eb(k, m2) / T)) /
#             (ϵ + (p0 - 2 * Eb(k, m2))^2)^2 +
#             (2 * Eb(k, m2) * sinh(Eb(k, m2) / T)) / (p0 + 2 * Eb(k, m2))^2 -
#             ((p0 - 2 * Eb(k, m2)) * (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))) /
#             (T * (ϵ + (p0 - 2 * Eb(k, m2))^2)) +
#             (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T)) / (p0 * T + 2 * T * Eb(k, m2))
#         )
#     ) / (48 * pi^2 * Eb(k, m2)^4)
# end


function ReFb2(p0, k, m2, T, ϵ = 0.1, ϵ2 = 1e-7)
    (k^4 * coth(Eb(k, m2) / (2 * T)) * PrincipalIntdeltadx(p0, 2 * Eb(k, m2),ϵ)) /
    (12 * pi^2 * Eb(k, m2)^3) +
    (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        PrincipalIntdelta(p0, 2 * Eb(k, m2),ϵ) *
        (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))
    ) / (48 * pi^2 * T * Eb(k, m2)^4)
    # +
    # (k^4*ϵ2^2*csch(Eb(k,m2)/(2*T))^2*(k^5*ϵ2^3*coth(Eb(k,m2)/(2*T))+Eb(k,m2)*(k^3*T*ϵ2*(-2*p0^2+ϵ2*(2*k+ϵ2))+Eb(k,m2)*(k^3*ϵ2*(-p0^2+k*ϵ2)*coth(Eb(k,m2)/(2*T))+p0^2*Eb(k,m2)*(k*T*(-4*k+ϵ2)-k^2*coth(Eb(k,m2)/(2*T))*Eb(k,m2)+2*T*Eb(k,m2)^2)))))/(48*pi^2*T^2*Eb(k,m2)*(k*ϵ2+Eb(k,m2)^2)^2*(k^2*ϵ2^2-p0^2*Eb(k,m2)^2)^2)
end


function ReFb2A(p0, k, m2, T, ϵ = 0.1, ϵ2 = 1e-7)
    (k^4 * coth(Eb(k, m2) / (2 * T)) * PrincipalIntdeltadx(p0, 2 * Eb(k, m2),ϵ)) /
    (12 * pi^2 * Eb(k, m2)^3)
    # +
    # (k^4*ϵ2^2*csch(Eb(k,m2)/(2*T))^2*(k^5*ϵ2^3*coth(Eb(k,m2)/(2*T))+Eb(k,m2)*(k^3*T*ϵ2*(-2*p0^2+ϵ2*(2*k+ϵ2))+Eb(k,m2)*(k^3*ϵ2*(-p0^2+k*ϵ2)*coth(Eb(k,m2)/(2*T))+p0^2*Eb(k,m2)*(k*T*(-4*k+ϵ2)-k^2*coth(Eb(k,m2)/(2*T))*Eb(k,m2)+2*T*Eb(k,m2)^2)))))/(48*pi^2*T^2*Eb(k,m2)*(k*ϵ2+Eb(k,m2)^2)^2*(k^2*ϵ2^2-p0^2*Eb(k,m2)^2)^2)
end

function ReFb2B(p0, k, m2, T, ϵ = 0.1, ϵ2 = 1e-7)
    (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        PrincipalIntdelta(p0, 2 * Eb(k, m2),ϵ) *
        (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))
    ) / (48 * pi^2 * T * Eb(k, m2)^4)
    # +
    # (k^4*ϵ2^2*csch(Eb(k,m2)/(2*T))^2*(k^5*ϵ2^3*coth(Eb(k,m2)/(2*T))+Eb(k,m2)*(k^3*T*ϵ2*(-2*p0^2+ϵ2*(2*k+ϵ2))+Eb(k,m2)*(k^3*ϵ2*(-p0^2+k*ϵ2)*coth(Eb(k,m2)/(2*T))+p0^2*Eb(k,m2)*(k*T*(-4*k+ϵ2)-k^2*coth(Eb(k,m2)/(2*T))*Eb(k,m2)+2*T*Eb(k,m2)^2)))))/(48*pi^2*T^2*Eb(k,m2)*(k*ϵ2+Eb(k,m2)^2)^2*(k^2*ϵ2^2-p0^2*Eb(k,m2)^2)^2)
end



# (
#             k^4 *
#             csch(Eb(k, m2) / (2 * T))^2 *
#             (
#                 -((p0 - ϵ) * (p0 + ϵ) * (p0^2 + ϵ^2)^2 * Eb(k, m2)) +
#                 4 * (3 * p0^4 - 2 * p0^2 * ϵ^2 + 3 * ϵ^4) * Eb(k, m2)^3 +
#                 48 * (-p0^2 + ϵ^2) * Eb(k, m2)^5 +
#                 64 * Eb(k, m2)^7 -
#                 T * (p0 - ϵ) * (p0 + ϵ) * (p0^2 + ϵ^2)^2 * sinh(Eb(k, m2) / T) +
#                 4 *
#                 T *
#                 (5 * p0^4 - 14 * p0^2 * ϵ^2 + 5 * ϵ^4) *
#                 Eb(k, m2)^2 *
#                 sinh(Eb(k, m2) / T) +
#                 112 * T * (-p0^2 + ϵ^2) * Eb(k, m2)^4 * sinh(Eb(k, m2) / T) +
#                 192 * T * Eb(k, m2)^6 * sinh(Eb(k, m2) / T)
#             )
#         ) / (
#             12 *
#             pi^2.0f0 *
#             T *
#             Eb(k, m2)^3 *
#             ((p0^2 + ϵ^2)^2 + 8 * Eb(k, m2)^2 * (-p0^2 + ϵ^2 + 2 * Eb(k, m2)^2))^2
#         )

function ReFb2(p0::Float32, k::Float32, m2::Float32, T::Float32, ϵ = 0.1f0, ϵ2 = 1f-8)
    # (
    #     k^4 *
    #     csch(Eb(k, m2) / (2 * T))^2 *
    #     (
    #         -((p0 - ϵ) * (p0 + ϵ) * (p0^2 + ϵ^2)^2 * Eb(k, m2)) +
    #         4 * (3 * p0^4 - 2 * p0^2 * ϵ^2 + 3 * ϵ^4) * Eb(k, m2)^3 +
    #         48 * (-p0^2 + ϵ^2) * Eb(k, m2)^5 +
    #         64 * Eb(k, m2)^7 -
    #         T * (p0 - ϵ) * (p0 + ϵ) * (p0^2 + ϵ^2)^2 * sinh(Eb(k, m2) / T) +
    #         4 *
    #         T *
    #         (5 * p0^4 - 14 * p0^2 * ϵ^2 + 5 * ϵ^4) *
    #         Eb(k, m2)^2 *
    #         sinh(Eb(k, m2) / T) +
    #         112 * T * (-p0^2 + ϵ^2) * Eb(k, m2)^4 * sinh(Eb(k, m2) / T) +
    #         192 * T * Eb(k, m2)^6 * sinh(Eb(k, m2) / T)
    #     )
    # ) / (
    #     12 *
    #     pi^2.0f0 *
    #     T *
    #     Eb(k, m2)^3 *
    #     ((p0^2 + ϵ^2)^2 + 8 * Eb(k, m2)^2 * (-p0^2 + ϵ^2 + 2 * Eb(k, m2)^2))^2
    # )
    (k^4 * coth(Eb(k, m2) / (2 * T)) * PrincipalIntdeltadx(p0, 2 * Eb(k, m2),ϵ)) /
    (12 * pi^2 * Eb(k, m2)^3) +
    (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        PrincipalIntdelta(p0, 2 * Eb(k, m2),ϵ) *
        (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))
    ) / (48 * pi^2.0f0 * T * Eb(k, m2)^4)
    # +
    # (k^4*ϵ2^2*csch(Eb(k,m2)/(2*T))^2*(k^5*ϵ2^3*coth(Eb(k,m2)/(2*T))+Eb(k,m2)*(k^3*T*ϵ2*(-2*p0^2+ϵ2*(2*k+ϵ2))+Eb(k,m2)*(k^3*ϵ2*(-p0^2+k*ϵ2)*coth(Eb(k,m2)/(2*T))+p0^2*Eb(k,m2)*(k*T*(-4*k+ϵ2)-k^2*coth(Eb(k,m2)/(2*T))*Eb(k,m2)+2*T*Eb(k,m2)^2)))))/(48*pi^2.0f0*T^2*Eb(k,m2)*(k*ϵ2+Eb(k,m2)^2)^2*(k^2*ϵ2^2-p0^2*Eb(k,m2)^2)^2)
end

#we add the ϵ by analitical continuation
function ReFb2positive_ac(ω, k, m2, T, ϵ = 0.1)
    -1 / 6 * (
        k^4 * (
            (ϵ^2 + (ω - 2 * Eb(k, m2))^2) *
            Eb(k, m2) *
            (ϵ^2 - ω^2 + 4 * Eb(k, m2)^2) *
            (ϵ^2 + (ω + 2 * Eb(k, m2))^2) +
            T *
            (
                (ϵ - ω) * (ϵ + ω) * (ϵ^2 + ω^2)^2 +
                4 * (5 * ϵ^4 - 14 * ϵ^2 * ω^2 + 5 * ω^4) * Eb(k, m2)^2 +
                112 * (ϵ - ω) * (ϵ + ω) * Eb(k, m2)^4 +
                192 * Eb(k, m2)^6
            ) *
            sinh(Eb(k, m2) / T)
        )
    ) / (
        pi^2 *
        T *
        (-1 + cosh(Eb(k, m2) / T)) *
        Eb(k, m2)^3 *
        ((ϵ^2 + ω^2)^2 + 8 * (ϵ - ω) * (ϵ + ω) * Eb(k, m2)^2 + 16 * Eb(k, m2)^4)^2
    )
end



function ReFb2positiveEu(p0, k, m2, T)
    -1 / 12 * (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        (
            Eb(k, m2) * (p0^2 + 4 * Eb(k, m2)^2) +
            T * (p0^2 + 12 * Eb(k, m2)^2) * sinh(Eb(k, m2) / T)
        )
    ) / (pi^2 * T * Eb(k, m2)^3 * (p0^2 + 4 * Eb(k, m2)^2)^2)
end



function ImFb2(p0, k, m2, T, ϵ = 0.1, ϵ2 = 1e-14)
    (
        k^4 *
        coth(Eb(k, m2) / (2 * T)) *
        (deltadxfun(p0 - 2 * Eb(k, m2), ϵ) + deltadxfun(p0 + 2 * Eb(k, m2), ϵ))
    ) / (12 * pi * Eb(k, m2)^3) +
    (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        (deltafun(p0 - 2 * Eb(k, m2), ϵ) - deltafun(p0 + 2 * Eb(k, m2), ϵ)) *
        (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))
    ) / (48 * pi * T * Eb(k, m2)^4)
end

function ImFb2(p0::Float32, k::Float32, m2::Float32, T::Float32, ϵ = 0.1f0, ϵ2 = 1f-14)
    (
        k^4 *
        coth(Eb(k, m2) / (2 * T)) *
        (deltadxfun(p0 - 2 * Eb(k, m2), ϵ) + deltadxfun(p0 + 2 * Eb(k, m2), ϵ))
    ) / (12.0f0 * pi * Eb(k, m2)^3) +
    (
        k^4 *
        csch(Eb(k, m2) / (2 * T))^2 *
        (deltafun(p0 - 2 * Eb(k, m2), ϵ) - deltafun(p0 + 2 * Eb(k, m2), ϵ)) *
        (Eb(k, m2) + 2 * T * sinh(Eb(k, m2) / T))
    ) / (48.0f0 * pi * T * Eb(k, m2)^4)
end



# function ImFb2(p0, k, m2, T, ϵ = 0.1)
#     sign(p0) * ImFb2positive(abs(p0), k, m2, T, ϵ)
# end

# function ReFb2(p0, k, m2, T, ϵ = 0.1)
#     ReFb2positive(abs(p0), k, m2, T, ϵ)
# end

function ReFb2_ac(p0, k, m2, T, ϵ = 0.1)
    if p0 >= 0.0
        ReFb2positive_ac(p0, k, m2, T, ϵ)
    elseif p0 < 0.0
        ReFb2positive_ac(-p0, k, m2, T, ϵ)
    end
end

function ReFb2Eu(p0, k, m2, T)
    if p0 >= 0.0
        ReFb2positiveEu(p0, k, m2, T)
    elseif p0 < 0.0
        ReFb2positiveEu(-p0, k, m2, T)
    end
end


end
