################################################################################
#  The simplified computation of Im parts. we have integrated out qs & cos(θ)  #
################################################################################
@doc raw"""
    dkVIm(p0, ps, q0, k, m, T, Npi, lam4pik)

compute $\tilde{\partial_k}\mathrm{Im}V_k(q_0)$.

`dkVIm` only contains $V(q_0)$, no $V(-q_0)$

`dkVIm` contains type-1 delta function

`dkVIm` doesn't contains type-2 delta function
"""
function dkVIm(p0, ps, q0, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    π *
    (
        3 * (
            dkF1All(p0 - q0, ps, k, m, T) +
            dkF1All(p0 + q0, ps, k, m, T) +
            dkF2All(p0 - q0, ps, k, m, T) +
            dkF2All(p0 + q0, ps, k, m, T)
        )
        # +
        # (2 + Npi) * (
        #     dkF1All(1e-8 - 1e-14, 1e-8, k, m, T) +
        #     dkF2All(1e-8 - 1e-14, 1e-8, k, m, T)
        # )
    )
end








@doc raw"""
    dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)

compute $\int_0^{qsmax}dq_s qs^2\int_{-1}^{1}d\cos\theta \tilde{\partial_k}\mathrm{Im}V(q_0)$.

`dkVImintqs` only contains $V(q_0)$, for $-q_0$, we have $\int d\cos\theta V(q_0)=\int d\cos\theta V(-q_0)$,
so we need an extra $2$ at somewhere.

`dkVImintqs` doesn't have any delta function contribution, we include the type-1 delta function in `VImintqs_delta1` separately.

# Arguments
- `qsmax`: we integrate $q_s$ from $0$ to $k$, `qsmax` will set to `k` when we do the integration $dk'$, it should be distinguished from $k'$
- `m`: mass square, it will be $m(k')$ when we do the integration $dk'$.
- `lam4pik`: $\lambda_{4\pi}$, it will be $\lambda_{4\pi}(k')$ when we do the integration $dk'$ .
"""
function dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    π *
    (
        3 * (
            dkF1Allintqs(p0 - q0, ps, qsmax, k, m, T)
            +
            dkF1Allintqs(p0 + q0, ps, qsmax, k, m, T)
            +
            dkF2Allintqs(p0 - q0, ps, qsmax, k, m, T)
            +
            dkF2Allintqs(p0 + q0, ps, qsmax, k, m, T)
        )
    )
end


function dkVImintqs_gpu(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    π *
    (
        3 * (
            dkF1Allintqs_gpu(p0 - q0, ps, qsmax, k, m, T)
            +
            dkF1Allintqs_gpu(p0 + q0, ps, qsmax, k, m, T)
            +
            dkF2Allintqs_gpu(p0 - q0, ps, qsmax, k, m, T)
            +
            dkF2Allintqs_gpu(p0 + q0, ps, qsmax, k, m, T)
        )
    )
end



function dkVImintqs_F4(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 *
    (2 + Npi) *
    π *
    (
        3 * (
             dkF4Allintqs(p0 - q0, ps, qsmax, k, m, T)
             +
             dkF4Allintqs(p0 + q0, ps, qsmax, k, m, T)
        )
    )
end


function dkVImintqs_zeroMomentum(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)
    lam4pik^2 * (2 + Npi) * π * (Npi + 2) * 2 / 3 *
    qsmax^3 *
    (dkF1All(k, m, T) + dkF2All(k, m, T))
end



@doc raw"""
    VImintqs(p0, ps, k, T, Npi,UVScale,mfun::Function,lamfun::Function)

compute $\int_0^{k}dq_s qs^2\int_{-1}^{1}d\cos\theta \mathrm{Im}V(q_0,k)$.
In our code, we perform integration over `kprim`, `q0` & `qs` does not involved,
so `qs=k`, `q0=Epi(k, mfun(k))`.


`VImintqs` don't have any delta function contribution, we include the type-1 delta function in `VImintqs_delta1` separately.

# Arguments
- `mfun::Function`: $m^2(k)$, input from zero momentum result
- `lampifun::Function`: $\lambda_{4\pi}(k)$, input from zero momentum result.
"""
function VImintqs(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -hquadrature(
        kprim -> dkVImintqs(
            p0,
            ps,
            Epi(k, mfun(k)),
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim),
        ),
        k,
        UVScale,
        rtol = 1e-4,
        atol = 1e-4,
        initdiv=200,
    )[1]
end



function VImintqs(p0, ps, q0,k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -hquadrature(
        kprim -> dkVImintqs(
            p0,
            ps,
            q0,
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim),
        ),
        k,
        UVScale,
        rtol = 1e-4,
        atol = 1e-4,
        initdiv=200,
    )[1]
end


function VImintqs_zeroMomentum(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -hquadrature(
        kprim -> dkVImintqs_zeroMomentum(
            p0,
            ps,
            Epi(k, mfun(k)),
            k,
            kprim,
            mfun(kprim),
            T,
            Npi,
            lamfun(kprim),
        ),
        k,
        UVScale,
        rtol = 1e-2,
        atol = 1e-2,
    )[1]
end





@doc raw"""
    VImintqs_delta1(p0, ps, k, T, Npi,IRScale,UVScale, mfun::Spline1D, lamfun::Spline1D)

compute type-1 delta function contribution in $\mathrm{Im}V_k$, we do the triple integral analytically.
"""
function VImintqs_delta1(p0, ps, k, T, Npi,IRScale,UVScale, mfun, lamfun)
    -(2 + Npi) *
    π *
    3 *
    (
        deltasumkAll(p0 + Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun)
        +
        deltasumkAll(p0 - Epi(k, mfun(k)), ps, k, T, Npi,IRScale, UVScale, mfun, lamfun)
    )
end




@doc raw"""
    Coeffgamm2(k, T, mfun::Spline1D)

The coefficient before the two point function flow.
"""
function Coeffgamm2(k, T, mfun)
    (
        k * (
            -(coth(Epi(k, mfun(k)) / (2 * T)) / Epi(k, mfun(k))^3) -
            csch(Epi(k, mfun(k)) / (2 * T))^2 / (2 * T * Epi(k, mfun(k))^2)
        )
    ) / (16 * pi^2)
end








function propImintqs(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return -hquadrature(
        k ->
            2 *
            VImintqs(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun) *
            Coeffgamm2(k, T, mfun),
        IRScale,
        UVScale,
        rtol = 1e-4,
        atol = 1e-4,
        maxevals = 2000
    )[1]
end



function propImintqs_zeroMomentum(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return -hquadrature(
        k ->
            2 *
            VImintqs_zeroMomentum(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun) *
            Coeffgamm2(k, T, mfun),
        IRScale,
        UVScale,
        rtol = 1e-4,
        atol = 1e-4,
        maxevals = 2000
    )[1]
end

function propImintqs_delta1(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return -hquadrature(
        k ->
            2 *
            VImintqs_delta1(p0, ps, k, T, Npi, IRScale, UVScale, mfun, lamfun) *
            Coeffgamm2(k, T, mfun),
        IRScale,
        UVScale,
        rtol = 1e-6,
        atol = 1e-6,
        maxevals = 2000,
    )[1]
end

function propImintqs_All(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    return propImintqs(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun) +
           propImintqs_delta1(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)+propImintqs_zeroMomentum(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
end



function fastpropImintqs(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return -hcubature(
        x ->
            (x[1] - UVScale) *
            2 *
            dkVImintqs(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end



function fastpropImintqs_F4(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return -hcubature(
        x ->
            (x[1] - UVScale) *
            2 *
            dkVImintqs_F4(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end


function fastpropImintqs_zeroMomentum(
    p0,
    ps,
    T,
    IRScale,
    UVScale,
    Npi,
    mfun,
    lamfun;
    kwarg...,
)
    return -hcubature(
        x ->
            (x[1] - UVScale) *
            2 *
            dkVImintqs_zeroMomentum(
                p0,
                ps,
                Epi(x[1], mfun(x[1])),
                x[1],
                x[2] * (x[1] - UVScale) + UVScale,
                mfun(x[2] * (x[1] - UVScale) + UVScale),
                T,
                Npi,
                lamfun(x[2] * (x[1] - UVScale) + UVScale),
            ) *
            Coeffgamm2(x[1], T, mfun),
        [IRScale, 0.0],
        [UVScale, 1.0];
        kwarg...,
    )[1]
end


function fastpropImintqs_All(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun;kwarg...)
    return fastpropImintqs(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun;kwarg...) +
           propImintqs_delta1(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)+fastpropImintqs_F4(
               p0,
               ps,
               T,
               IRScale,
               UVScale,
               Npi,
               mfun,
               lamfun;
               kwarg...,
           )
end







function propImintqs_F4simple(p0, ps, T, IRScale, UVScale, Npi, mfun, lamfun)
    if p0^2<=mfun(IRScale)
        return 0.0
    elseif p0^2>mfun(IRScale)
        return (5*(2+Npi)*(-mfun(IRScale)+p0^2)^3*csch(p0/(2*T))^4*lamfun(sqrt(-mfun(IRScale)+p0^2))^2*sinh(p0/T)*(p0+T*sinh(p0/T)))/(4096*p0^4*pi^2*pi*T)
    end
end









#integrate delta function in F1 , we integrate out qs, cos(θ) and k' this will be used in Im part calculation
function deltasumImprop1(p0, ps, T,IRScale,UVScale, Npi, mfun, lamfun)
    #find the location of k0 where p0==2Epi(k,m)
    deltaf(x) = 3 * Epi(x, mfun(x)) - p0
    if deltaf(IRScale) * deltaf(UVScale) >= 0
        return 0.0
    else
        k0 = find_zero(deltaf, (IRScale, UVScale))
        # println(" k0=", k0)
        # k0 should lies between k~Λ, when k>k0 you will get 0
        if IRScale > k0
            return 0.0
            #δ function only appears in p<2k so we have the following division
        elseif IRScale <= k0
            if ps < k0
                return (
                    3 *
                    k0^2 *
                    (2 + Npi) *
                    coth(p0 / (6 * T)) *
                    csch(p0 / (6 * T))^2 *
                    lamfun(k0)^2 *
                    (
                        3075 * k0^4 - 1900 * k0^2 * ps^2 +
                        137 * ps^4 +
                        60 * ps^2 * (-20 * k0^2 + ps^2) * log(k0 / ps)
                    ) *
                    sign(2 * k0 + derivative(mfun, k0)) *
                    (p0 + 3 * T * sinh(p0 / (3 * T)))
                ) / (
                    25600 * p0^4 * pi^3 * T * (2 * k0 + derivative(mfun, k0))^2
                )
            elseif k0 < ps <= 2 * k0
                return (
                    3 *
                    k0^2 *
                    (2 + Npi) *
                    coth(p0 / (6 * T)) *
                    csch(p0 / (6 * T))^2 *
                    lamfun(k0)^2 *
                    (
                        (3 * k0 - ps) * (
                            132 * k0^4 + 3099 * k0^3 * ps - 2047 * k0^2 * ps^2 -
                            9 * k0 * ps^3 + 137 * ps^4
                        ) +
                        60 *
                        (2 * k0 - ps)^3 *
                        (4 * k0^2 + 6 * k0 * ps + ps^2) *
                        log(k0 / (2 * k0 - ps))
                    ) *
                    sign(2 * k0 + derivative(mfun, k0)) *
                    (p0 + 3 * T * sinh(p0 / (3 * T)))
                ) / (
                    51200 *
                    p0^4 *
                    pi^3 *
                    ps *
                    T *
                    (2 * k0 + derivative(mfun, k0))^2
                )
            elseif 2*k0 < ps <= 3 * k0
                return (
                    -3 *
                    k0^2 *
                    (2 + Npi) *
                    coth(p0 / (6 * T)) *
                    csch(p0 / (6 * T))^2 *
                    lamfun(k0)^2 *
                    (
                        -(
                            (3 * k0 - ps) * (
                                132 * k0^4 + 3099 * k0^3 * ps -
                                2047 * k0^2 * ps^2 - 9 * k0 * ps^3 + 137 * ps^4
                            )
                        ) +
                        60 *
                        (2 * k0 - ps)^3 *
                        (4 * k0^2 + 6 * k0 * ps + ps^2) *
                        log(-2 + ps / k0)
                    ) *
                    sign(2 * k0 + derivative(mfun, k0)) *
                    (p0 + 3 * T * sinh(p0 / (3 * T)))
                ) / (
                    51200 *
                    p0^4 *
                    pi^3 *
                    ps *
                    T *
                    (2 * k0 + derivative(mfun, k0))^2
                )
            elseif ps > 3 * k0
                return 0.0
            end
        end
    end
end


function deltasumImprop2(p0, ps, T,IRScale,UVScale, Npi, mfun, lamfun)
    #find the location of k0 where p0==2Epi(k,m)
    deltaf(x) = Epi(x, mfun(x)) - p0
    if deltaf(IRScale) * deltaf(UVScale) >= 0
        return 0.0
    else
        k0 = find_zero(deltaf, (IRScale, UVScale))
        # println(" k0=", k0)
        # k0 should lies between k~Λ, when k>k0 you will get 0
        if IRScale > k0
            return 0.0
            #δ function only appears in p<2k so we have the following division
        elseif IRScale <= k0
            if ps < k0
                return (k0^2*(2+Npi)*csch(p0/(2*T))^4*lamfun(k0)^2*(3075*k0^4-1900*k0^2*ps^2+137*ps^4+60*ps^2*(-20*k0^2+ps^2)*log(k0/ps))*sign(2*k0+derivative(mfun,k0))*sinh(p0/T)*(p0+T*sinh(p0/T)))/(153600*p0^4*pi^3*T*(2*k0+derivative(mfun,k0))^2)
            elseif k0 < ps <= 2 * k0
                return (k0^2*(2+Npi)*csch(p0/(2*T))^4*lamfun(k0)^2*((3*k0-ps)*(132*k0^4+3099*k0^3*ps-2047*k0^2*ps^2-9*k0*ps^3+137*ps^4)+60*(2*k0-ps)^3*(4*k0^2+6*k0*ps+ps^2)*log(k0/(2*k0-ps)))*sign(2*k0+derivative(mfun,k0))*sinh(p0/T)*(p0+T*sinh(p0/T)))/(307200*p0^4*pi^3*ps*T*(2*k0+derivative(mfun,k0))^2)
            elseif 2*k0 < ps <= 3 * k0
                return (k0^2*(2+Npi)*csch(p0/(2*T))^4*lamfun(k0)^2*((3*k0-ps)*(132*k0^4+3099*k0^3*ps-2047*k0^2*ps^2-9*k0*ps^3+137*ps^4)+60*(2*k0-ps)^3*(4*k0^2+6*k0*ps+ps^2)*log(k0/(2*k0-ps)))*sign(2*k0+derivative(mfun,k0))*sinh(p0/T)*(p0+T*sinh(p0/T)))/(307200*p0^4*pi^3*ps*T*(2*k0+derivative(mfun,k0))^2)
            elseif ps > 3 * k0
                return 0.0
            end
        end
    end
end
