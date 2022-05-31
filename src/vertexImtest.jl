@doc raw"""
    Coeffgamm2Simple(k, T, m)

The coefficient before the two point function flow.
"""
function Coeffgamm2Simple(k, T, m)
    (
        k * (
            -(coth(Epi(k, m) / (2 * T)) / Epi(k, m)^3) -
            csch(Epi(k, m) / (2 * T))^2 / (2 * T * Epi(k, m)^2)
        )
    ) / (16 * pi^2)
end


@doc raw"""
    VImSimple(p0, ps, q0, k, m, T, Npi, lam4pik)

compute $\mathrm{Im}V_k(q_0)$, the `k` dependence of $\lambda_{4\pi}$ and $m$ are neglected.

`VImSimple` only contains $V(q_0)$, no $V(-q_0)$

`VImSimple` contains type-1 delta function

`VImSimple` doesn't contains type-2 delta function
"""
function VImSimple(p0, ps, q0, k, m, T, Npi, lam4pik, UVScale)
    lam4pik^2 *
    (2 + Npi) *
    π *
    (
        3 * (
            F1All(p0 - q0, ps, k, m, T) - F1All(p0 - q0, ps, UVScale, m, T) +
            F1All(p0 + q0, ps, k, m, T) - F1All(p0 + q0, ps, UVScale, m, T)
            +
            F2All(p0 - q0, ps, k, m, T) - F2All(p0 - q0, ps, UVScale, m, T) +
            F2All(p0 + q0, ps, k, m, T) - F2All(p0 + q0, ps, UVScale, m, T)
        )
        # +
        # (2 + Npi) * (
        #     F1All(1e-8 - 1e-14, 1e-8, k, m, T) -
        #     F1All(1e-8 - 1e-14, 1e-8, UVScale, m, T) +
        #     F2All(1e-8 - 1e-14, 1e-8, k, m, T) -
        #     F2All(1e-8 - 1e-14, 1e-8, UVScale, m, T)
        # )
    )
end

@doc raw"""
    VImintqsSimple(p0, ps, k, T, Npi, m, lamda,UVScale)

compute $\int_0^kq_s^2dqs\int_{-1}^{1}d\cos\theta\mathrm{Im}V_k(q_0)$, the `k` dependence of $\lambda_{4\pi}$ and $m$ are neglected.

`VImintqsSimple` only contains $V(q_0)$, no $V(-q_0)$

`VImintqsSimple` contains type-1 delta function

`VImintqsSimple` doesn't contains type-2 delta function
"""
function VImintqsSimple(p0, ps, k, T, Npi, m, lamda,UVScale)
    hcubature(
        x ->
            x[1]^2 * VImSimple(
                p0,
                sqrt(x[1]^2 + ps^2 + 2 * x[1] * x[2] * ps),
                Epi(k, m),
                k,
                m,
                T,
                Npi,
                lamda,UVScale,
            ),
        [0.0, -1.0],
        [k, 1.0],atol=1e-6,rtol=1e-6,initdiv=10,maxevals=8000,
    )[1]
end




function VImintqsSimple(p0, ps,q0, k, T, Npi, m, lamda,UVScale)
    hcubature(
        x ->
            x[1]^2 * VImSimple(
                p0,
                sqrt(x[1]^2 + ps^2 + 2 * x[1] * x[2] * ps),
                q0,
                k,
                m,
                T,
                Npi,
                lamda,UVScale,
            ),
        [0.0, -1.0],
        [k, 1.0],atol=1e-2,rtol=1e-2,initdiv=40,maxevals=16000,
    )[1]
end


@doc raw"""
    propImSimple(p0, ps, T,IRScale,UVScale, Npi, m, lamda)
`x[1]` is `qs`, `x[2]` is `costh`, `x[3]` is `k`

# Arguments
- `m`: mass square, it's a constant number.
- `lamda`: $\lambda_{4\pi}$, it's a constant number.
"""
function propImSimple(p0, ps, T, IRScale, UVScale, Npi, m, lamda)
    -hquadrature(
        x ->
            2*Coeffgamm2Simple(x, T, m) *
            VImintqsSimple(p0, ps, x, T, Npi, m, lamda,UVScale),
        IRScale,
        UVScale,
        atol = 1e-4,
        rtol = 1e-4,maxevals=200,
    )[1]
end
