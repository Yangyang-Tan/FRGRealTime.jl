using FRGRealTime.MathFuns

function massflow1(k, m2, T, lambda, Nf=4)
    (1 / 6) * Nf * lambda * TF.Fb1(k, m2, T)
end


function massflow2(k, m2, T, lambda, Nf=4)
    (1 / 6) * (Nf + 2) * lambda * TF.Fb1(k, m2, T)
end

function lambdaflow1(k, m2, T, lambda)
    -(1 / 3) * lambda^2 * TF.Fb2(k, m2, T)
end


function lambdaflow2(k, m2, T, lambda, Nf=4)
    -(1 / 6) * lambda^2 * (Nf + 8) * TF.Fb2(k, m2, T)
end



function Imlambdaflow1(p0, k, m2, lambda, T, Imlambda, Relambda, ϵ=0.1)
    -(1 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0, k, m2, T, ϵ) +
        2 * TF.ReFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        TF.ImFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Relambdaflow1(p0, k, m2, lambda, T, Imlambda, Relambda, ϵ=0.1)
    -(1 / 3) * (
        -Imlambda^2 * TF.ReFb2(p0, k, m2, T, ϵ) -
        2 * TF.ImFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        TF.ReFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Imlambdaflow2(p0, k, m2, lambda, T, Imlambda, Relambda, ϵ=0.1)
    -(2 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0, k, m2, T, ϵ) +
        2 * TF.ReFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        TF.ImFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Relambdaflow2(
    p0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    Relambda,
    ϵ=0.1,
    Nf=4,
)
    -(2 / 3) * (
        -Imlambda^2 * TF.ReFb2(p0, k, m2, T, ϵ) -
        2 * TF.ImFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        TF.ReFb2(p0, k, m2, T, ϵ) * Relambda^2
    ) - (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end



function Imlambdaflow3(p0, q0, k, m2, lambda, T, Imlambda, Relambda, ϵ=0.1)
    -(1 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0 + q0, k, m2, T, ϵ) +
        2 *
        TF.ReFb2(p0 + q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda +
        TF.ImFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2
    ) -
    (1 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0 - q0, k, m2, T, ϵ) +
        2 *
        TF.ReFb2(p0 - q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda +
        TF.ImFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2
    )
end

function Imlambdaflow3_ac(p0, q0, k, m2, lambda, T, Imlambda, Relambda, ϵ=0.1)
    -(1 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0 + q0, k, m2, T, ϵ) +
        2 *
        TF.ReFb2_ac(p0 + q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda +
        TF.ImFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2
    ) -
    (1 / 3) * (
        -Imlambda^2 * TF.ImFb2(p0 - q0, k, m2, T, ϵ) +
        2 *
        TF.ReFb2_ac(p0 - q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda +
        TF.ImFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2
    )
end
function Relambdaflow3_ac(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    Relambda,
    ϵ=0.1,
    Nf=4,
)
    -(1 / 3) * (
        -Imlambda^2 * TF.ReFb2_ac(p0 + q0, k, m2, T, ϵ) -
        2 *
        TF.ImFb2(p0 + q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda + TF.ReFb2_ac(p0 + q0, k, m2, T, ϵ) * Relambda^2
    ) -
    (1 / 3) * (
        -Imlambda^2 * TF.ReFb2_ac(p0 - q0, k, m2, T, ϵ) -
        2 *
        TF.ImFb2(p0 - q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda + TF.ReFb2_ac(p0 - q0, k, m2, T, ϵ) * Relambda^2
    ) - (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end


function Relambdaflow3(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    Relambda,
    ϵ=0.1,
    Nf=4,
)
    -(1 / 3) * (
        -Imlambda^2 * TF.ReFb2(p0 + q0, k, m2, T, ϵ) -
        2 *
        TF.ImFb2(p0 + q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda + TF.ReFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2
    ) -
    (1 / 3) * (
        -Imlambda^2 * TF.ReFb2(p0 - q0, k, m2, T, ϵ) -
        2 *
        TF.ImFb2(p0 - q0, k, m2, T, ϵ) *
        Relambda *
        Imlambda + TF.ReFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2
    ) - (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end

function ReFb2(p0, k, m2, T, ϵ = 0.1)
    if p0>=0.0
        ReFb2positive(p0, k, m2, T, ϵ)
    elseif p0<0.0
        ReFb2positive(-p0, k, m2, T, ϵ)
    end
end




function Imlambdaflow4(p0, q0, k, m2, lambda, T, Imlambda, Relambda, ϵ = 0.1)
    -(1 / 3) * (TF.ImFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 3) * (TF.ImFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2)
end

function Imlambdaflow4_ac(p0, q0, k, m2, lambda, T, Imlambda, Relambda, ϵ = 0.1)
    -(1 / 3) * (TF.ImFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 3) * (TF.ImFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2)
end


function Relambdaflow4(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    Relambda,
    ϵ = 0.1,
    Nf = 4,
)
    -(1 / 3) * (TF.ReFb2(p0 + q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 3) * (TF.ReFb2(p0 - q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end



function Relambdaflow4_ac(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    Relambda,
    ϵ = 0.1,
    Nf = 4,
)
    -(1 / 3) * (TF.ReFb2_ac(p0 + q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 3) * (TF.ReFb2_ac(p0 - q0, k, m2, T, ϵ) * Relambda^2) -
    (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end


# function twopointflow(p0, k, m2, T, lambda)
#     -1 / 48 * (
#         k^4 *
#         csch(sqrt(k^2 + m2) / (2 * T))^2 *
#         (
#             lambda(k, abs(sqrt(k^2 + m2) - p0)) +
#             lambda(k, abs(sqrt(k^2 + m2) + p0))
#         ) *
#         (sqrt(k^2 + m2) + T * sinh(sqrt(k^2 + m2) / T))
#     ) / ((k^2 + m2)^(3 / 2) * pi^2 * T)
# end


function twopointflow2(p0, k, m2, T, lambda,Nf=4)
    1 / 144 *(Nf)* (
        k^4 *
        csch(sqrt(k^2 + m2) / (2 * T))^2 *
        lambda(k) *
        (sqrt(k^2 + m2) + T * sinh(sqrt(k^2 + m2) / T))
    ) / ((k^2 + m2)^(3 / 2) * pi^2 * T)
end


function twopointflow3(p0, k, m2, T, lambda,Nf=4)
    1 / 144 *(Nf)* (
        k^4 *
        csch(sqrt(k^2 + m2) / (2 * T))^2 *
        lambda *
        (sqrt(k^2 + m2) + T * sinh(sqrt(k^2 + m2) / T))
    ) / ((k^2 + m2)^(3 / 2) * pi^2 * T)
end



function RelambdaflowMktest(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    Imlambda,
    RelambdaMk,
    ϵ = 0.1,
    Nf = 4,
)
    -(1 / 3) * (ReFb2Mkpositivetest(p0, k, m2, T, ϵ) * RelambdaMk^2)
end

function ReFb2Mkpositivetest(p0, k, m2, T, ϵ = 0.1)
    -1 / 48 * (k^4 * (ϵ - 4 * p0 * Eb(k, m2) + 8 * Eb(k, m2)^2)) / (
        pi^2 *
        Eb(k, m2)^3 *
        (p0 + 2 * Eb(k, m2)) *
        (p0^2 + ϵ + 4 * Eb(k, m2) * (-p0 + Eb(k, m2)))
    )
end
