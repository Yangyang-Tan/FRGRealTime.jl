function massflow1(k,m2,T,lambda,Nf=4)
    (1/6)*Nf*lambda*ThresholdFunctions.Fb1(k, m2, T)
end


function massflow2(k,m2,T,lambda,Nf=4)
    (1/6)*(Nf+2)*lambda*ThresholdFunctions.Fb1(k, m2, T)
end

function lambdaflow1(k,m2,T,lambda)
    -(1/3) *lambda^2*ThresholdFunctions.Fb2(k, m2, T)
end


function lambdaflow2(k,m2,T,lambda,Nf=4)
    -(1/6) *lambda^2*(Nf+8)*ThresholdFunctions.Fb2(k, m2, T)
end



function Imlambdaflow1(p0, k, m2,lambda, T, Imlambda,Relambda, ϵ = 0.1)
    -(1 / 3) * (
        -Imlambda^2 * ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) +
        2 * ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Relambdaflow1(p0, k, m2, lambda,T, Imlambda,Relambda, ϵ = 0.1)
    -(1 / 3) * (
        -Imlambda^2 * ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) -
        2 * ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Imlambdaflow2(p0, k, m2,lambda, T, Imlambda,Relambda,ϵ = 0.1)
    -(2 / 3) * (
        -Imlambda^2 * ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) +
        2 * ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) * Relambda^2
    )
end

function Relambdaflow2(p0, k, m2,lambda, T, Imlambda,Relambda,ϵ = 0.1,Nf=4)
    -(2 / 3) * (
        -Imlambda^2 * ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) -
        2 * ThresholdFunctions.ImFb2(p0, k, m2, T, ϵ) * Relambda * Imlambda +
        ThresholdFunctions.ReFb2(p0, k, m2, T, ϵ) * Relambda^2
    )-(1/6) *lambda^2*(Nf+4)*ThresholdFunctions.Fb2(k, m2, T)
end

function twopointflow(p0, k, m2, T,lambda)
    -1 / 48 * (k^4 * csch(sqrt(k^2 + m2) / (2 * T))^2 * (lambda(k, abs(sqrt(k^2 + m2) - p0)) + lambda(k, abs(sqrt(k^2 + m2) + p0))) * (sqrt(k^2 + m2) + T * sinh(sqrt(k^2 + m2) / T))) / ((k^2 + m2)^(3 / 2) * pi^2 * T)
end


function twopointflow2(p0, k, m2, T, lambda)
    -1 / 24 * (
        k^4 *
        csch(sqrt(k^2 + m2) / (2 * T))^2 *
        lambda(k) *
        (sqrt(k^2 + m2) + T * sinh(sqrt(k^2 + m2) / T))
    ) / ((k^2 + m2)^(3 / 2) * pi^2 * T)
end
