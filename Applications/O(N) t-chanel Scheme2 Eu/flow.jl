using FRGRealTime.MathFuns

function massflow2(k, m2, T, lambda, Nf=4)
    (1 / 6) * (Nf + 2) * lambda * TF.Fb1(k, m2, T)
end


#zero momentum lambda flow
function lambdaflow0(k, m2, T, lambda, Nf=4)
    -(1 / 6) * lambda^2 * (Nf + 8) * TF.Fb2(k, m2, T)
end


function RelambdaflowEu(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    RelambdaEu,
    Nf = 4,
)
    -(1 / 3) * (TF.ReFb2Eu(p0, k, m2, T) * RelambdaEu^2) -
    (1 / 3) * (TF.ReFb2Eu(p0, k, m2, T) * RelambdaEu^2) -
    (1 / 6) * lambda^2 * (Nf + 4) * TF.Fb2(k, m2, T)
end

function RelambdaflowEutest(
    p0,
    q0,
    k,
    m2,
    lambda,
    T,
    RelambdaEu,
    Nf = 4,
)
    -(1 / 3) * (ReFb2Eutest(p0, k, m2, T) * RelambdaEu^2)
end




function ReFb2Eutest(p0, k, m2, T)
    -1 / 12 * k^4 / (pi^2 * Eb(k, m2)^2 * (p0^2 + 4 * Eb(k, m2)^2))
end
