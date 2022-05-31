using Interpolations
tsolparallel1

range(1.0, stop = 400.0, length = 2000)

lambda_k = 1:0.2:1000


lambda_omega = range(1.0, stop = 2000.0, length = 4000)


f(x1, x2) = log(x1 + x2)

function getImlambdaData(k, omega)
    n = round(Int64, (omega - 1) / step(lambda_omega)) + 1
    tsolparallel1[n](k)[1]
end
function getRelambdaData(k, omega)
    n = round(Int64, (omega - 1) / step(lambda_omega)) + 1
    tsolparallel1[n](k)[2]
end


ImA = [getImlambdaData(x1, x2) for x1 in lambda_k, x2 in lambda_omega]
ReA = [getRelambdaData(x1, x2) for x1 in lambda_k, x2 in lambda_omega]


Imitp = interpolate(ImA, BSpline(Linear()))
Imsitp = scale(Imitp, lambda_k, lambda_omega)
Imetpf = extrapolate(Imsitp, Flat())

Reitp = interpolate(ReA, BSpline(Linear()))
Resitp = scale(Reitp, lambda_k, lambda_omega)
Reetpf = extrapolate(Resitp, Flat())




plot(k -> -getImlambdaData(k, lambda_omega[800]), 1, 800)

plot!(k -> -Imetpf(k, lambda_omega[1000] + 100), 1, 800)

m2fun(x) = sol1(x)[1]

twopointflow(100.0, 10.0, m2fun(10.0), 145.0, Imetpf)

using HCubature

function Imgamma2(p0)
    -hcubature(
        x -> twopointflow(p0, x[1], m2fun(x[1]), 145.0, Imetpf),
        [1.0],
        [800.0],
        atol = 1e-4,
        rtol = 1e-4,
    )[1]
end

function Regamma2(p0)
    -hcubature(
        x -> twopointflow(p0, x[1], m2fun(x[1]), 145.0, Reetpf),
        [1.0],
        [800.0],
        atol = 1e-4,
        rtol = 1e-4,
    )[1] - p0^2 + m2fun(800.0)
end

function bosonicspec(p0)
    -2 * Imgamma2(p0) / (Imgamma2(p0)^2 + Regamma2(p0)^2)
end

plot(p0 -> Imgamma2(p0), 1.0, 400.0)
plot(p0 -> Regamma2(p0), 1.0, 400.0)

plot(p0 -> -bosonicspec(p0), 1.0, 400.0, yaxis = :log)
