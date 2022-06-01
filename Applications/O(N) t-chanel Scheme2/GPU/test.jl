using BenchmarkTools
Imouttest = CuArray(zeros(2041, 2041))
Reouttest = CuArray(sin.(p0grid) .^ 2 .+ q0grid' .^ 2)
Imouttest2 = CuArray(randn(2041, 2041))
Reouttest2 = CuArray(randn(2041, 2041))
p0grid = CuArray(-102:0.1:102)
q0grid = CuArray(-102:0.1:102)
testLambdaGrid = LambdaGridini(Float64, 2041)
reverse2!(
    testLambdaGrid,
    Imouttest,
    Reouttest,
    lengthy=2041,
    dq0=0.1,
    q0min=102.0,
    Ek=10.02,
    threads=1024,
)

@benchmark begin
    ImRelambdaflow_Launch(
        $Imouttest,
        $Reouttest,
        $p0grid,
        $q0grid,
        $testLambdaGrid,
        500.0,
        10.0,
        145.0,
        1.0,
        4.0,
        threads=512,
    )
    synchronize()
end

@benchmark begin
    ImRelambdaflow_Launchxy(
        $Imouttest,
        $Reouttest,
        $p0grid,
        $q0grid,
        $testLambdaGrid,
        500.0,
        10.0,
        145.0,
        1.0,
        4.0,
        threads=(16, 16),
    )
    synchronize()
end



ImRelambdaflow_Launch(
    Imouttest,
    Reouttest,
    p0grid,
    q0grid,
    testLambdaGrid,
    500.0,
    10.0,
    145.0,
    1.0,
    4.0,
    threads=512,
)
ImRelambdaflow_Launchxy(
    Imouttest2,
    Reouttest2,
    p0grid,
    q0grid,
    testLambdaGrid,
    500.0,
    10.0,
    145.0,
    1.0,
    4.0,
    threads=(32, 16),
)
Imouttest2 == Imouttest
Reouttest2 .- Reouttest


Imouttest


Imouttest2
Reouttest

Reouttest[7, end-1]
Reouttest[7, 2]



Reouttest2[7, end-1]
Reouttest2[7, 2]


Reouttest


Relambdaflow_GPU(
    200.0,
    200.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    10.0,
    10.0,
    10.0,
    2.0,
    1.0,
    100.0,
    10.0,
    145.0,
    10.0,
    4.0,
)


Relambdaflow_GPU(
    -200.0,
    200.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    10.0,
    10.0,
    10.0,
    2.0,
    1.0,
    100.0,
    10.0,
    145.0,
    10.0,
    4.0,
)

testLambdaGrid.Imlambdap0q0


using CUDA

function quadspline(x, n, c::AbstractArray, t::AbstractArray)
    (
        ((x - t[1+n]) * (x * (c[1] - c[2]) + c[2] * t[-1+n] - c[1] * t[1+n])) /
        (t[-1+n] - t[1+n]) +
        ((x - t[n]) * (x * (-c[2] + c[3]) - c[3] * t[n] + c[2] * t[2+n])) /
        (t[n] - t[2+n])
    ) / (t[n] - t[1+n])
end

testa = CuArray(randn(10))


testn = BSplineKit.Splines.knot_interval(knots(itp), 3.0)

testn |> println

testa = randn(1000, 1000)


@benchmark BSplineKit.Splines.spline_kernel.(
    Ref(coefficients(itp)),
    Ref(knots(itp)),
    3,
    $testa,
    Ref(BSplineOrder(3)),
)

@benchmark BSplineKit.Splines.spline_kernel_alt.(
    Ref(coefficients(itp)),
    Ref(knots(itp)),
    3,
    $testa,
    Ref(BSplineOrder(3)),
)

@benchmark quadspline.($testa, 3, Ref(coefficients(itp)), Ref(knots(itp)))

using BSplineKit
xs = 1.0:1.0:3.0
ys = sin.(xs)
itp = interpolate(xs, ys, BSplineOrder(3))
1


import Base.+, Base.*, Base./

+(f, g) = x -> f(x) + g(x)
*(t::Number, g) = x -> t * g(x)

myinterpolate(a, b) = t -> (1 - t) * a + t * b
b1(p1, p2) = myinterpolate(p1, p2)
b2(p1, p2, p3) = myinterpolate(b1(p1, p2), b1(p2, p3))
b3(p1, p2, p3, p4) = myinterpolate(b2(p1, p2, p3), b2(p2, p3, p4))


b2(coefficients(itp)...)((2.3 - 1) / 2)((2.3 - 1) / 2)

quadspline.(1.1:0.1:2.9, 3, Ref(coefficients(itp)), Ref(knots(itp))) ≈ itp.(1.1:0.1:2.9)






using Base.Cartesian

b1(coefficients(itp)...)(0.2)

ys[1] / 2 + ys[2] / 2

itp(1.2)





BSplineKit.Splines.knot_interval(knots(itp), 2.0)
coefficients(itp)
BSplineOrder(3)
using BSplineKit
xdata = 1:1.0:3  # points don't need to be uniformly distributed
ydata = rand(length(xdata))
itp = interpolate(xdata, ydata, BSplineOrder(3))
vgpu = CuArray(1:0.1:2)
itp(vgpu)
B = BSplineBasis(BSplineOrder(3), xdata)
B(10)

BasisFunction(B, 2)
