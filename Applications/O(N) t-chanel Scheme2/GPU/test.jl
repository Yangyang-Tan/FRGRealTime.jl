using BenchmarkTools
Imouttest=CuArray(zeros(2041,2041))
Reouttest=CuArray(sin.(p0grid).^2 .+q0grid'.^2)
Imouttest2=CuArray(randn(2041,2041))
Reouttest2=CuArray(randn(2041,2041))
p0grid=CuArray(-102:0.1:102)
q0grid=CuArray(-102:0.1:102)
testLambdaGrid=LambdaGridini(Float64, 2041)
reverse2!(
            testLambdaGrid,
            Imouttest,
            Reouttest,
            lengthy = 2041,
            dq0 = 0.1,
            q0min = 102.0,
            Ek = 10.02,
            threads = 1024,
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
    threads = 512,
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
    threads = (16,16),
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
    threads = 512,
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
    threads = (32,16),
    )
Imouttest2 == Imouttest
Reouttest2 .- Reouttest


Imouttest


Imouttest2
Reouttest

Reouttest[7,end-1]
Reouttest[7,2]



Reouttest2[7,end-1]
Reouttest2[7,2]


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

testLambdaGrid.Relambdaq0Ek[1:1020]==reverse(testLambdaGrid.Relambdaq0Ek[1022:end])
expand_index
using CUDA

BSplineKit.Splines.spline_kernel(coefficients(itp),knots(itp),3,2.0,BSplineOrder(3))

using BSplineKit
xs=1.0:1.0:2.0
ys=sin.(xs)
itp=interpolate(xs, ys, BSplineOrder(2))



import Base.+, Base.*, Base./

+(f, g) = x -> f(x) + g(x)
*(t::Number, g) = x -> t * g(x)

myinterpolate(a, b) = t ->(1 - t) * a + t * b
b1(p1, p2) = myinterpolate(p1, p2)
b2(p1, p2, p3) = myinterpolate(b1(p1, p2), b1(p2, p3))
b3(p1, p2, p3, p4) = myinterpolate(b2(p1, p2, p3), b2(p2, p3, p4))

T((1 - α) * $d_p + α * $d_j)

@generated function foo(x)
           Core.println(x)
           return :(x * x)
       end


foo(1.0)




b1(coefficients(itp)...)(1.5)

ys[1]/2+ys[2]/2

itp(1.5)





BSplineKit.Splines.knot_interval(knots(itp),2.0)
coefficients(itp)
BSplineOrder(3)
using BSplineKit
xdata = 1:1.0:3  # points don't need to be uniformly distributed
ydata = rand(length(xdata))
itp=interpolate(xdata, ydata, BSplineOrder(3))
vgpu=CuArray(1:0.1:2)
itp(vgpu)
B=BSplineBasis(BSplineOrder(3),xdata)
B(10)

BasisFunction(B,2)
