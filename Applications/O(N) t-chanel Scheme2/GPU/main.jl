path = @__DIR__
include(joinpath(path, "iniGPU.jl"))
include(joinpath(path, "GPUsolver.jl"))
config_spec.atol = 1e-8
config_spec.rtol = 1e-8

@time solselfGPU_log = tchanelSolveFourPointGPU(
    # 9.4685f0,
    145.835,
    kmin = 20.0,
    ϵ = 5.0,
    pgridmax = 100.0,
    dp0 = 0.5,
    dq0 = 0.5,
    Nf = 4.0,
    config = config_spec,
)

plot(solselfGPU_log.sol.u[end][:, 102-20, 1] + solselfGPU_log.sol.u[end][:, 102+20, 1])
1

@time solselfGPU_AB = tchanelSolveFourPointGPU_AB(
    # 9.4685f0,
    145.835,
    kmin = 20.0,
    ϵ = 2.0,
    pgridmax = 100.0,
    dp0 = 0.125 / 0.125,
    dq0 = 0.125 / 0.125,
    Nf = 4.0,
    config = config_spec,
)

plot(solselfGPU_log.sol.u[end][:,102-20,1]+solselfGPU_log.sol.u[end][:,102+20,1])

plot(solselfGPU_log.sol.u[end][:,102-20,2]+solselfGPU_log.sol.u[end][:,102+20,2])


plot(solselfGPU_AB.sol.u[end][:,102-20,1]+solselfGPU_AB.sol.u[end][:,102+20,1])
plot!(solselfGPU_AB.sol.u[end][:,102-20,2]+solselfGPU_AB.sol.u[end][:,102+20,2]+solselfGPU_AB.sol.u[end][:,102-20,3]+solselfGPU_AB.sol.u[end][:,102+20,3])

plot(solselfGPU_AB.sol.u[end][:,102-20,2]+solselfGPU_AB.sol.u[end][:,102+20,2])
plot(abs.(solselfGPU_AB.sol.u[end][:,102-20,3]+solselfGPU_AB.sol.u[end][:,102+20,3]))


plot(solselfGPU_log.p0,solselfGPU_log.spec,yaxis=:log)
plot!(solselfGPU_AB.p0,solselfGPU_AB.spec,yaxis=:log)

solselfGPU_AB.spec




plot(vcat([-reverse(solselfGPU_log.p0),solselfGPU_log.p0[2:end]]...),Array(solselfGPU_log.sol.u[end][:,11,1]+solselfGPU_log.sol.u[end][:,end-11,2]))



solselfGPU_log.m2[end]

plot(vcat([-reverse(solselfGPU_log.p0),solselfGPU_log.p0[2:end]]...),Array(solselfGPU_log.sol.u[end][:,1,1]+solselfGPU_log.sol.u[end][:,end-1,1]))

solselfGPU_log64 = tchanelSolveFourPointGPU_log(
    # 9.4685f0,
    135.2155,
    kmin = 1.0,
    ϵ = 1.0,
    pgridmax = 80.0,
    dp0 = 0.125 / 2,
    dq0 = 0.125 / 2,
    Nf = 4.0,
    config = config_spec,
)


solselfGPU_log264 = tchanelSolveFourPointGPU_log2(
    # 9.4685f0,
    # 128.232,
    128.532,
    kmin = 1.0,
    ϵ = 1.0,
    pgridmax = 20.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 4.0,
    config = config_spec,
)


plot(vcat([-reverse(solselfGPU_log264.p0),solselfGPU_log264.p0[2:end]]...),Array(solselfGPU_log264.sol.u[end][:,1,1]+solselfGPU_log264.sol.u[end][:,end-1,1]))
solselfGPU_log264.m2[end:-1:end-10]
plot(solselfGPU_log264.k[end:-1:end-10],solselfGPU_log264.m2[end:-1:end-10])

plot(vcat([-reverse(solselfGPU_log264.p0),solselfGPU_log264.p0[2:end]]...),Array(solselfGPU_log264.sol.u[end][:,1,2]+solselfGPU_log264.sol.u[end][:,end-1,2]))


solselfGPU_log264.p0

solselfGPU_log264_N1 = tchanelSolveFourPointGPU_log2(
    # 9.4685f0,
    107.55092,
    m2ini = Float64(-0.18216546875*0.5 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 1.0,
    pgridmax = 10.0,
    pgridmin=2*1e-3,
    gridlength=1400,
    Nf = 1.0,
    config = config_spec,
)


solselfGPU_log264_N120 = tchanelSolveFourPointGPU_log2(
    # 9.4685f0,
    138.85925,
    m2ini = Float64(-0.18216546875*2.6 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 1.0,
    pgridmax = 10.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 20.0,
    config = config_spec,
)

solselfGPU_log232 = tchanelSolveFourPointGPU_log2(
    # 9.4685f0,
    138.8598f0,
    m2ini = Float32(-0.18216546875*2.6 * lpakrang.UV^2),
    kmin = 1.0f0,
    ϵ = 1.0f0,
    pgridmax = 10.0f0,
    pgridmin=1f-3,
    gridlength=1000,
    Nf = 20.0f0,
    config = config_spec,
)
solselfGPU_log264.m2[end]
solselfGPU_log264_N1.m2[end]
Array(solselfGPU_log264.sol.u[end][:,:,1])




heatmap(Array(solselfGPU_log264.sol.u[end][:,:,1]))



solselfGPU_log232.m2[end]
solselfGPU_log264.spec
plot(solselfGPU_log264.p0[2:end-10], solselfGPU_log264.spec[2:end-10],yaxis=:log,xaxis=:log)
plot!(solselfGPU_log264_N1.p0[2:end], solselfGPU_log264_N1.spec[2:end],yaxis=:log,xaxis=:log)

plot!(solselfGPU_log264_N120.p0[2:end], solselfGPU_log264_N120.spec[2:end],yaxis=:log,xaxis=:log)


solselfGPU.p0


plot(solselfGPU.p0, solselfGPU.spec,yaxis=:log)

plot!(solselfGPU_log64.p0, solselfGPU_log64.spec,yaxis=:log)

plot!(solselfGPU_log264.p0, solselfGPU_log264.spec,yaxis=:log)


plot(solselfGPU_log64.k, solselfGPU_log64.λ0)
plot!(solselfGPU_log264.k, solselfGPU_log264.λ0)

sp1=Spline1D(log.(solselfGPU_log264_N1.p0[2:end]), log.(solselfGPU_log264_N1.spec[2:end]))

sp3=Spline1D(log.(solselfGPU_log264_N120.p0[2:end]), log.(solselfGPU_log264_N120.spec[2:end]))

xdata=log.(solselfGPU_log264_N120.p0[2:end])[290:end-630]
ydata=log.(solselfGPU_log264_N120.spec[2:end])[290:end-630]

xdata=log.(solselfGPU_log264_N1.p0[2:end])[290:end-630]
ydata=log.(solselfGPU_log264_N1.spec[2:end])[290:end-630]

fit3 = curve_fit(Polynomial, xdata, ydata, 1)
fit1 = curve_fit(Polynomial, xdata, ydata, 1)

1
function f(x)
    I=0
    for i in 1:5
        I=I+sin(i*x)
    end
    return I
end

f(2.0)

plot(fit3,-5.2,-1.5)
plot(fit1,-6.2,2.5)
plot!(x->sp1(x),-6.2,2.5)
plot!(x->sp3(x),-5.2,-1.5)

plot!(x->sp3(x),-4.2,-3.5)

plot!(x->derivative(sp3,x),-4.2,-3.5)


plot(1:2000,exp.(range(log(1e-3),log(400),2000)))




1.5/step(range(1.0,2.5,3))

2.5/(range(1.0,2.5,17)[2])

collect(range(1.0,2.5,17))


using LinearAlgebra
A=randn(Float32,4000,4000)
B=randn(Float32,4000,4000)
C=randn(Float32,4000,4000)
A=CuArray(randn(Float32,4000,4000))
B=CuArray(randn(Float32,4000,4000))
C=CuArray(randn(Float32,4000,4000))
@benchmark CUDA.@sync mul!(C,A,B)

@benchmark @sync mul!(C,A,B)



145.0

2.0

plot(solselfGPU_log.k, solselfGPU_log.λ0)
plot!(solselfGPU.k, solselfGPU.λ0)

plot(solselfGPU_log.k, solselfGPU.m2)



solselfGPU1.m2[end]
solselfGPU.m2[end]

plot!(solselfGPU.p0, solselfGPU.ReΓ2)

plot(solselfGPU.p0, solselfGPU.ImΓ2)
plot(solselfGPU.p0, solselfGPU.ImΓ2ini)


plot(solselfGPU.p0[2:end], solselfGPU.spec[2:end], xaxis = :log, yaxis = :log)
plot!(solselfGPU1.p0[2:end], solselfGPU1.spec[2:end], xaxis = :log, yaxis = :log)



fit1 = Spline1D(log.(solselfGPU.p0[2:end]), log.(solselfGPU.spec[2:end]))
log.(solselfGPU.p0[2:end])

plot(x -> derivative(fit1, x), -0.0, 2.0)

function reverse!(A::CuArray, B::CuArray; lengthy::Int)
    @tullio A[x, y] = B[x, lengthy-y]
end



testv1 = CUDA.randn(Float64, 2000, 2000)
testv1s = CUDA.randn(Float64, 2000, 2000)

@benchmark testRe .= Relambdaflow_GPU(testp0, testq0, test1, 10.0, 0.0, 145.0, 1.0, 4.0)

@benchmark testfun1()
# @benchmark

@benchmark CUDA.@sync ImRelambdaflow_Launch(testIm, testRe, testp0, testq0, test1, 10.0, 0.0, 145.0, 1.0, 4.0,threads=100)


@benchmark CUDA.@sync ImRelambdaflow_Launch(testIm, testRe, testp0, testq0, test1, 10.0f0, 0.0f0, 145.0f0, 1.0f0, 4.0f0,threads=100)


# test1
t1=1.0
testp0=CuArray(randn(typeof(t1),4000))
testq0=CuArray(randn(typeof(t1),4000)')

# @benchmark myreverse!(testv1,testv1s,lengthy=2000)

# @benchmark myreverse2!(testv1,testv1s,lengthy=2000)
test1=LambdaGridini(typeof(t1),4000)




testIm=CuArray(fill(0.001*t1,4000,4000))
testRe=CuArray(fill(24.0*t1,4000,4000))
# @benchmark reverse2!(test1,testIm,testRe,lengthy=100,dq0=1.0f0,q0min=-10.0f0,Ek=5.0f0,)
# @benchmark reverse2!(test1,testIm,testRe,lengthy=4000,dq0=1.0,q0min=-10.0,Ek=5.0,)

@benchmark CUDA.@sync reverse2!(test1,testIm,testRe,lengthy=4000,dq0=1.0,q0min=-10.0,Ek=5.2,)

@benchmark CUDA.@sync reverse2!(test1,testIm,testRe,lengthy=4000,dq0=1.0f0,q0min=-50.0f0,Ek=5.1f0,)



reverse2!(test1,testIm,testRe,lengthy=4000,dq0=1.0,q0min=-10.0,Ek=5.2,)

test1.Imlambdap0Ek

# @benchmark testRe .=testfun.(testIm)


# @benchmark testRe .=testfun.(testIm)


solselfGPU_log264 = tchanelSolveFourPointGPU_log(
    # 9.4685f0,
    128.232,
    kmin = 1.0,
    ϵ = 3.0,
    pgridmax = 20.0,
    pgridmin=1e-3,
    gridlength=999,
    Nf = 4.0,
    config = config_spec,
)


solselfGPU_log264 = tchanelSolveFourPointGPU_log(
    # 9.4685f0,
    # 128.22675,
    107.54732005,
    m2ini = Float64(-0.18216546875*0.5 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 4.0,
    pgridmax = 20.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 1.0,
    config = config_spec,
)

solselfGPU_log264_20 = tchanelSolveFourPointGPU_log(
    # 9.4685f0,
    # 128.22675,
    138.866062,
    m2ini = Float64(-0.18216546875*2.6 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 4.0,
    pgridmax = 20.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 20.0,
    config = config_spec,
)



solselfGPU_log264_100 = tchanelSolveFourPointGPU_log(
    # 9.4685f0,
    # 128.22675,
    159.0176405,
    m2ini = Float64(-0.18216546875*4.24 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 4.0,
    pgridmax = 20.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 100.0,
    config = config_spec,
)

solselfGPU_log264_100.m2[end]


solselfGPU_log264_20.m2[end]


solselfGPU_log264.m2[end]
plot(solselfGPU_log264.p0[2:end],solselfGPU_log264.spec[2:end],xaxis=:log,yaxis=:log)

plot!(solselfGPU_log264_20.p0[2:end],solselfGPU_log264_20.spec[2:end],xaxis=:log,yaxis=:log)

fitfun=Spline1D(solselfGPU_log264.p0[2:end].|>log,solselfGPU_log264.spec[2:end].|>log)
plot(x->derivative(fitfun,x),-6.9,0.0)



plot!(solselfGPU_log264.p0[2:end].|>log,solselfGPU_log264.spec[2:end].|>log,label="N=1,slop=0.9968")
plot!(solselfGPU_log264_20.p0[2:end].|>log,solselfGPU_log264_20.spec[2:end].|>log,label="N=20,slop=0.9990")
plot!(solselfGPU_log264_100.p0[2:end].|>log,solselfGPU_log264_100.spec[2:end].|>log,label="N=100,slop=0.9984")

xdata1=(solselfGPU_log264.p0[2:end].|>log)[10:200]
ydata1=(solselfGPU_log264.spec[2:end].|>log)[10:200]
xdata2=(solselfGPU_log264_20.p0[2:end].|>log)[10:200]
ydata2=(solselfGPU_log264_20.spec[2:end].|>log)[10:200]
xdata3=(solselfGPU_log264_100.p0[2:end].|>log)[10:200]
ydata3=(solselfGPU_log264_100.spec[2:end].|>log)[10:200]

fit1 = curve_fit(Polynomial, xdata1, ydata1, 1)
fit2 = curve_fit(Polynomial, xdata2, ydata2, 1)
fit2 = curve_fit(Polynomial, xdata3, ydata3, 1)


tchanelSolveFourPointGPU(
    # 9.4685f0,
    # 128.22675,
    107.54732005,
    m2ini = Float64(-0.18216546875*0.5 * lpakrang.UV^2),
    kmin = 1.0,
    ϵ = 4.0,
    pgridmax = 100.0,
    pgridmin=1e-3,
    gridlength=1000,
    Nf = 1.0,
    config = config_spec,
)
