using Distributed
using SharedArrays
addprocs(8)
nprocs()
@everywhere using FRGRealTime, DelimitedFiles, Dierckx, Plots
@everywhere kdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/CTPcode/Zero-ODEVersion/realtime_zero_k.dat")[:,1]|>reverse
@everywhere lamdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/CTPcode/Zero-ODEVersion/realtime_zero_lamdak.dat")[:,1]|>reverse
@everywhere m2data=readdlm("/home/tyy/Documents/CTP-fRG-Test/CTPcode/Zero-ODEVersion/realtime_zero_m2k.dat")[:,1]|>reverse
@everywhere lamfun=Spline1D(kdata,-0.5*lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)
v1=SharedArray(collect(400.1:6.5:500.0))
@everywhere FRGRealTime.propReintqs(4.0, 10.0, 145.0,1.0,400.0, 4.0, m2fun, lamfun)


plot(x->m2fun(x),1.0,800.0)
plot(x->lamfun(x),1.0,800.0)

FRGRealTime.propReintqs(
    500.0000001,
    0.000005,
    35.05,
    780.0,
    800.0,
    4.0,
    m2fun,
    lamfun,
    rtol = 1e-4,
    atol = 1e-4,
    maxevals = 100,
)




outv1=pmap(
    p0 -> FRGRealTime.propReintqs(
        p0,
        0.05,
        145.0,
        600.0,
        800.0,
        4.0,
        m2fun,
        lamfun,rtol = 1e-6,
        atol = 1e-6,
        maxevals = 200,
    ),v1
)


plot!(v1,outv1)
