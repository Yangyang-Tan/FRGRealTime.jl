using Distributed
using SharedArrays
addprocs(8)
nprocs()
@everywhere using FRGRealTime, DelimitedFiles,Dierckx,Plots
@everywhere kdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/k_pde_break.dat")[:,1]|>reverse
@everywhere lamdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/lamdak_pde_break.dat")[:,1]|>reverse
@everywhere m2data=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/m2k_pde_break.dat")[:,1]|>reverse
@everywhere lamfun=Spline1D(kdata,lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)
v1=SharedArray(collect(0.1:3.0:200.0))
@everywhere FRGRealTime.propReintqs(4.0, 10.0, 145.0,1.0,400.0, 4.0, m2fun, lamfun)


FRGRealTime.propReintqs(100.0, 0.05, 145.0,1.0,400.0, 4.0, m2fun, lamfun,maxevals=500)

outv1=pmap(
    p0 -> FRGRealTime.propReintqs(
        p0,
        0.05,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),v1
)

plot(v1,outv1)
