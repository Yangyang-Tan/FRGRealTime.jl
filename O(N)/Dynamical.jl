using FRGRealTime
using DelimitedFiles
using Dierckx
using Plots
kdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/k_pde_break.dat")[:,1]|>reverse
lamdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/lamdak_pde_break.dat")[:,1]|>reverse
m2data=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/m2k_pde_break.dat")[:,1]|>reverse
lamfun=Spline1D(kdata,lamdata)
m2fun=Spline1D(kdata,m2data)

plot(p0->FRGRealTime.propImsimpleintqs(p0, 10.0, 145.0,1.0,400.0, 4.0, m2fun, lamfun),1.0,100.0)
