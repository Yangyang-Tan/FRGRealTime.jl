
plot(v1,2*outvIm./(outvIm.^2+outvRe.^2),xaxis="p0",yaxis=:log)
using Dierckx

Spline1D(readdlm("$dir3/Imp0.dat")[:,1],readdlm("$dir3/ImGamma2.dat")[:,1])
