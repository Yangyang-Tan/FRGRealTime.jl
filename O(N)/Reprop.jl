using Distributed
using Plots
addprocs(24)
nprocs()
using SharedArrays
@everywhere using FRGRealTime, DelimitedFiles, Dierckx

kdata=sol.t[1:1:end]|>reverse
lamdata=sol[1,:][1:1:end]|>reverse
m2data=-sol[2,:][1:1:end]|>reverse

@eval @everywhere kdata=$kdata
@eval @everywhere lamdata=$lamdata
@eval @everywhere m2data=$m2data

@everywhere kdata=collect(1:1.0:1000)
@everywhere lamdata=-10*one.(kdata)
@everywhere m2data=10.0^2 *one.(kdata)


@everywhere lamfun=Spline1D(kdata,lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)
v1=SharedArray(collect(0.1:4:50.0))



@everywhere kdata=readdlm("/home/wjfu1/tyy/CTP-fRG-Test/Real-Time-data/realtime_zero_k.dat")[:,1]|>reverse
@everywhere lamdata=readdlm("/home/wjfu1/tyy/CTP-fRG-Test/Real-Time-data/realtime_zero_lamdak.dat")[:,1]|>reverse
@everywhere m2data=readdlm("/home/wjfu1/tyy/CTP-fRG-Test/Real-Time-data/realtime_zero_m2k.dat")[:,1]|>reverse
@everywhere lamfun=Spline1D(kdata,0.5*lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)






plot(x->m2fun(x),1.0,800.0)
plot(x->lamfun(x),1.0,800.0)


p0v1=SharedArray(collect(0.01:0.5:400.0))
psv1=SharedArray(collect(0.01:2.0:300.0))


Rev1=[(i,j) for j in psv1 for i in p0v1]



@time outvRe2=pmap(
    pv -> FRGRealTime.fastpropReintqs_All(
        pv[1],
        pv[2],
        temp_T,
        1.0,
        800.0,
        4.0,
        m2fun,
        lamfun,
        rtol = 1e-5,
        atol = 1e-10,
        maxevals=100000000,
        initdiv=20,
    ),Imv1
)
# interrupt(collect(1:200))
# p0v1=readdlm("$dir/p0.dat")[:,1]
# outvRe=readdlm("$dir/Re.dat")[:,1]


FRGRealTime.fastpropReintqs_All(
        1.01,
        0.000001,
        20.44,
        1.0,
        800.0,
        4.0,
        m2fun,
        lamfun,
        rtol = 1e-4,
        atol = 1e-5,
        maxevals=10000000,
        initdiv=40,
    )


plot(p0v1,outvRe2)

interrupt(collect(1:100))



plot!(p0v1,p0v1.^2)


interrupt(collect(1:400))

dir="/home/wjfu1/tyy/FRG-RealTime-data/Gamma2Im"
dir2="/home/wjfu1/tyy/FRG-RealTime-data/Re/ps=0"

# writedlm("$dir/Rep0.dat",v1)
# writedlm("$dir/ReGamm2_ps=0.dat",outvRe)

Imgamm2fun_All=Spline1D(
    readdlm("$dir/Imp0.dat")[:, 1],
    readdlm("$dir/propImI1p0+q0.dat")[:, 1] +
    readdlm("$dir/propImI1p0-q0.dat")[:, 1] +
    readdlm("$dir/propImI2p0+q0.dat")[:, 1] +
    readdlm("$dir/propImI2p0-q0.dat")[:, 1],
)


Imgamm2fun_1=Spline1D(
    readdlm("$dir/Imp0.dat")[:, 1],
    readdlm("$dir/propImI1p0+q0.dat")[:, 1]
)

Imgamm2fun_2=Spline1D(
    readdlm("$dir/Imp0.dat")[:, 1],
    readdlm("$dir/propImI1p0-q0.dat")[:, 1]
)
Imgamm2fun_3=Spline1D(
    readdlm("$dir/Imp0.dat")[:, 1],
    readdlm("$dir/propImI2p0+q0.dat")[:, 1]
)
Imgamm2fun_4=Spline1D(
    readdlm("$dir/Imp0.dat")[:, 1],
    readdlm("$dir/propImI2p0-q0.dat")[:, 1]
)



Regamm2fun=Spline1D(readdlm("$dir2/Rep0.dat")[:, 1],
    readdlm("$dir2/ReGamm2_ps=0.dat")[:, 1],
)

plot(x->Imgamm2fun_1(x),1.0,50.0)




plot(x->-(Imgamm2fun_3(x)+Imgamm2fun_4(x))/((Imgamm2fun_3(x)+Imgamm2fun_4(x))^2+Regamm2fun(x)^2),1.0,400.0,yaxis=:log)


rho(x)=Imgamm2fun_All(x) / (Imgamm2fun_All(x)^2 + Regamm2fun(x)^2)

mv1=collect(1.0:0.1:400.0)

mv2=collect(202.9:0.1:400.0)

plot(
    mv1,rho.(mv1),
    yaxis = :log,
)

plot!(
    mv2,sign.(rho.(mv2)).*log.(abs.(rho.(mv2)))
    # yaxis = :log,
)



sign(x)*log(1+abs(x)/10^C)


plot(v1,outvRe/145^2,legend=:topleft)
