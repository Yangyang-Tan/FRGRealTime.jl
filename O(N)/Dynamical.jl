using Distributed
using SharedArrays
addprocs(8)
nprocs()
@everywhere using FRGRealTime, DelimitedFiles,Dierckx,Plots
@everywhere kdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_zero_k.dat")[:,1]|>reverse
@everywhere lamdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_zero_lamdak.dat")[:,1]|>reverse
@everywhere m2data=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_zero_m2k.dat")[:,1]|>reverse
@everywhere lamfun=Spline1D(kdata,lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)

m2fun
p0out=SharedArray(collect(0.0001:0.0001:0.1))

p0max=1e-3
p0min=1e-5
psmax=30.0
psmin=1.0
p0out=SharedArray(exp.(collect(log(p0min):(log(p0max)-log(p0min))/20:log(p0max))))
psout=SharedArray(exp.(collect(log(psmin):(log(psmax)-log(psmin))/80:log(psmax))))



dyv1=[(i,j) for j in psout for i in p0out]
outvdy = pmap(
    pv -> FRGRealTime.fastpropImintqs_All(
        pv[1],
        pv[2],
        20.43287,
        1.0,
        800.0,
        4.0,
        m2fun,
        lamfun,
        rtol=1e-10,atol=1e-10,maxevals=20000000,initdiv=180,
    ),
    dyv1,
)


plot(p0out,outvdy)


writedlm(path,data)

plot(psout,-(outvdy)./psout.^2,xaxis=:log,yaxis=:log,legend=:topleft)
plot(p0out,outvdy)


dir="/home/wjfu1/tyy/FRG-RealTime-data/DynamicalCE"
writedlm("$dir/dyp0_T=21.43287.dat",p0out)

writedlm("$dir/dyps_T=20.43287.dat",psout)
writedlm("$dir/dyp0_T=20.43287.dat",p0out)

writedlm("$dir/dy_T=20.43287.dat",outvdy)



dyfun=Spline1D(p0out,outvdy)
plot(p0out,outvdy)

plot(x->derivative(dyfun,x),p0min,p0max)



p0out=SharedArray(collect(0.0001:0.0001:0.001))

psout=SharedArray(exp.(collect(log(0.002):(log(0.1)-log(0.002))/80:log(0.1))))

@everywhere myfun(p0, ps) =
    FRGRealTime.fastpropImintqs_All(
        p0,
        ps,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ) + FRGRealTime.fastpropImintqs_All_dq0(
        p0,
        ps,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,rtol=1e-4,atol=1e-4
    )


propp0out4=pmap(ps->myfun.(p0out, ps),psout)

writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Λ=400_p0_0.001-0.01_ps.dat",p0out)
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Λ=400_p0_ps_0.002-1.0.dat",psout)
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Λ=400_Gamma2p0_0.001-0.01_ps_0.002-1.0_all.dat",propp0out4)



outvIm = pmap(
    p0 -> FRGRealTime.fastpropImintqs_All(
        p0,
        0.000001,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),
    v1,
)

plot(v1,outvIm)



plot(v1,(p0->FRGRealTime.deltasumImprop1(
    p0,
    0.05,
    145.0,
    1.0,
    400.0,
    4.0,
    m2fun,
    lamfun,
)).(v1).+(p0->FRGRealTime.deltasumImprop2(
    p0,
    0.05,
    145.0,
    1.0,
    400.0,
    4.0,
    m2fun,
    lamfun,
)).(v1))



outv1|>maximum


plot(p0->FRGRealTime.flowpm(p0,10.0,20.0,-10.0,145.0),0.0,10.0)

plot(p0->FRGRealTime.flowpp_intcostheqs(p0,0.00,100.0,100.0,m2fun(100.0),145.0),1.0,200.0)
