

v1=SharedArray(collect(0.01:0.01:10.0))

v2=SharedArray(collect(20.0:5:60))


p0max=400.0
p0min=0.01

psmax=300.0
psmin=0.01

v1=SharedArray(exp.(collect(log(p0min):(log(p0max)-log(p0min))/100:log(p0max))))
v2=SharedArray(exp.(collect(log(psmin):(log(psmax)-log(psmin))/75:log(psmax))))


Imv1=[(i,j) for j in v2 for i in v1]


outvIm3 = pmap(
    p0 -> FRGRealTime.fastpropImintqs_All(
        pv,
        0.000001,
        30.0,
        1.0,
        800.0,
        4.0,
        m2fun,
        lamfun,
        rtol=1e-2,atol=1e-2,maxevals=80,initdiv=8,
    ),
    v1,
)

temp_T

outvIm4 = pmap(
    pv -> FRGRealTime.fastpropImintqs_All(
        pv[1],
        pv[2],
        temp_T,
        1.0,
        800.0,
        4.0,
        m2fun,
        lamfun,
        rtol=1e-10,atol=1e-10,maxevals=10000000,initdiv=40,
    ),
    Imv1,
)


plot(v1,outvIm4)

plot(v1[120:140],(outvIm4./(outvIm4.^2+outvRe.^2))[120:140],yaxis=:log)


outvIm4=pmap(p0->FRGRealTime.propImintqs_F4simple(p0, 0.0001, 40.0, 1.0, 800.0, 4.0, m2fun, lamfun),v1)



using Plots



plot(v1,outvIm3)
dir="/home/wjfu1/tyy/FRG-RealTime-data/Gamma2Im"
writedlm("$dir/propImI2p0+q0.dat",outvIm3)
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/ImI2p0+q0.dat",outvIm3)

vI1pmq=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/ImI1p0-q0.dat")
vI1ppq=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/ImI1p0+q0.dat")
vI2pmq=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/ImI2p0-q0.dat")
vI2ppq=readdlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Im/ImI2p0+q0.dat")

plot(v1,vI1pmq,label="I1 p-q",legend=:topleft)
plot!(v1,vI1ppq,label="I1 p+q",legend=:topleft)
plot!(v1,vI2pmq,label="I2 p-q",legend=:topleft)
plot!(v1,vI2ppq,label="I2 p+q",legend=:topleft)

plot!(v1,vI1pmq+vI1ppq+vI2pmq+vI2ppq,label="Total",legend=:topleft)




FRGRealTime.fastpropImintqs_F4(
    10.00,
    0.000001,
    145.0,
    1.0,
    800.0,
    4.0,
    m2fun,
    lamfun,
    rtol=1e-2,atol=1e-2,maxevals=10,initdiv=200,
)









outvIm4=pmap(
    p0 -> FRGRealTime.fastpropImintqs_All(
        p0,
        0.000001,
        100.0,
        20.0,
        100.0,
        4.0,
        m2fun,
        lamfun,rtol=1e-3,atol=1e-3
    ),
    v2,
)-pmap(
    p0 -> FRGRealTime.fastpropImintqs_F4(
        p0,
        0.000001,
        100.0,
        20.0,
        100.0,
        4.0,
        m2fun,
        lamfun,rtol=1e-3,atol=1e-3
    ),
    v2
)


outvIm5=pmap(
    p0 -> FRGRealTime.propImSimple_dq0(
        p0,
        0.000001,
        100.0,
        20.0,
        100.0,
        4.0,
        m2fun(200.0),
        -10.0
    ),
    v2,
)
FRGRealTime.



outvIm6=pmap(
    p0 -> FRGRealTime.fastpropImintqs_All_dVdq0(
        p0,
        0.000001,
        100.0,
        20.0,
        100.0,
        4.0,
        m2fun,
        lamfun,rtol=1e-3,atol=1e-3
    ),
    v2,
)



FRGRealTime.propImSimple_dq0(
    40.0,
    0.000001,
    100.0,
    40.0,
    100.0,
    4.0,
    m2fun(200.0),
    -10.0
)


plot(k->3 * Epi(k, m2fun(k)),20,100)


using Plots

plot(v2,outvIm4)
plot(v2,outvIm5)
plot!(v2,outvIm6)



FRGRealTime.deltasumkfix(100.0, 0.0001, 100.0, 100.0, 4.0,1.0, 800.0, m2fun, lamfun)

FRGRealTime.propImSimple_dq0(
    10.0,
    0.000001,
    100.0,
    10.0,
    200.0,
    4.0,
    10.0^2,
    -10.0
)




FRGRealTime.VImintqsSimple_dq0(10.0, 0.0001, 10.0, 100.0, 4.0, 100.0, -10.0,200.0)

FRGRealTime.VImSimple(10.0, 0.0001, 20.0, 10.0, 100.0, 100.0, 4.0, -10.0, 200.0)


FRGRealTime.loopfunppSmooth(10.0, 0.0001, 20.0, 100.0, 100.0)



plot(v1,outvIm3)
plot!(v1,outvIm3+outvIm4)
plot!(v1,outvIm4)


pmap(
    p0 -> FRGRealTime.propImintqs(
        p0,
        0.000001,
        100.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),
    v1,
)

using Plots

outvIm3=map(
    p0 -> FRGRealTime.propImintqs_F4simple(
        p0,
        0.000001,
        100.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),
    v1,
)


outvIm4=map(
    p0 -> FRGRealTime.fastpropImintqs_F4(
        p0,
        0.000001,
        100.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,rtol=1e-3,atol=1e-3
    ),
    v1,
)


FRGRealTime.fastpropImintqs_F4(
    15.0,
    0.000001,
    100.0,
    1.0,
    50.0,
    4.0,
    m2fun,
    lamfun,rtol=1e-4,atol=1e-4,initdiv=200,maxevals=100000000
)

FRGRealTime.propImintqs_F4simple(
    15.0,
    0.000001,
    100.0,
    1.0,
    50.0,
    4.0,
    m2fun,
    lamfun,
)




using Plots


plot(x->sin(x),1.0,10.0)

plot(v1, outvIm3, xaxis = "p0", label = "T=145,Λ=400")
plot!(v1, outvIm4, xaxis = "p0", label = "T=145,Λ=400")



plot(v1, outvIm2+outvIm, xaxis = "p0", label = "T=145,Λ=400")


FRGRealTime.fastpropImintqs_All(
    600.0,
    0.000001,
    100.0,
    1.0,
    1600.0,
    4.0,
    m2fun,
    lamfun,atol=1e-3,rtol=1e-3
)


plot(k->FRGRealTime.VImintqs(
    20.0,
    0.0001,
    k,
    100.0,
    4.0,
    1.0,
    1600.0,
    m2fun,
    lamfun,
),1.0,1600.0)



FRGRealTime.VImintqsSimple(
    335.0,
    0.0001,
    10.0,
    100.0,
    4.0,
    m2fun,
    lamfun,
    1600.0,
)
