using FRGRealTime,Plots

u0 = [-8.0, 116585.3]

sol=FRGRealTime.ZeroLPASolve(
    145.0,
    4.0,
    1.0,
    800.0,
    u0 = [-8.0, 116585.3],
)

-sol[end][2]


writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_lpa_m2k.dat",-hcat(sol.u...)[2,:])
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_lpa_lamdak.dat",2*hcat(sol.u...)[1,:])
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_lpa_k.dat",sol.t)
