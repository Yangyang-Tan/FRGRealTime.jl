using Distributed
using SharedArrays,Plots
addprocs(48)
nprocs()
@everywhere using FRGRealTime,QuadGK,Dierckx,DelimitedFiles



# u0 = [-8.0, 188498.0]
-20.0, 136 * 188498.0
v1 = SharedArray(collect(25:1.0:400.0))
v2 = SharedArray(collect(145:1.0:400.0))

@everywhere function mreturnfun(T)
    sol_temp = FRGRealTime.ZeroSolve(
        T,
        4.0,
        1.0,
        800.0,
        u0=[-8.0, 116585.3],
        # dtmax = 5.0,
        atolint=1e-4,
        rtolint=1e-4,
    )
    println("T=", T)
    return sqrt(-sol_temp[end][2])
end



@everywhere function mreturnfun2(T)
    sol_temp = FRGRealTime.ZeroLPASolve(
        T,
        4.0,
        1.0,
        800.0,
        u0=[-8.0, 116585.3],
    )
    println("T=", T)
    return sqrt(-sol_temp[end][2])
end

outv1 = pmap(mreturnfun, v1)
outv2 = pmap(mreturnfun2, v2)

plot!(v1,outv1)

dir = "/home/wjfu1/tyy/FRG-RealTime-data/Zero/ThermalMass"
writedlm("$dir/realtime_T.dat",v1)
writedlm("$dir/LPA_T.dat",v2)
writedlm("$dir/realtime_m2.dat",outv1)
writedlm("$dir/LPA_m2.dat",outv2)

116585.3
13108
20.43287
temp_T=80.0
sol = FRGRealTime.ZeroSolve(
    temp_T,
    4.0,
    1.0,
    800.0,
    u0=[-8.0, 116585.3],
    # dtmax = 5.0,
    # dtmax=0.2,
    atolint=1e-6,
    rtolint=1e-6,
)


plot(sol.t,-3*sol[1,:])
plot(sol.t,-sol[2,:])

<<<<<<< HEAD
plot(sol.t,-sol[1,:])

plot(sol.t,-sol[2,:])
=======
-sol[2,:][end]

-sol[2,:][end]
>>>>>>> c1b5f5caeab5b4fe6926368d82a33f6b833eefa9

plot(k->sqrt(-sol(k)[2]+k^2),1.0,2.0)


sqrt(-sol[2,:][end])

sol = FRGRealTime.ZeroLPASolve(
        145.0,
        4.0,
        1.0,
        800.0,
        u0=[-8.0, 116585.3],
)


<<<<<<< HEAD
sol=FRGRealTime.ZeroLPASolve(
    145.0,
    4.0,
    1.0,
    800.0,
    u0 = [-8.0, 116585.3],
)

sol=FRGRealTime.ZeroSolve(
    100.0,
    4.0,
    1.0,
    800.0,
    u0 = [-8.0, 116585.3],
    # dtmax = 5.0,
    atolint = 1e-4,
    rtolint=1e-4,
)




sqrt(-sol[end][2])
=======
sol = FRGRealTime.ZeroLPASolve(
        145.0,
        4.0,
        1.0,
        800.0,
        u0=[-14.0, 116585.3],
)

>>>>>>> c1b5f5caeab5b4fe6926368d82a33f6b833eefa9


-sol[2,:][end]

1
sol2 =
    FRGRealTime.ZeroSolve(145.0, 4.0, 1.0, 800.0, u0=sol(800.0), dtmax=0.5)

plot(sol.t,-sol[1,:])

<<<<<<< HEAD
-quadgk(
    k -> FRGRealTime.propReZeroflow(k, 60.0, 4.0, sol, 800.0, atol = 1e-6, rtol = 1e-6),
    1.0,
    800.0,
    atol = 1e-8,
    rtol = 1e-8,
)[1]+188498.0
xaxis
=======
plot!(v1,outv1)
plot!(v1,outv2)
>>>>>>> c1b5f5caeab5b4fe6926368d82a33f6b833eefa9

plot(sol.t,sol[1,:])

<<<<<<< HEAD
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Zero/realtime_zero_T=100_m2k.dat",-hcat(sol.u...)[2,:])
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Zero/realtime_zero_T=100_lamdak.dat",2*hcat(sol.u...)[1,:])
writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Zero/realtime_zero_T=100_k.dat",sol.t)
=======
dir3
dir3 = "/home/wjfu1/tyy/FRG-RealTime-data/SpectralFunction/ps=0/LPAT=145"
writedlm("$dir3/m2k.dat",-hcat(sol.u...)[2,:])
writedlm("$dir3/lamdak.dat",2 * hcat(sol.u...)[1,:])
writedlm("$dir3/k.dat",sol.t)
writedlm("$dir3/Imp0_small.dat",v1)
writedlm("$dir3/ImGamma2_small.dat",outvIm4)
writedlm("$dir3/Rep0.dat",p0v1)
writedlm("$dir3/Reps.dat",psv1)
writedlm("$dir3/ReGamma2.dat",outvRe2)

>>>>>>> c1b5f5caeab5b4fe6926368d82a33f6b833eefa9

writedlm("$dir3/Imp0.dat",v1)
writedlm("$dir3/ImI2p0-q0.dat",outvIm4)





plot!(sol2.t,-sol2[2,:])
plot!(sol2.t,sol2[1,:])


sol = FRGRealTime.ZeroLPASolve(
    100.0,
    4.0,
    1.0,
    800.0,
    u0=[-8.0, 0.8423 * 116585.3],
)

sqrt(-sol[end][2])




m2fun(x) = -sol2(x)[2]
lamdafun(x) = sol2(x)[1]



k_temp = collect(1:1.0:800)
qs_temp = collect(collect(1:1.0:800)')

Voutvec = FRGRealTime.lamdabar_Zero.(
    Epi.(k_temp, m2fun.(k_temp)),
    qs_temp,
    k_temp,
    145.0,
    4.0,
    m2fun,
    lamdafun,
    800.0;
    atol=1e-6,
    rtol=1e-6,
)


k_temp = 800.0
qs_temp = 0.00001

using FRGRealTime

FRGRealTime.lamdabar_Zero(
    Epi.(k_temp, m2fun.(k_temp)),
    qs_temp,
    k_temp,
    145.0,
    4.0,
    m2fun,
    lamdafun,
    800.0;
    atol=1e-6,
    rtol=1e-6,
)


writedlm("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/realtime_zero_vertex.dat",Voutvec)


@eval @everywhere sol2 = $sol2
@everywhere m2fun(x) = -sol2(x)[2]
@everywhere lamdafun(x) = sol2(x)[1]

@everywhere k_temp = collect(1:20.0:800)
qs_temp = SharedArray(collect(collect(1:20.0:800)'))


Voutvec = map((k, qs) -> FRGRealTime.lamdabar_Zero(
    Epi.(k, m2fun.(k)),
    qs,
    k,
    145.0,
    4.0,
    m2fun,
    lamdafun,
    800.0;
    atol=1e-6,
    rtol=1e-6,
),k_temp,qs_temp)

plot(k_temp,outv1)





k -> lam4piflow(k, m2fun(k), 145.0, 4.0, lamdafun(k))
plot(sol2.t[end - 800:end],sol2[1,:][end - 800:end])
plot(sol2.t,sol2[1,:])
