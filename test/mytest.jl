using FRGRealTime
using HCubature
using Plots
δk=0.02
ppfunps_delta(p0, k, kprim, m, T)

plot(p0 -> 2*ppfunps_delta(p0, 10.0, 10.0, -2.0, 145.0), 19.7, 19.9)
plot(
    p0 -> FRGRealTime.delta2_intcosthqs(p0, 0.001, 10.0, 10.0, -2.0, 145.0),
    19.7, 19.9
)
