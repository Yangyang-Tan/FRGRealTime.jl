function Coeffgamm2Simple_dq0(k, T, m)
    (
        k * (coth(Epi(k, m) / (2 * T)) / Epi(k, m)^2)
    ) / (16 * pi^2)
end


function propImSimple_dq0(p0, ps, T, IRScale, UVScale, Npi, m, lamda)
    -hquadrature(
        x ->
            2*Coeffgamm2Simple_dq0(x, T, m) *
            central_fdm(5,1)(z->VImintqsSimple(p0, ps,z, k, T, Npi, m, lamda,UVScale),Epi(k,m)),
        IRScale,
        UVScale,
        atol = 1e-3,
        rtol = 1e-3,initdiv=10,maxevals=1600,
    )[1]
end
