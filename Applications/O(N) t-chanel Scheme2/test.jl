solfoursimple=tchanelSolveFourPointParallelSimple(
    0.0,
    Tself,
    lpakrang,
    solmassfun,
    sollambdafun,
    Imlambdaflow3,
    Relambdaflow3,
    0.1,
    u0 = [0.00000, 6 * 8.0],
    config=config_spec,
)
config_spec

-solself.u[end][1]

q0grid[1]

solfoursimplefun2=Spline1D(reverse(solfoursimple.t),reverse(hcat(solfoursimple.u...)[2,:]))

lines(0.1:0.1:800.0,sollambdafun.(0.1:0.1:800.0).-solfoursimplefun2(0.1:0.1:800.0))
lines(solself.t,hcat(solself.u...)[2,:])
solself.t
lines!(solfoursimple.t,hcat(solfoursimple.u...)[2,:])
current_figure()

solmassfun(0.0)

tv1=0.1:0.1:800.0
tv2=similar(tv1)
vmapntt(solmassfun,tv1)
@tturbo for i ∈ eachindex(tv1)
    tv2[i]=solmassfun(tv1[i])
end


testf1(x,y)=Relambdaflow3(0.0, x, 100.0, -10.0, 10.0, 200.0, 0.0, y, ϵ5)

testf1(x,y)=sign(x-y)
vmapntt!(flowf, lambdaq0dk, q0grid, lambdaq0)
vmapntt(testf1, tv2,tv1)
