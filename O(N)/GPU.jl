using CUDA
using FRGRealTime
using Plots
vx1=collect(1:0.1:800)
vx2=collect(0.0:0.001:1.0)
cu(vx2)
using JuliaZH
sign(0)


FRGRealTime.dkVImintqs_gpu.(
    10.0,
    0.0001,
    Epi.(vx1, m2fun.(vx1))|>cu,
    vx1|>cu,
    vx2 .* (vx1 .- 800.0) .+ 800.0|>cu,
    m2fun.(vx2 .* (vx1 .- 800.0) .+ 800.0)|>cu,
    145.0,
    4.0,
    lamfun.(vx2 .* (vx1 .- 800.0) .+ 800.0)|>cu,
)

using FastGaussQuadrature

gausslegendre(1000)
cuv1,cuweight1=cu(gausslegendre(10000))
cuv2,cuweight2=cu.(collect.(transpose.(gausslegendre(10000))))

function domaintrans!(v1,v2,a,b)
    v1.=((b-a)/2).*v1.+((b-a)/2);
    v2.=((b-a)/2).*v2;
end

domaintrans!(cuv1,cuweight1,1.0,800.0)
domaintrans!(cuv2,cuweight2,0.0,1.0)

v2trans=cuv2 .* (cuv1 .- 800.0) .+ 800.0
lamvec=lamfun.(Array(v2trans))|>cu
mvec2=m2fun.(Array(v2trans))|>cu
mvec1=m2fun.(Array(cuv1))|>cu

function myfun5(w1,w2,v1,v2,mv1,mv2,lamv,p0,ps,Npi,T,UVScale)
    w1*w2*(v1 - UVScale) *2 *FRGRealTime.dkVImintqs_gpu(
        p0,
        ps,
        Epi(v1, mv1),
        v1,
        v2,
        mv2,
        T,
        Npi,
        lamv,
    )
end

function sumfun2(p0,ps)::Float32
    tempv=myfun5.(
        cuweight1,
        cuweight2,
        cuv1,
        v2trans,
        mvec1,
        mvec2,
        lamvec,
        p0,
        ps,
        4.0f0,
        145.0f0,
        800.0f0)
    return sum(tempv)
end

cup0vec=collect(1.0f0:0.1f0:4.9f0)
sumfun2.(cup0vec,0.001f0)

CUDA.memory_status()
CUDA.reclaim()


cpuweight1=cuweight1|>Array
cpuweight2=cuweight2|>Array
cpuv1=cuv1|>Array
cpuv2trans=v2trans|>Array
cpumvec1=mvec1|>Array
cpumvec2=mvec2|>Array
cpulamvec=lamvec|>Array

@btime(myfun5.(
    cpuweight1,
    cpuweight2,
    cpuv1,
    cpuv2trans,
    cpumvec1,
    cpumvec2,
    cpulamvec,
    10.0f0,
    0.001f0,
    4.0f0,
    145.0f0,
    800.0f0,
))

@time(1+1)

vx2 .* (vx1 .- 800.0) .+ 800.0|>cu

FRGRealTime.dkF1Allintqs_gpu.(
    100.0,
    10.001,
    10*vx2 |> cu,
    vx1 |> cu,
    m2fun.(vx1) |> cu,
    145.0,
)|>sum

function testfun2(x,y)::Float32
    sin(x)-y
end
cuvx1=vx1|>cu
cuvx2=vx2|>cu

function int4()
    mapreduce(x->testfun.(x,cuvx2),+,cuvx1)
end

mapreduce(x->testfun2.(x,cuvx2),+,cuvx1)


mapreduce(x->int3(x),+,cuvx2)


int2(1.0f0)

ps_temp=7.0
plot!(
    p0 -> FRGRealTime.loopfunpm(
        sqrt(p0^2+ps_temp^2-5^2),
        ps_temp,
        10.0,
        -80.0,
        10.0,
    ),
    0.0,
    100.0,
)
