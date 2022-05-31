using FastGaussQuadrature
x, w = gausslegendre(3)

function ImRewarp(
    Imdst,
    Redst,
    imtex,
    retexA,
    retexB,
    p0min,
    q0min,
    kmin,
    kmax,
    dp0,
    dq0,
    dk,
    p0grid,
    q0grid,
    k,
    m2,
    Ek,
    Temper,
    ϵ = 1.0f0,
    Nf = 4.0f0,
)
    tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    I = CartesianIndices(Imdst)
    @inbounds if tid <= length(I)
        i, j, l= Tuple(I[tid])
        Imdst[i, j, l] = w[l]*Imlambdaflow_Texture_3d(p0grid[i], q0grid[j],p0min,q0min,kmin,dp0,dq0,dk, k[l], m2,Ek, Temper, imtex, retexA,retexB, ϵ, Nf)
        Redst[i, j,l] = w[l]*Relambdaflow_Texture_3d(p0grid[i], q0grid[j],p0min,q0min,kmin,dp0,dq0, dk,k[l], m2,Ek, Temper, imtex, retexA,retexB, ϵ, Nf)
    end
    return
end



testv=CUDA.randn(800,800,100)
testv2=CUDA.randn(800,800,100)

for i in 2:100
    testv2[:,:,i].=testv2[:,:,i].+testv2[:,:,i-1]
end


tx,tw=gausslobatto(10)

gausslobatto(10)|>println
testf(x)=x^2+cos(x)
tw[2:end]' *testf.(tx[2:end])
