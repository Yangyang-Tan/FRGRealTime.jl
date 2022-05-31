100

(100 +160) / 0.125 + 1.0f0

Array(solselfGPU[end][:,2081,2])
plot(-160:0.125:160,Array(solselfGPU[end][:,2081,2]), labels="Full")
plot!(-160:0.125:160,Array(sol_zero[end][:,2081,2]),labels="ReadOff",title="Relambda(p_0,Ek)")

plot(-160:0.125:160,Array(solselfGPU[end][:,2081,1]), labels="Full")
plot!(-160:0.125:160,Array(sol_zero[end][:,2081,1]),labels="ReadOff",title="Imlambda(p_0,Ek)")




plot(-800:0.5:800,Array(sol_zero[end][:,1001,2]))


plot(-100:0.25:100,Array(solselfGPU[end][:,100,1]))
plot!(-100:0.25:100,-Array(solselfGPU[end][:,end-100,1]))


plot(-200:0.125:200,Array(solselfGPU[end][100,:,2]))

plot(0:0.125:200,Array(solselfGPU[end][400,1:1601,2])|>reverse)
plot(0:0.5:400,Array(solselfGPU[end][1300,801:end,2])-(Array(solselfGPU[end][1300,1:801,2])|>reverse))

plot(-400:4.0:400,Array(solselfGPU[end][500,801:end,1])+(Array(solselfGPU[end][end-500,1:801,1])|>reverse))




plot(-400:0.5:400,Array(solselfGPU3[end][:,801-200,1]))


plot(-100:1.0:100,Array(solselfGPU3[end][:,1,2]))



plot!(-400:0.5:400,Array(solselfGPU2[end][end,:,2]))


plot(0:0.5:400,Array(solselfGPU[end][end,1:801,1])|>reverse)

plot(0:0.5:400,Array(solselfGPU2[end][end,801:end,1]))
plot!(0:0.5:400,Array(solselfGPU[end][end,801:end,1]))

plot!(0:0.5:400,-Array(solselfGPU[end][1,801:end,1]))



plot(0:4.0:400,Array(solselfGPU[end][end,201:end,2])+(Array(solselfGPU[end][end,1:201,2])|>reverse))

20/0.125

plot(-200:0.125:200,-Array(solselfGPU[end][:,1600-160,1])-Array(solselfGPU[end][:,1600+160,1]))


plot(-200:0.125:200,Array(solselfGPU[end][:,1600-160,2])+Array(solselfGPU[end][:,1600+160,2]))





plot!(-400:0.5:400,Array(solselfGPU[end][end-100,:,1]))


plot!(-200:0.125:200,Array(solselfGPU[end][end-100,:,1]))


heatmap((solselfGPU.u[end][:,:,1]|>Array)')

p0min=-800.0f0
p0max=800.0f0
q0min=-800.0f0
q0max=800.0f0
dp0=80.0f0
dq0=80.0f0
k=790.0f0
m2=0.0f0
Ek=k
Temper=145.0f0
p0grid = collect(p0min:dp0:p0max) |> cu
q0grid = collect(q0min:dq0:q0max)' |> cu
u0 = cu(init_fourpoint_2d_test(p0grid, q0grid,margin=0))
ImlambdaTextureMem = CuTextureArray(u0[:, :, 1])
RelambdaTextureMem = CuTextureArray(u0[:, :, 2])



ImlambdaTextureMem = CuTextureArray(u0[:, :, 1])
RelambdaTextureMem = CuTextureArray(u0[:, :, 2])


Imlambdatex =CuTexture(ImlambdaTextureMem; interpolation = CUDA.LinearInterpolation())
Relambdatex =CuTexture(RelambdaTextureMem; interpolation = CUDA.LinearInterpolation())

dkImlambda=similar(u0[:, :, 1])
dkRelambda=similar(u0[:, :, 2])

ImRewarp2!(
            dkImlambda,
            dkRelambda,
            Imlambdatex,
            Relambdatex,
            p0min,
            q0min,
            dp0,
            dq0,
            p0grid,
            q0grid,
            k,
            0.0f0,
            Ek,
            Temper,
            5.0,
        )

plot(Array(reverse(dkRelambda[1,1:11])-dkRelambda[1,11:end]))



plot(Array(reverse(dkImlambda[1,1:11])+dkImlambda[end,11:end]))


heatmap(Array(dkRelambda))
heatmap(Array(dkImlambda))

TF.ImFb2(-2*Ek, k, 0.0f0, 145.0f0, 4.0f0)
TF.ImFb2(2*Ek, k, 0.0f0, 145.0f0, 4.0f0)

TF.ImFb2(0.0f0, k, 0.0f0, 145.0f0, 4.0f0)


plot(Array(reverse(u0[1,1:11,1])+u0[end,11:end,1]))
heatmap(Array(u0[:,:,1]))
