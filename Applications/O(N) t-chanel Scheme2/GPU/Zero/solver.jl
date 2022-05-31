function tchanelSolveFourPoint_Zero_Zero(
    k0::T,
    Temper::T,
    krang::RGscale,
    ϵ = 1.0;
    Nf = 4,
    u0,
    config::ODEConfig,
) where {T<:Number}
    function flow_tZero_Zero(du,u, p, k)
        Ek = k
        # Texture Interpolation start
        println("k=", k)
        Relambdaflow_Zero(0.0, 0.0, k, 0.0, Temper, u, Nf, ϵ)
        Twopointflow_Zero(lambdamp0Ek,lambdamp0mEk,lambdamp0p0, k, m2, T, Nf = 4.0f0)
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tZero_Zero, u0, (krang.UV, k0))
    return solve(
        prob,
        alg,
        # dt=0.5,
        # dtmax = 0.5f0,
        abstol = config.atol,
        reltol = config.rtol,
        # dt=0.0000001,
        dtmin = config.dtmin,
        # adaptive = config.adaptive,
        progress = config.progress,
    )
end

config_spec.atol = 1e-10
config_spec.rtol = 1e-10
sol_zero_zero = tchanelSolveFourPoint_Zero_Zero(
    1.0,
    145.0,
    lpakrangGPU,
    1.0,
    u0 = 24.0,
    config = config_spec,
)

plot!(sol_zero_zero.t, sol_zero_zero.u)





function tchanelSolveFourPoint_Zero(
    k0::T,
    Temper::T,
    krang::RGscale,
    # lambdafun,
    ϵ = 1.0f0;
    dp0 = 0.125f0,
    dq0 = 0.125f0,
    Nf = 4.0f0,
    m2ini = -0.18216546875*lpakrang.UV^2,
    margin::Int = 0,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -krang.UV / 1.0f0 - margin * dp0
    p0max = krang.UV / 1.0f0 + margin * dp0
    q0min = -krang.UV / 1.0f0 - margin * dq0
    q0max = krang.UV / 1.0f0 + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = cu(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
    Imlambdamp0p0=similar(p0grid)
    Relambdamp0p0=similar(p0grid)
    idxzero = Int((length(p0grid) + 1) / 2)
    mytimer=1
    function flow_tZero(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]

        m2 = -ReGamm2[idxzero]
        Relambda0=Relambda[idxzero,idxzero]


        Ek = Eb(k,m2)
        idxEk = round(Int32, max(min((Ek - q0min) / dq0 + 2.0f0, length(q0grid)+1), 1))
        idxmEk = round(Int32, max(min((-Ek - q0min) / dq0 + 2.0f0, length(q0grid)+1), 1))



        mImlambdamp0Ek = @view u[:, idxmEk, 1]
        mImlambdamp0mEk = @view u[:, idxEk, 1]
        Imlambdamp0p0 .= -diag(Imlambda)
        Relambdamp0Ek =@view u[:, idxEk, 2]
        Relambdamp0mEk = @view u[:, idxmEk, 2]
        Relambdamp0p0 .= diag(Relambda)



        dkImGamm2 .= Twopointflow_Zero.(-mImlambdamp0Ek,-mImlambdamp0mEk,Relambda0, k, m2, Temper, Nf)
        dkReGamm2 .= Twopointflow_Zero.(Relambdamp0Ek,Relambdamp0mEk,Relambda0, k, m2, Temper, Nf )


        dkImlambda .= Imlambdaflow_Zero.(p0grid, q0grid, k, m2, Temper, Relambda0, Nf, ϵ)
        dkRelambda .= Relambdaflow_Zero.(p0grid, q0grid, k, m2, Temper, Relambda0, Nf, ϵ)

        # Texture Interpolation start
        println("k=", k,"idxzero=",idxzero)
        # Relambdacpu=Array(Relambda)[:, idxEk]
        # plot(-400.0:0.5:400.0, @. RelambdaAcpu^2/RelambdaBcpu^2)
        # plot(p0min:dp0:p0max, @. RelambdaAcpu/(RelambdaBcpu))|>display
        # display(plot!(-400.0:0.5:400.0, @. RelambdaAcpu^2/RelambdaBcpu^2,ylims =(-1,5000)))
        # plot(p0min:dp0:p0max, Imlambdacpu1)
        # display(plot!(p0min:dp0:p0max, Imlambdacpu2))
        # Imlambdacpu1 = Array(Relambda)[:, 1]
        # Imlambdacpu2 = Array(Relambda)[:, end]
        # plot(p0min:dp0:p0max, Imlambdacpu1)
        # display(plot!(p0min:dp0:p0max, Imlambdacpu2))
        mytimer=mytimer+1
        if mod1(mytimer,10)==1
           display(plot(Array(Relambdamp0Ek)))
        end
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tZero, u0, (krang.UV, k0))
    return solve(
        prob,
        alg,
        # dt=0.5,
        # dtmax = 0.5f0,
        abstol = config.atol,
        reltol = config.rtol,
        calck = false,
        # dt=0.0000001,
        dtmin = config.dtmin,
        # adaptive = config.adaptive,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_end = config.save_end,
        progress = config.progress,
    )
end




config_spec.atol = 1e-5
config_spec.rtol = 1e-5

lambdafun_zero=Spline1D(sol_zero_zero.t|>reverse, sol_zero_zero.u|>reverse)

sol_zero=tchanelSolveFourPoint_Zero(
    1.0f0,
    145.0f0,
    lpakrangGPU,
    1.0f0,
    dp0 = 0.5f0,
    dq0 = 0.5f0,
    Nf = 4.0f0,
    margin = 0,
    config = config_spec,
)

plot(p0->Relambdaflow_Zero(p0,0.0f0,10.0,-10.0,145.0f0,24.0f0,4.0f0,1.0f0),-100.0,100.0)



plot(-10.0:0.5f0:10.0,Relambdaflow_Zero.(-10.0:0.5f0:10.0,0.0f0,10.0,-10.0,145.0f0,24.0f0,4.0f0,0.1f0))

plot(-10.0:0.5f0:10.0,TF.ReFb2.(-10.0:0.5f0:10.0,10.0f0,-10.0f0,145.0f0,0.1f0))



TF.ReFb2(0.0f0,10.0f0,-10.0f0,145.0f0,0.1f0)

TF.ReFb2(0.01f0,10.0f0,-10.0f0,145.0f0,0.1f0)



plot(-800:0.5:800,Array(sol_zero[end][:,1,2]))


img2=Array(sol_zero[end][1602:end,1,1])
reg2=Array(sol_zero[end][1602:end,1,2])
spec = @. img2/(img2^2+reg2^2)


plot(spec[1:1400], yaxis=:log)

plot(img2[1:400])


-sol_zero[end][1601,1,2]







heatmap(Array(sol_zero[end][:,:,2])')

heatmap(-200:0.125:200,-200:0.125:200,Array(sol_zero[end][:,:,2])')
plot()
