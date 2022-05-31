
function thetafun(p0::Float32, Ek::Float32)
    if p0 > Ek && p0 > 200.0f0
        return 0.0f0
    else
        return 1.0f0
    end
end

# thetafun(abs(p0grid[i])+2*dp0,Ek)*thetafun(abs(q0grid[j])+2*dq0,Ek)*
# thetafun(abs(p0grid[i])+2*dp0,Ek)*thetafun(abs(q0grid[j])+2*dq0,Ek)*
function ImRewarp(
    Imdst,
    Redst,
    imtex,
    retex,
    p0min,
    q0min,
    dp0,
    dq0,
    p0grid,
    q0grid,
    k,
    m2,
    Ek,
    Temper,
    ϵ,
    Nf,
)
    tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    I = CartesianIndices(Imdst)
    @inbounds if tid <= length(I)
        i, j = Tuple(I[tid])
        Imdst[i, j] = Imlambdaflow_Texture(
            p0grid[i],
            q0grid[j],
            p0min,
            q0min,
            dp0,
            dq0,
            k,
            m2,
            Ek,
            Temper,
            imtex,
            retex,
            ϵ,
            Nf,
        )
        Redst[i, j] = Relambdaflow_Texture(
            abs(p0grid[i]),
            abs(q0grid[j]),
            p0min,
            q0min,
            dp0,
            dq0,
            k,
            m2,
            Ek,
            Temper,
            imtex,
            retex,
            ϵ,
            Nf,
        )
    end
    return
end

function ImReTwowarp(
    Imdst,
    Redst,
    imtex,
    retex,
    # lambda0,
    p0min,
    q0min,
    dp0,
    dq0,
    p0grid,
    k,
    m2,
    Ek,
    Temper,
    Nf,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    @inbounds if i <= length(Imdst)
        Imdst[i] = Twopointflow_Texture(
            p0grid[i],
            p0min,
            q0min,
            dp0,
            dq0,
            k,
            m2,
            Ek,
            Temper,
            imtex,
            # 0.0f0,
            Nf,
        )
        Redst[i] = Twopointflow_Texture(
            p0grid[i],
            p0min,
            q0min,
            dp0,
            dq0,
            k,
            m2,
            Ek,
            Temper,
            retex,
            # lambda0,
            Nf,
        )
    end
    return
end


function ImRewarp2!(
    Imdst,
    Redst,
    imtex,
    retex,
    p0min,
    q0min,
    dp0,
    dq0,
    p0grid,
    q0grid,
    k,
    m2,
    Ek,
    Temper,
    ϵ,
    Nf,
)
    kernel = @cuda launch = false ImRewarp(
        Imdst,
        Redst,
        imtex,
        retex,
        p0min,
        q0min,
        dp0,
        dq0,
        p0grid,
        q0grid,
        k,
        m2,
        Ek,
        Temper,
        ϵ,
        Nf,
    )
    config = launch_configuration(kernel.fun)
    threads = Base.min(length(Imdst), config.threads)
    blocks = cld(length(Imdst), threads)
    CUDA.@sync begin
        kernel(
            Imdst,
            Redst,
            imtex,
            retex,
            p0min,
            q0min,
            dp0,
            dq0,
            p0grid,
            q0grid,
            k,
            m2,
            Ek,
            Temper,
            ϵ,
            Nf,
            ;
            threads,
            blocks,
        )
    end
end









function ImReTwowarp2!(
    Imdst,
    Redst,
    imtex,
    retex,
    # lambda0,
    p0min,
    q0min,
    dp0,
    dq0,
    p0grid,
    k,
    m2,
    Ek,
    Temper,
    Nf;
    threads = 1024,
    blocks,
)
    @cuda threads = threads blocks = blocks ImReTwowarp(
        Imdst,
        Redst,
        imtex,
        retex,
        # lambda0,
        p0min,
        q0min,
        dp0,
        dq0,
        p0grid,
        k,
        m2,
        Ek,
        Temper,
        Nf,
    )
end



# function tchanelSolveFourPointGPU(
#     Temper::T;
#     kmin::T,
#     krang::RGscale = lpakrang,
#     dp0::T = 0.125f0,
#     dq0::T = 0.125f0,
#     ϵ::T = 1.0f0,
#     Nf::T = 4.0f0,
#     m2ini::T = Float32(-0.18216546875 * lpakrang.UV^2),
#     margin::Int = 0,
#     pgridmax::T = lpakrang.UV / 8.0f0,
#     config::ODEConfig,
# ) where {T<:Number}
#     GC.gc(true)
#     CUDA.reclaim()
#     p0min = -pgridmax - margin * dp0
#     p0max = pgridmax + margin * dp0
#     q0min = -pgridmax - margin * dq0
#     q0max = pgridmax + margin * dq0
#     p0gridcpu = collect(p0min:dp0:p0max)
#     q0gridcpu = collect(q0min:dq0:q0max)'
#     p0grid = p0gridcpu |> cu
#     q0grid = q0gridcpu |> cu
#     u0 = cu(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
#     # u0=sol_zero[end]
#     ImlambdaTextureMem = CuTextureArray(u0[:, 2:end, 1])
#     RelambdaTextureMem = CuTextureArray(u0[:, 2:end, 2])
#     idxzero = Int((length(p0grid) + 1) / 2)
#     threads = 1024
#     blocks = cld(length(p0grid), threads)
#     saved_values = SavedValues(Float32, Vector{Float32})
#     cb = SavingCallback(
#         (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
#         saved_values,
#     )
#     function flow_tgpu(du, u, p, k)
#         # k*0.5
#         ImGamm2 = @view u[:, 1, 1]
#         ReGamm2 = @view u[:, 1, 2]
#         Imlambda = @view u[:, 2:end, 1]
#         Relambda = @view u[:, 2:end, 2]
#         dkImGamm2 = @view du[:, 1, 1]
#         dkReGamm2 = @view du[:, 1, 2]
#         dkImlambda = @view du[:, 2:end, 1]
#         dkRelambda = @view du[:, 2:end, 2]

#         m2 = -ReGamm2[idxzero]
#         Relambda0 = Relambda[idxzero, idxzero]
#         Ek = Eb(k, m2)
#         idxEk = round(Int32, max(min((Ek - q0min) / dq0 + 1.0f0, length(q0grid)), 1))
#         idxmEk = round(Int32, max(min((-Ek - q0min) / dq0 + 1.0f0, length(q0grid)), 1))

#         copyto!(ImlambdaTextureMem, Imlambda)
#         copyto!(RelambdaTextureMem, Relambda)
#         # copyto!(RelambdaTextureMemB, RelambdaB)
#         Imlambdatex =
#             CuTexture(ImlambdaTextureMem; interpolation = CUDA.LinearInterpolation())
#         Relambdatex =
#             CuTexture(RelambdaTextureMem; interpolation = CUDA.LinearInterpolation())


#         println("k=", k)
#         # Relambdacpu=Array(Relambda)[:, idxEk]
#         # plot(-400.0:0.5:400.0, @. RelambdaAcpu^2/RelambdaBcpu^2)
#         # plot(p0min:dp0:p0max, @. RelambdaAcpu/(RelambdaBcpu))|>display
#         # display(plot!(-400.0:0.5:400.0, @. RelambdaAcpu^2/RelambdaBcpu^2,ylims =(-1,5000)))
#         # plot(p0min:dp0:p0max, Imlambdacpu1)
#         # display(plot!(p0min:dp0:p0max, Imlambdacpu2))
#         # Imlambdacpu1=Array(dkImlambda)[:, idxEk]
#         # Imlambdacpu2=Array(dkImlambda)[:, idxmEk]
#         # # plot(p0min:dp0:p0max, Imlambdacpu1)
#         # display(plot(p0min:dp0:p0max, Imlambdacpu1+Imlambdacpu2))
#         ImRewarp2!(
#             dkImlambda,
#             dkRelambda,
#             Imlambdatex,
#             Relambdatex,
#             p0min,
#             q0min,
#             dp0,
#             dq0,
#             p0grid,
#             q0grid,
#             k,
#             m2,
#             Ek,
#             Temper,
#             ϵ,
#             Nf,
#         )
#         ImReTwowarp2!(
#             dkImGamm2,
#             dkReGamm2,
#             Imlambdatex,
#             Relambdatex,
#             # Relambda0,
#             p0min,
#             q0min,
#             dp0,
#             dq0,
#             p0grid,
#             k,
#             m2,
#             Ek,
#             Temper,
#             Nf,
#             threads = threads,
#             blocks = blocks,
#         )

#         # Imlambdacpu = Array(Imlambda[1:2:end,1:2:end])'
#         # display(heatmap(Imlambdacpu))
#         # display(plot(0:0.125:200,Array(dkRelambda[100,1601:end])-(Array(dkRelambda[100,1:1601])|>reverse)))
#         if any(isnan, dkImlambda)
#             println("Imdk=", dkImlambda[end, 1], "Im=", Imlambda[end, 1])
#         end
#         # sleep(2.0)
#         # println("imdk[0,0]=",dkImlambda[1601,1601])

#         # broadcast!(dkImlambda, p0grid, q0grid) do p0, q0
#         #     Imlambdaflow(p0, q0, k, 0.0f0, Ek, Temper, Imlambdafun, Relambdafun, ϵ, Nf)
#         # end
#         # broadcast!(dkRelambda, p0grid, q0grid) do p0, q0
#         #     Relambdaflow(p0, q0, k, 0.0f0, Ek, Temper, Imlambdafun, Relambdafun, ϵ, Nf)
#         # end
#     end
#     alg = Tsit5()
#     prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
#     sol = solve(
#         prob,
#         alg,
#         # dt=0.5,
#         # dtmax = 0.5f0,
#         abstol = config.atol,
#         reltol = config.rtol,
#         calck = false,
#         # dt=0.0000001,
#         dtmin = config.dtmin,
#         # adaptive = config.adaptive,
#         dense = config.dense,
#         save_on = config.save_on,
#         save_start = config.save_start,
#         save_end = config.save_end,
#         callback = cb,
#         progress = config.progress,
#     )
#     img2 = Array(sol.u[end][idxzero:end, 1, 1])
#     reg2 = Array(sol.u[end][idxzero:end, 1, 2])
#     spec = @. img2 / (img2^2 + reg2^2)
#     return SpecSolution(
#         sol,
#         u0,
#         Array(u0[idxzero:end, 1, 1]),
#         Array(u0[idxzero:end, 1, 2]),
#         saved_values.t,
#         hcat(saved_values.saveval...)[1, :],
#         hcat(saved_values.saveval...)[2, :],
#         p0gridcpu[idxzero:end],
#         q0gridcpu[idxzero:end],
#         Array(sol.u[end][idxzero:end, 1, 1]),
#         Array(sol.u[end][idxzero:end, 1, 2]),
#         spec,
#     )
# end



#Log grid

function tchanelSolveFourPointGPU(
    Temper::T;
    kmin::T,
    krang::RGscale = lpakrang,
    dp0::T = 0.125,
    dq0::T = 0.125,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    margin::Int = 0,
    pgridmax::T = lpakrang.UV / 8.0,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -pgridmax - margin * dp0
    p0max = pgridmax + margin * dp0
    q0min = -pgridmax - margin * dq0
    q0max = pgridmax + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
    lengthy = length(p0gridcpu)
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = Int((length(p0grid) + 1) / 2)
    threads = 1024
    blocks = cld(length(p0grid), threads)
    saved_values = SavedValues(Float32, Vector{Float32})
    cb = SavingCallback(
        (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
        saved_values,
    )
    function flow_tgpu(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]

        m2 = -ReGamm2[idxzero]
        Relambda0 = Relambda[idxzero, idxzero]
        Ek = Eb(k, m2)
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            Relambda,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync ImRelambdaflow_Launchxy(
            dkImlambda,
            dkRelambda,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = (16,16),
        )
        println("k=", k)
        CUDA.@sync ImReTwopointflow_Launch(
            dkImGamm2,
            dkReGamm2,
            lambdagrid,
            k,
            m2,
            Temper,
            Nf,
            threads = threads,
        )
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
    sol = solve(
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
        callback = cb,
        progress = config.progress,
    )
    img2 = Array(sol.u[end][idxzero:end, 1, 1])
    reg2 = Array(sol.u[end][idxzero:end, 1, 2])
    spec = @. img2 / (img2^2 + reg2^2)
    return SpecSolution(
        Array(sol.u[end]),
        Array(u0),
        Array(u0[idxzero:end, 1, 1]),
        Array(u0[idxzero:end, 1, 2]),
        saved_values.t,
        hcat(saved_values.saveval...)[1, :],
        hcat(saved_values.saveval...)[2, :],
        p0gridcpu[idxzero:end],
        q0gridcpu[idxzero:end],
        Array(sol.u[end][idxzero:end, 1, 1]),
        Array(sol.u[end][idxzero:end, 1, 2]),
        spec,
    )
end



function tchanelSolveFourPointGPU_log(
    Temper::T;
    kmin::T,
    krang::RGscale = lpakrang,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    gridlength::Int = 1000,
    pgridmax::T = lpakrang.UV / 8.0,
    pgridmin::T = 1e-3,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()

    p0min = pgridmin
    p0max = pgridmax
    q0min = pgridmin
    q0max = pgridmax
    p0min_log = log(pgridmin)
    p0max_log = log(pgridmax)
    q0min_log = log(pgridmin)
    q0max_log = log(pgridmax)
    p0gridcpu =
        T.(
            vcat(
                -exp.(range(p0max_log, p0min_log, gridlength)),
                0,
                exp.(range(p0min_log, p0max_log, gridlength)),
            ),
        )
    q0gridcpu = T.(
            vcat(
                -exp.(range(q0max_log, q0min_log, gridlength)),
                0,
                exp.(range(q0min_log, q0max_log, gridlength)),
            ),
        )
    dq0_log=step(range(q0min_log, q0max_log, gridlength))

    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
    lengthy = 2*gridlength+1
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = gridlength+1
    threads = 1024
    blocks = cld(length(p0grid), threads)
    saved_values = SavedValues(T, Vector{T})
    cb = SavingCallback(
        (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
        saved_values,
    )
    function flow_tgpu(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]

        m2 = -ReGamm2[idxzero]
        Relambda0 = Relambda[idxzero, idxzero]
        Ek = Eb(k, m2)
        println("k=", k)
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            Relambda,
            lengthy = lengthy,
            dq0 = dq0_log,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
            Interpfun=Interpolate_Log_Linear,
        )
        CUDA.@sync ImRelambdaflow_Launchxy(
            dkImlambda,
            dkRelambda,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = (16,16),
        )
       CUDA.@sync ImReTwopointflow_Launch(
            dkImGamm2,
            dkReGamm2,
            lambdagrid,
            k,
            m2,
            Temper,
            Nf,
            threads = threads,
        )
    end
    alg = VCABM()
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
    sol = solve(
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
        callback = cb,
        progress = config.progress,
    )
    img2 = Array(sol.u[end][idxzero:end, 1, 1])
    reg2 = Array(sol.u[end][idxzero:end, 1, 2])
    spec = @. img2 / (img2^2 + reg2^2)
    return SpecSolution(
        Array(sol.u[end]),
        Array(u0),
        Array(u0[idxzero:end, 1, 1]),
        Array(u0[idxzero:end, 1, 2]),
        saved_values.t,
        hcat(saved_values.saveval...)[1, :],
        hcat(saved_values.saveval...)[2, :],
        p0gridcpu[idxzero:end],
        q0gridcpu[idxzero:end],
        Array(sol.u[end][idxzero:end, 1, 1]),
        Array(sol.u[end][idxzero:end, 1, 2]),
        spec,
    )
end






function tchanelSolveFourPointGPU_ini(
    Temper::T;
    maxiters=1,
    dp0::T = 0.125,
    dq0::T = 0.125,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    margin::Int = 0,
    pgridmax::T = lpakrang.UV / 8.0,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -pgridmax - margin * dp0
    p0max = pgridmax + margin * dp0
    q0min = -pgridmax - margin * dq0
    q0max = pgridmax + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
    lengthy = length(p0gridcpu)
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = Int((length(p0grid) + 1) / 2)
    threads = 1024
    blocks = cld(length(p0grid), threads)
    function flow_tgpu(du, u, k)
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]

        m2 = m2ini
        Ek = Eb(k, m2)
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            Relambda,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync ImRelambdaflow_Launch(
            dkImlambda,
            dkRelambda,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = 500,
        )
        Relambda0 = Relambda[idxzero, idxzero]
        dkRelambda .=T(RelambdaUV/Relambda0)*dkRelambda
        u.=du
        println("k=", k)
    end
    du0=similar(u0);
    for i in 1:Int32(maxiters)
        flow_tgpu(du0, u0, T(lpakrang.UV))
    end
    return u0
end


function tchanelSolveFourPointGPU_invRe(
    Temper::T;
    kmin::T,
    krang::RGscale = lpakrang,
    dp0::T = 0.125,
    dq0::T = 0.125,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    margin::Int = 0,
    pgridmax::T = lpakrang.UV / 8.0,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -pgridmax - margin * dp0
    p0max = pgridmax + margin * dp0
    q0min = -pgridmax - margin * dq0
    q0max = pgridmax + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d_InvRe(p0gridcpu, q0gridcpu, m2ini))
    lengthy = length(p0gridcpu)
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = Int((length(p0grid) + 1) / 2)
    threads = 1024
    blocks = cld(length(p0grid), threads)
    saved_values = SavedValues(Float32, Vector{Float32})
    cb = SavingCallback(
        (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
        saved_values,
    )
    function flow_tgpu(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        invRelambda = @view u[:, 2:end, 2]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        invdkRelambda = @view du[:, 2:end, 2]

        m2 = -ReGamm2[idxzero]
        # Relambda0 = Relambda[idxzero, idxzero]
        Ek = Eb(k, m2)
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            invRelambda,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync ImRelambdaflow_Launch(
            dkImlambda,
            invdkRelambda,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = 500,
            Imlambdaflow_fun = Imlambdaflow_GPU_InvRe,
            Relambdaflow_fun = Relambdaflow_GPU_InvRe,
        )
        # invdkRelambda .= -invdkRelambda.*invRelambda.^2
        println("k=", k)
        CUDA.@sync ImReTwopointflow_Launch(
            dkImGamm2,
            dkReGamm2,
            lambdagrid,
            k,
            m2,
            Temper,
            Nf,
            threads = threads,
            ImTwopointflow_fun = ImTwopointflow_GPU,
            ReTwopointflow_fun = ReTwopointflow_GPU_InvRe,
        )
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
    sol = solve(
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
        callback = cb,
        progress = config.progress,
    )
    img2 = Array(sol.u[end][idxzero:end, 1, 1])
    reg2 = Array(sol.u[end][idxzero:end, 1, 2])
    spec = @. img2 / (img2^2 + reg2^2)
    return SpecSolution(
        sol,
        u0,
        Array(u0[idxzero:end, 1, 1]),
        Array(u0[idxzero:end, 1, 2]),
        saved_values.t,
        hcat(saved_values.saveval...)[1, :],
        hcat(saved_values.saveval...)[2, :],
        p0gridcpu[idxzero:end],
        q0gridcpu[idxzero:end],
        Array(sol.u[end][idxzero:end, 1, 1]),
        Array(sol.u[end][idxzero:end, 1, 2]),
        spec,
    )
end


function tchanelSolveFourPointGPU_AB(
    Temper::T;
    kmin::T,
    krang::RGscale = lpakrang,
    dp0::T = 0.125,
    dq0::T = 0.125,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    margin::Int = 0,
    pgridmax::T = lpakrang.UV / 8.0,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -pgridmax - margin * dp0
    p0max = pgridmax + margin * dp0
    q0min = -pgridmax - margin * dq0
    q0max = pgridmax + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d_AB(p0gridcpu, q0gridcpu, m2ini))
    lengthy = length(p0gridcpu)
    lambdagridA = LambdaGridini(T, lengthy)
    lambdagridB = LambdaGridini(T, lengthy)
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = Int((length(p0grid) + 1) / 2)
    threads = 1024
    blocks = cld(length(p0grid), threads)
    saved_values = SavedValues(Float32, Vector{Float32})
    cb = SavingCallback(
        (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
        saved_values,
    )
    function flow_tgpu(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        RelambdaA = @view u[:, 2:end, 2]
        RelambdaB = @view u[:, 2:end, 3]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambdaA = @view du[:, 2:end, 2]
        dkRelambdaB = @view du[:, 2:end, 3]
        Relambda=RelambdaA+RelambdaB
        # m2 = -ReGamm2[idxzero]
        m2=0.0
        # Relambda0 = Relambda[idxzero, idxzero]
        Ek = Eb(k, m2)
        CUDA.@sync reverse2!(
            lambdagridA,
            Imlambda,
            RelambdaA,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync reverse2!(
            lambdagridB,
            Imlambda,
            RelambdaB,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            Relambda,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        #Aflow
        CUDA.@sync ImRelambdaflow_Launch(
            dkImlambda,
            dkRelambdaA,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = 500,
            Imlambdaflow_fun = Imlambdaflow_GPU,
            Relambdaflow_fun = Relambdaflow_GPU_A,
        )
        #Bflow
        CUDA.@sync ImRelambdaflow_Launch(
            dkImlambda,
            dkRelambdaB,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            4*ϵ,
            Nf,
            threads = 500,
            Imlambdaflow_fun = Imlambdaflow_GPU,
            Relambdaflow_fun = Relambdaflow_GPU_B,
        )
        println("k=", k)
        CUDA.@sync ImReTwopointflow_Launch(
            dkImGamm2,
            dkReGamm2,
            lambdagridB,
            k,
            m2,
            Temper,
            Nf,
            threads = threads,
        )
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
    sol = solve(
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
        callback = cb,
        progress = config.progress,
    )
    img2 = Array(sol.u[end][idxzero:end, 1, 1])
    reg2 = Array(sol.u[end][idxzero:end, 1, 2])
    spec = @. img2 / (img2^2 + reg2^2)
    return SpecSolution(
        sol,
        u0,
        Array(u0[idxzero:end, 1, 1]),
        Array(u0[idxzero:end, 1, 2]),
        saved_values.t,
        hcat(saved_values.saveval...)[1, :],
        hcat(saved_values.saveval...)[2, :],
        p0gridcpu[idxzero:end],
        q0gridcpu[idxzero:end],
        Array(sol.u[end][idxzero:end, 1, 1]),
        Array(sol.u[end][idxzero:end, 1, 2]),
        spec,
    )
end








function tchanelSolveFourPointGPU2(
    Temper::T;
    kmin::T,
    krang::RGscale = lpakrang,
    dp0::T = 0.125,
    dq0::T = 0.125,
    ϵ::T = 1.0,
    Nf::T = 4.0,
    m2ini::T = T(-0.18216546875 * lpakrang.UV^2),
    margin::Int = 0,
    pgridmax::T = lpakrang.UV / 8.0,
    config::ODEConfig,
) where {T<:Number}
    GC.gc(true)
    CUDA.reclaim()
    p0min = -pgridmax - margin * dp0
    p0max = pgridmax + margin * dp0
    q0min = -pgridmax - margin * dq0
    q0max = pgridmax + margin * dq0
    p0gridcpu = collect(p0min:dp0:p0max)
    q0gridcpu = collect(q0min:dq0:q0max)'
    p0grid = p0gridcpu |> cu
    q0grid = q0gridcpu |> cu
    u0 = CuArray(init_fourpoint_2d(p0gridcpu, q0gridcpu, m2ini))
    lengthy = length(p0gridcpu)
    lambdagrid = LambdaGridini(T, lengthy)
    # u0=sol_zero[end]
    idxzero = Int((length(p0grid) + 1) / 2)
    threads = 1024
    blocks = cld(length(p0grid), threads)
    saved_values = SavedValues(Float32, Vector{Float32})
    cb = SavingCallback(
        (u, t, integrator) -> [-u[idxzero, 1, 2], u[idxzero, idxzero+1, 2]],
        saved_values,
    )
    function flow_tgpu(du, u, p, k)
        ImGamm2 = @view u[:, 1, 1]
        ReGamm2 = @view u[:, 1, 2]
        Imlambda = @view u[:, 2:end, 1]
        Relambda = @view u[:, 2:end, 2]
        dkImGamm2 = @view du[:, 1, 1]
        dkReGamm2 = @view du[:, 1, 2]
        dkImlambda = @view du[:, 2:end, 1]
        dkRelambda = @view du[:, 2:end, 2]

        m2 = -ReGamm2[idxzero]
        Relambda0 = Relambda[idxzero, idxzero]
        Ek = Eb(k, m2)
        CUDA.@sync reverse2!(
            lambdagrid,
            Imlambda,
            Relambda,
            lengthy = lengthy,
            dq0 = dq0,
            q0min = q0min,
            Ek = Ek,
            threads = threads,
        )
        CUDA.@sync ImRelambdaflow_Launchxy(
            dkImlambda,
            dkRelambda,
            p0grid,
            q0grid,
            lambdagrid,
            k,
            m2,
            Temper,
            ϵ,
            Nf,
            threads = (16,16),
        )
        println("k=", k)
        CUDA.@sync ImReTwopointflow_Launch(
            dkImGamm2,
            dkReGamm2,
            lambdagrid,
            k,
            m2,
            Temper,
            Nf,
            threads = threads,
        )
    end
    alg = Tsit5()
    prob = ODEProblem(flow_tgpu, u0, (krang.UV, kmin))
    sol = solve(
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
        callback = cb,
        progress = config.progress,
    )
    return
    # img2 = Array(sol.u[end][idxzero:end, 1, 1])
    # reg2 = Array(sol.u[end][idxzero:end, 1, 2])
    # spec = @. img2 / (img2^2 + reg2^2)
    # return SpecSolution(
    #     Array(sol.u[end]),
    #     Array(u0),
    #     Array(u0[idxzero:end, 1, 1]),
    #     Array(u0[idxzero:end, 1, 2]),
    #     saved_values.t,
    #     hcat(saved_values.saveval...)[1, :],
    #     hcat(saved_values.saveval...)[2, :],
    #     p0gridcpu[idxzero:end],
    #     q0gridcpu[idxzero:end],
    #     Array(sol.u[end][idxzero:end, 1, 1]),
    #     Array(sol.u[end][idxzero:end, 1, 2]),
    #     spec,
    # )
end
