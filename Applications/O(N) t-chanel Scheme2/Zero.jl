# function flowslefZero(du, u, p, k)
#     # global testi+=1
#     m2 = -u[1]
#     lambdaq0 = @view u[2:end]
#     lambdaq0dk = @view du[2:end] # za
#     Ek = Epi(k, m2)
#     idx = max(round(Int64, 10* Ek),2)
#     lambdaq0fun = Spline1D(q0grid[idx:idx+4], lambdaq0[idx:idx+4])
#     # lambdaq0fun = CubicSplineInterpolation(q0grid[idx]:q0grid[idx+4], lambdaq0[idx:idx+4])
#     # itp = interpolate(lambdaq0[idx:idx+4], BSpline(Cubic(Line(OnGrid()))))
#     # lambdaq0fun = scale(itp, q0grid[idx:idx+4])
#     # lambdaq0fun = interpolate((q0grid[idx:idx+4],), lambdaq0[idx:idx+4], Gridded(Linear()))
#     lambdaEk = lambdaq0fun(Ek)
#     du[1] =
#         twopointflow3(0.0, k, m2, p, lambdaEk, 4) +
#         twopointflow3(0.0, k, m2, p, u[2], 2)
#     # @time lambdaq0dk .= Relambdaflow3.(0.0, q0grid, k, m2, u[2], p, 0.0, lambdaq0,0.1)
#     # ThreadsX.foreach(
#     #     referenceable(lambdaq0dk),
#     #     q0grid,
#     #     lambdaq0,
#     #     basesize = 25,
#     # ) do z, x, y
#     #     z[] = Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5)
#     # end
#     # vmapntt!((x,y)->Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5), lambdaq0dk, q0grid, lambdaq0)
#     # setvalue(lambdaq0dk, v1,v2) = @tullio lambdaq0dk[i]=Relambdaflow3(0.0, v1[i], k, m2, u[2], p, 0.0, v2[i], ϵ5)
#     # setvalue(lambdaq0dk, q0grid,lambdaq0,k,m2,u[2],p)
#     @tullio threads=168 avx=4 grad=false lambdaq0dk[i]=Relambdaflow3(0.0, q0grid[i], k, m2, u[2], p, 0.0, lambdaq0[i], ϵ5)
#       # println("du1=",du[2]," ","duend=",du[end])
# end

# probself = ODEProblem(
#     flowslefZero,
#     collect(vcat([116585.9, fill(6 * 8.0, length(q0grid))]...)),
#     (800.0,1e-18),
#     # (800.0,400.0),
#     Tself+0.001939,
# )
# 0.00728
# lines(sol.t,hcat(-sol.u...)[1,:])
# hcat(sol1.u...)
# lines!(sol1.t,hcat(sol1.u...)[1,:])
#
# @time sol.t
#
# lines(q0grid,sol.u[end][2:end])
#
# lines!(q0grid,sol.u[1][2:end])
#
# current_figure()

# function flowslefZeroNf0(du, u, p, k)
#     # global testi+=1
#     m2 = -u[1]
#     lambdaq0 = @view u[2:end]
#     lambdaq0dk = @view du[2:end] # za
#     Ek = Epi(k, m2)
#     idx = max(round(Int64, 10* Ek),2)
#     lambdaq0fun = Spline1D(q0grid[idx:idx+4], lambdaq0[idx:idx+4])
#     # lambdaq0fun = CubicSplineInterpolation(q0grid[idx]:q0grid[idx+4], lambdaq0[idx:idx+4])
#     # itp = interpolate(lambdaq0[idx:idx+4], BSpline(Cubic(Line(OnGrid()))))
#     # lambdaq0fun = scale(itp, q0grid[idx:idx+4])
#     # lambdaq0fun = interpolate((q0grid[idx:idx+4],), lambdaq0[idx:idx+4], Gridded(Linear()))
#     lambdaEk = lambdaq0fun(Ek)
#     du[1] = twopointflow3(0.0, k, m2, p, u[2], 2)
#     # @time lambdaq0dk .= Relambdaflow3.(0.0, q0grid, k, m2, u[2], p, 0.0, lambdaq0,0.1)
#     # ThreadsX.foreach(
#     #     referenceable(lambdaq0dk),
#     #     q0grid,
#     #     lambdaq0,
#     #     basesize = 25,
#     # ) do z, x, y
#     #     z[] = Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5)
#     # end
#     # vmapntt!((x,y)->Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5), lambdaq0dk, q0grid, lambdaq0)
#     # setvalue(lambdaq0dk, v1,v2) = @tullio lambdaq0dk[i]=Relambdaflow3(0.0, v1[i], k, m2, u[2], p, 0.0, v2[i], ϵ5)
#     # setvalue(lambdaq0dk, q0grid,lambdaq0,k,m2,u[2],p)
#     @tullio threads=168 avx=4 grad=false lambdaq0dk[i]=Relambdaflow3(0.0, q0grid[i], k, m2, u[2], p, 0.0, lambdaq0[i], ϵ5,0.0)
#       # println("du1=",du[2]," ","duend=",du[end])
# end
# # 0.35093892265
#
#
#
# probselfNf0 = ODEProblem(
#     flowslefZeroNf0,
#     collect(vcat([0.3509389062993*116585.9, fill(6 * 8.0, length(q0grid))]...)),
#     (800.0,0.1),
#     # (800.0,400.0),
#     Tself+0.001939,
# )


function flowslefZeroNf(du, u, p, k)
    # global testi+=1
    temper, Nf, batch, spacing=p
    m2 = -u[1]
    lambdaq0 = @view u[2:end]
    lambdaq0dk = @view du[2:end] # za
    Ek = Epi(k, m2)
    idx = max(round(Int64, spacing* Ek),2)
    lambdaq0fun = Spline1D(q0grid[idx:idx+4], lambdaq0[idx:idx+4])
    # lambdaq0fun = CubicSplineInterpolation(q0grid[idx]:q0grid[idx+4], lambdaq0[idx:idx+4])
    # itp = interpolate(lambdaq0[idx:idx+4], BSpline(Cubic(Line(OnGrid()))))
    # lambdaq0fun = scale(itp, q0grid[idx:idx+4])
    # lambdaq0fun = interpolate((q0grid[idx:idx+4],), lambdaq0[idx:idx+4], Gridded(Linear()))
    lambdaEk = lambdaq0fun(Ek)
    du[1] =
        twopointflow3(0.0, k, m2, temper, lambdaEk, Nf) +
        twopointflow3(0.0, k, m2, temper, u[2], 2)
    # @time lambdaq0dk .= Relambdaflow3.(0.0, q0grid, k, m2, u[2], p, 0.0, lambdaq0,0.1)
    # ThreadsX.foreach(
    #     referenceable(lambdaq0dk),
    #     q0grid,
    #     lambdaq0,
    #     basesize = 25,
    # ) do z, x, y
    #     z[] = Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5)
    # end
    # vmapntt!((x,y)->Relambdaflow3(0.0, x, k, m2, u[2], p, 0.0, y, ϵ5), lambdaq0dk, q0grid, lambdaq0)
    # setvalue(lambdaq0dk, v1,v2) = @tullio lambdaq0dk[i]=Relambdaflow3(0.0, v1[i], k, m2, u[2], p, 0.0, v2[i], ϵ5)
    # setvalue(lambdaq0dk, q0grid,lambdaq0,k,m2,u[2],p)
    @tullio threads=batch avx=4 grad=false lambdaq0dk[i]=Relambdaflow3(0.0, q0grid[i], k, m2, u[2], temper, 0.0, lambdaq0[i], ϵ5,Nf)
      # println("du1=",du[2]," ","duend=",du[end])
end

# 0.03390691

# probselfNf = ODEProblem(
#     flowslefZeroNf,
#     collect(vcat([116585.9*sqrt(16.5/6), fill(6 * 8.0*sqrt(12/16.5), length(q0grid))]...)),
#     (800.0,0.1),
#     # (800.0,400.0),
#     [145.0,10.0],
# )
# Tself+0.001939-120.0

function tchanelZeroSelfSolve(
    Temper,
    krang::RGscale,
    flowfun;
    u0,
    Nf,
    config::ODEConfig,
)
    probselfNf1 = ODEProblem(
        flowfun,
        collect(vcat([-u0[1], fill(u0[2], length(q0grid))]...)),
        krang.kspan,
        # (800.0,400.0),
        (
            Temper,
            Nf,
            ceil(Int, length(q0grid) / Threads.nthreads()),#batch
            1 / (q0grid[2] - q0grid[1]), #spacing,
        ),
    )
    solselfNf1 = solve(
        probselfNf1,
        config.alg,
        dense = config.dense,
        save_on = config.save_on,
        save_start = config.save_start,
        save_timeseries = false,
        save_end = config.save_end,
        reltol = config.rtol,
        abstol = config.atol,
        progress = config.progress,
        dtmax = config.dtmax,
        dtmin = config.dtmin,
    )
    -solselfNf1.u[end][1][1]
end


function tchanelZeroSelfSolve2(
    Temper,
    krang::RGscale,
    flowfun;
    u0,
    Nf,
    config::ODEConfig,
)
    probselfNf1 = ODEProblem(
        flowfun,
        collect(vcat([-u0[1], fill(u0[2], length(q0grid))]...)),
        krang.kspan,
        # (800.0,400.0),
        (
            Temper,
            Nf,
            ceil(Int, length(q0grid) / Threads.nthreads()),#batch
            1 / (q0grid[2] - q0grid[1]), #spacing,
        ),
    )
    solselfNf1 = solve(
        probselfNf1,
        config.alg,
        # dense = config.dense,
        # save_on = config.save_on,
        # save_start = config.save_start,
        # save_timeseries = false,
        # save_end = config.save_end,
        reltol = config.rtol,
        abstol = config.atol,
        progress = config.progress,
        dtmax = config.dtmax,
        dtmin = config.dtmin,
    )
end



# function tchanelZeroSelfSolve2(Temper, krang::RGscale, flowfun, u0, Nf)
#     probselfNf1 = ODEProblem(
#         flowfun,
#         collect(vcat([-u0[1], fill(u0[2], length(q0grid))]...)),
#         krang.kspan,
#         # (800.0,400.0),
#         [Temper, Nf],
#     )
#     solselfNf1 = solve(
#         probselfNf1,
#         Tsit5(),
#         # dense = false,
#         # save_on = false,
#         # save_start = false,
#         # save_timeseries=false,
#         # save_end = true,
#         reltol = 1e-6,
#         abstol = 1e-6,
#         progress = true,
#         dtmax=50.0,
#         dtmin=1e-14,
#     )
# end
