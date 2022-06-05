using Roots
function getmassKernel(T1, T2)
    local m1, m2
    @sync begin
        @async begin
            device!(0)
            tempsol1 = tchanelSolveFourPointGPU_log(
                # 9.4685f0,
                # 128.22675,
                T1,
                m2ini = Float64(-0.18216546875 * 0.5 * lpakrang.UV^2),
                kmin = 1.0,
                ϵ = 2.0,
                pgridmax = 20.0,
                pgridmin = 1e-3,
                gridlength = 1000,
                Nf = 1.0,
                config = config_spec,
            )
            m1 = tempsol1.m2[end]
        end
        @async begin
            device!(1)
            tempsol2 = tchanelSolveFourPointGPU_log(
                # 9.4685f0,
                # 128.22675,
                T2,
                m2ini = Float64(-0.18216546875 * 0.5 * lpakrang.UV^2),
                kmin = 1.0,
                ϵ = 2.0,
                pgridmax = 20.0,
                pgridmin = 1e-3,
                gridlength = 1000,
                Nf = 1.0,
                config = config_spec,
            )
            m2 = tempsol2.m2[end]
        end
    end
    return m1, m2
end


function getmassKernelN1(T1)
    tempsol1 = tchanelSolveFourPointGPU_log(
        # 9.4685f0,
        # 128.22675,
        T1,
        λini = 24 / 0.3591,
        kmin = 1.0,
        ϵ = 2.0,
        pgridmax = 20.0,
        pgridmin = 1e-3,
        gridlength = 1000,
        Nf = 1.0,
        config = config_spec,
    )
    m1 = tempsol1.m2[end]
end

function getmassKernelN10(T1)
    tempsol1 = tchanelSolveFourPointGPU_log(
        # 9.4685f0,
        # 128.22675,
        T1,
        kmin = 1.0,
        ϵ = 2.0,
        pgridmax = 20.0,
        pgridmin = 1e-3,
        gridlength = 1000,
        Nf = 10.0,
        λini = 24 / 2.4148,
        config = config_spec,
    )
    m1 = tempsol1.m2[end]
end


function getmassKernelN2(T1)
    tempsol1 = tchanelSolveFourPointGPU_log(
        # 9.4685f0,
        # 128.22675,
        T1,
        λini = 24 / 0.5620,
        kmin = 1.0,
        ϵ = 2.0,
        pgridmax = 20.0,
        pgridmin = 1e-3,
        gridlength = 1000,
        Nf = 2.0,
        config = config_spec,
    )
    m1 = tempsol1.m2[end]
end

function getmassKernelN4(T1)
    tempsol1 = tchanelSolveFourPointGPU_log(
        # 9.4685f0,
        # 128.22675,
        T1,
        kmin = 1.0,
        ϵ = 2.0,
        pgridmax = 20.0,
        pgridmin = 1e-3,
        gridlength = 1000,
        Nf = 4.0,
        λini = 24 / 1,
        config = config_spec,
    )
    m1 = tempsol1.m2[end]
end


config_spec.atol = 1e-10
config_spec.rtol = 1e-10
TciniN1=116.35
TciniN10=135.22
getmassKernelN10(TciniN1)
getmassKernelN1(TciniN10)
TciniN2=122.17
TciniN4=128.228
getmassKernelN2(122.17)
getmassKernelN4(128.228)


using Roots
function getminmass(fun1,fun2,T1, T2)
    local mass1, mass2
    @sync begin
        @async begin
            device!(0)
            mass1=find_zero(fun1, (T1, T1 + 1), Bisection(), atol = 1e-7, rtol = 1e-7)
        end
        @async begin
            device!(1)
            mass2=find_zero(fun2, (T2, T2 + 1), Bisection(), atol = 1e-7, rtol = 1e-7)
        end
    end
    return mass1, mass2
end

config_spec.atol = 1e-10
config_spec.rtol = 1e-10
minmassN1,minmassN10 = getminmass(getmassKernelN1,getmassKernelN10,TciniN1, TciniN10)
minmassN2,minmassN4 = getminmass(getmassKernelN2,getmassKernelN4,TciniN2, TciniN4)
