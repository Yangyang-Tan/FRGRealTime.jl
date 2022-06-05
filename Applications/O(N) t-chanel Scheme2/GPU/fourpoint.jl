function getsolN1(T1)
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
end

function getsolN10(T1)
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
end

function getsolN2(T1)
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
end

function getsolN4(T1)
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
end

specsolN1=getsolN1(minmassN1)
specsolN10=getsolN10(minmassN10)

device!(0)

specsolN2=getsolN2(minmassN2)
specsolN4=getsolN4(minmassN4)

specsolN4_0p01=getsolN4(minmassN4+0.01)
specsolN4_0p001=getsolN4(minmassN4+0.001)
specsolN4_0p0001=getsolN4(minmassN4+0.0001)
specsolN4_0p00005=getsolN4(minmassN4+0.00005)

plot(specsolN4_0p01.p0[2:end].|>log,specsolN4_0p01.spec[2:end].|>log)
plot!(specsolN4_0p001.p0[2:end].|>log,specsolN4_0p001.spec[2:end].|>log)

plot!(specsolN4_0p0001.p0[2:end].|>log,specsolN4_0p0001.spec[2:end].|>log)
plot!(specsolN4_0p00005.p0[2:end].|>log,specsolN4_0p00005.spec[2:end].|>log)


using DataFrames
using OrderedCollections, FileIO
dictN1 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN1), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 0.3591,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 1,
    "Temper" => minmassN1,
),);
dictN10 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN10), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 2.4148,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 10,
    "Temper" => minmassN10,
),);

dictN2 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN2), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 0.5620,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 2,
    "Temper" => minmassN2,
),);
dictN4 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN4), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 1,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 4,
    "Temper" => minmassN4,
),);

dictN4_0p01 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN4_0p01), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 1,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 4,
    "Temper" => minmassN4+0.01,
),);
dictN4_0p001 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN4_0p001), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 1,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 4,
    "Temper" => minmassN4+0.001,
),);
dictN4_0p0001 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN4_0p0001), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 1,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 4,
    "Temper" => minmassN4+0.0001,
),);
dictN4_0p00005 = OrderedDict((
    (
        string.(fieldnames(SpecSolution)) .=>
            getfield.(Ref(specsolN4_0p00005), fieldnames(SpecSolution))
    )...,
    "m2ini" => -0.18216546875 * lpakrang.UV^2,
    "λini"=>24 / 1,
    "ϵ" => 2,
    "pgridmax" => 20.0,
    "pgridmin" => 1e-3,
    "Nf" => 4,
    "Temper" => minmassN4+0.00005,
),);
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=1.h5", dictN1)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=10.h5", dictN10)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=2.h5", dictN2)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=4.h5", dictN4)

save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=4_0p01.h5", dictN4_0p01)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=4_0p001.h5", dictN4_0p001)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=4_0p0001.h5", dictN4_0p0001)
save("/home/tyy/Documents/CTP-fRG-Test/Real-Time-data/Critical/N=4_0p00005.h5", dictN4_0p00005)
