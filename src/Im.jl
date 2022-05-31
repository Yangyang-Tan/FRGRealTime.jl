
function deltasumkAll(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    if p0>=0.0
        deltasumkfix(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    elseif p0<0.0
        -deltasumkfix(-p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    end
end


function deltasumkAllSmooth(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    if p0>=0.0
        deltasumkfixSmooth(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    elseif p0<0.0
        -deltasumkfixSmooth(-p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    end
end


function deltasumkAll_dp0(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    if p0>=0.0
        deltasumkfix_dp0(p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    elseif p0<0.0
        deltasumkfix_dp0(-p0, ps, k, T,Npi,IRScale, UVScale, mfun, lamfun)
    end
end



@doc raw"""
    F1All(p0, ps, k, m, T)

 $F_1$ is an odd function about $p_0$

when `p0>0`,`F1All=loopfunpp`

when `p0<0`,`dkF1All=-loopfunpp(-p0)`

`F1All` doesn't contains the type-2 delta function
"""
function F1All(p0, ps, k, m, T)
    if p0>=0.0
        loopfunpp(p0, ps, k, m, T)
    elseif p0<0.0
        -loopfunpp(-p0, ps, k, m, T)
    end
end


@doc raw"""
    F2All(p0, ps, k, m, T)

 $F_2$ is an odd function about $p_0$

when `p0>0`,`F2All=loopfunpm`

when `p0<0`,`F2All=-loopfunpm(-p0)`
"""
function F2All(p0, ps, k, m, T)
    if p0>=0.0
        loopfunpm(p0, ps, k, m, T)
    elseif p0<0.0
        -loopfunpm(-p0, ps, k, m, T)
    end
end







@doc raw"""
    dkF1All(p0, ps, k, m, T)

 $F_1$ is an odd function about $p_0$

when `p0>0`,`dkF1All=flowpp`

when `p0<0`,`dkF1All=-flowpp(-p0)`

`dkF1All` doesn't contains the type-2 delta function
"""
function dkF1All(p0, ps, k, m, T,δ=0.02)
    if p0>=0.0
        flowpp(p0, ps, k, m, T,δ)
    elseif p0<0.0
        -flowpp(-p0, ps, k, m, T,δ)
    end
end


function dkF1All(k, m, T)
    return 0.0
end


@doc raw"""
    dkF2All(p0, ps, k, m, T)

 $F_2$ is an odd function about $p_0$

when `p0>0`,`dkF2All=flowpm`

when `p0<0`,`dkF2All=-flowpm(-p0)`
"""
function dkF2All(p0, ps, k, m, T)
    if p0>=0.0
        flowpm(p0, ps, k, m, T)
    elseif p0<0.0
        -flowpm(-p0, ps, k, m, T)
    end
end

function dkF2All(k, m, T)
    if m > 0.0
        return 0.0
    elseif m <= 0.0
        Ek = Epi(k, m)
        return k / (16 * Ek * pi^2 * T - 16 * Ek * pi^2 * T * cosh(Ek / T))
    end
end



@doc raw"""
    dkF1Allintqs(p0, ps, qsmax, k, m, T)

 $F_1$ is an odd function about $p_0$

when `p0>0`,`dkF1Allintqs=flowpp_intcostheqs`

when `p0<0`,`dkF1Allintqs=-flowpp_intcostheqs(-p0)`

`dkF1Allintqs` doesn't contains the type-1 delta function

`dkF1Allintqs` doesn't contains the type-2 delta function
"""
function dkF1Allintqs(p0, ps, qsmax, k, m, T)
    sign(p0)*flowpp_intcostheqs(abs(p0), ps, qsmax, k, m, T)
end

function dkF1Allintqs_gpu(p0, ps, qsmax, k, m, T)
    sign(p0)*flowpp_intcostheqs_gpu(abs(p0), ps, qsmax, k, m, T)
end


@doc raw"""
    flowpm(p0, ps, k, m, T)

 $F_2$ is an odd function about $p_0$

when `p0>0`,`dkF2Allintqs=flowpm_intcostheqs`

when `p0<0`,`dkF2Allintqs=-flowpm_intcostheqs(-p0)`
"""
function dkF2Allintqs(p0, ps, qsmax, k, m, T)
    sign(p0)*flowpm_intcostheqs(abs(p0), ps, qsmax, k, m, T)
end


function dkF2Allintqs_gpu(p0, ps, qsmax, k, m, T)
    sign(p0)*flowpm_intcostheqs_gpu(abs(p0), ps, qsmax, k, m, T)
end

# function dkF2Allintqs(p0, ps, qsmax, k, m, T)
#     if p0>=0.0
#         0.0
#     elseif p0<0.0
#         0.0
#     end
# end

function dkF4Allintqs(p0, ps, qsmax, k, m, T)
    if p0>=0.0
        F4_intcosthqs(p0, ps, qsmax, k, m, T)
    elseif p0<0.0
        -F4_intcosthqs(-p0, ps, qsmax, k, m, T)
    end
end
