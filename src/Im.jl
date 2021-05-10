F1 = loopfunpp
F2 = loopfunpm
dkF1(p0,ps,k,m,T)=flowpp2(p0,ps,k,m,T)
dkF2(p0,ps,k,m,T)=flowpm(p0,ps,k,m,T)
dkF1(p0, ps, qs, k, m, T)=flowpp_intcosthe(p0, ps, qs, k, m, T)
dkF2(p0, ps, qs, k, m, T)=flowpm_intcosthe(p0, ps, qs, k, m, T)

function deltasumkAll(p0, ps, k, T,Npi)
    if p0>=0.0
        deltasumk(p0, ps, k, T,Npi)
    elseif p0<0.0
        -deltasumk(-p0, ps, k, T,Npi)
    end
end


function dkF1Allintqs(p0, ps, qsmax, k, m, T)
    if p0>=0.0
        flowpp_intcostheqs(p0, ps, qsmax, k, m, T)
    elseif p0<0.0
        -flowpp_intcostheqs(-p0, ps, qsmax, k, m, T)
    end
end

function dkF2Allintqs(p0, ps, qsmax, k, m, T)
    if p0>=0.0
        flowpm_intcostheqs(p0, ps, qsmax, k, m, T)
    elseif p0<0.0
        -flowpm_intcostheqs(-p0, ps, qsmax, k, m, T)
    end
end
