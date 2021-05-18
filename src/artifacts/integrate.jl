function quadgk_PV(f, a, c, b; kwargs...)
    fc = f(c)
    g(x) = 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc)
    return quadgk(g, a, prevfloat(c); kwargs...)[1] +
           quadgk(g, nextfloat(c), b; kwargs...)[1] +
           fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
end

function hquadrature_PV(f, a, c, b; kwargs...)
    fc = f(c)
    g(x) = 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc)
    return hquadrature(g, a, prevfloat(c); kwargs...)[1] +
           hquadrature(g, nextfloat(c), b; kwargs...)[1] +
           fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
end

# function quadgk_PV2(f, c; kwargs...)
#     g(x) = (2*f(x))/(x)
#     return quadgk(g, 0.0, Inf; kwargs...)[1]
# end



# quadgk_PV(x->sin(x),10.0,15.0,20.0,atol=0.0,rtol=1e-4)
# PvdkF1Tildeps(p0, k, kprim, m, T)
