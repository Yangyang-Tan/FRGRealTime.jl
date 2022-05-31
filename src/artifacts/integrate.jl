function quadgk_PV(f, a, c, b; kwargs...)
    if b > c > a
        fc = f(c)
        g(x) = 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc)
        return quadgk(g, a, prevfloat(c); kwargs...)[1] +
               quadgk(g, nextfloat(c), b; kwargs...)[1] +
               fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
    else
        return quadgk(x -> 2 * x * (-c^2 + x^2)^-1 * f(x), a, b; kwargs...)[1]
    end
end

function quadgk_PV2(f, a, c, b; initdiv=1,kwargs...)
    fc = f(c)
    function g(x)
        if x == c
            2 * x * (-nextfloat(c)^2 + x^2)^-1 * (f(x) - fc)
        else
            return 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc)
        end
    end
    return quadgk(g, a, b; order=ceil(Int64,initdiv/5),kwargs...)[1] + fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
end




function quadgk_PV3(f, a, c, b;initdiv=1,kwargs...)
    if a < c < b
        fc = f(c)
        function g(x)
            if x == c
                2 * x * (-nextfloat(c)^2 + x^2)^-1 * (f(x) - fc)
            else
                return 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc)
            end
        end
        return quadgk(g, a, b; order=ceil(Int64,initdiv/5),kwargs...)[1]
        # + fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
    else
        fc2 = f(nextfloat(a))
        function g2(x)
            if x == c
                2 * x * (-nextfloat(c)^2 + x^2)^-1 * (f(x) - fc2)
            else
                return 2 * x * (-c^2 + x^2)^-1 * (f(x) - fc2)
            end
        end
        return quadgk(g2, a, b; order=ceil(Int64,initdiv/5),kwargs...)[1]
        # + fc * log(abs((b^2 - c^2) / (a^2 - c^2)))
    end
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
