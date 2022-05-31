module MathFuns

export Eb, deltafun, deltadxfun, dabs, PrincipalIntdelta, PrincipalIntdeltadx

"""
    Eb(x,m)

compute ``\\sqrt{(x^2+m)}``
"""
function Eb(x, m2)
    sqrt(x^2 + m2)
end


"""
    deltafun(x,dϵ=0.02)
compute ``\\delta`` function with the approximation ``\\frac{\\epsilon}{\\pi(\\epsilon^2+x^2)}``
and ``\\epsilon`` is set to be ``0.02`` by default
"""
function deltafun(x, ϵ = 1.0,a=1.0f0)
    ϵ / (pi * (ϵ^2 + x^2))
    # if abs(x)<=a*ϵ
    #     (-2 * x * ϵ) / (pi * (x^2 + ϵ^2)^2)
    # else
    #     return 0.0f0
    # end
end
# const Pi=3.141592653589
function deltadxfun(x, ϵ = 1.0f0,a=1.0f0)
    (-2 * x * ϵ) / (pi * (x^2 + ϵ^2)^2)
    # if abs(x)<=a*ϵ
    #     (-2 * x * ϵ) / (pi * (x^2 + ϵ^2)^2)
    # else
    #     return 0.0f0
    # end
    # -1/4*(x*exp(-1/4*x^2/ϵ))/(sqrt(pi*1.0f0)*ϵ^(3/2))
    # Float32(-1/4*(x*exp(-1/4*x^2/ϵ))/(sqrt(pi)*ϵ^(3/2)))
    # if 0.0f0>x>0.5f0*ϵ
    #     return -1/ϵ^2
    # elseif 0.5f0*ϵ>x>0.0f0
    #     return 1/ϵ^2
    # else
    #     return 0.0f0
    # end
end

function dabs(x)
    if x > 0.0
        return 1.0
    elseif x < 0.0
        return -1.0
    elseif x == 0.0
        return 0.0
    end
end

function PrincipalIntdelta(p0, x, ϵ)
    (-p0+x)/((p0-x)^2+ϵ^2)+(p0+x)/((p0+x)^2+ϵ^2)
    # 0.0f0
end

function PrincipalIntdeltadx(p0, x, ϵ)
    ((p0-x)^2-ϵ^2)/((p0-x)^2+ϵ^2)^2+((p0+x)^2-ϵ^2)/((p0+x)^2+ϵ^2)^2
    # (p0-x)^2/((p0-x)^2+ϵ^2)^2+(p0+x)^2/((p0+x)^2+ϵ^2)^2
    # (2*(p0^2-x^2)^2*(p0^2+x^2))/(p0^4+2*p0^2*(-x^2+ϵ^2)+(x^2+ϵ^2)^2)^2
    # 0.0f0
end


# function PrincipalIntdelta(p0, x, ϵ)
#     (2*x*(-p0^2+x^2+ϵ^2))/(p0^4+2*p0^2*(-x^2+ϵ^2)+(x^2+ϵ^2)^2)
# end


# function PrincipalIntdeltadx(p0, x, ϵ)
#     (2*(p0^6*pi+(x^2+ϵ^2)^2*(-2*x*ϵ+pi*(x-ϵ)*(x+ϵ))+p0^4*(-2*x*ϵ+pi*(-x^2+ϵ^2))-p0^2*(-4*x^3*ϵ+4*x*ϵ^3+pi*(x^4+10*x^2*ϵ^2+ϵ^4))))/(pi*(p0^4+2*p0^2*(-x^2+ϵ^2)+(x^2+ϵ^2)^2)^2)
# end


end
