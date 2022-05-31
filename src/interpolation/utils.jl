ratio(num::Integer, denom::Integer) = SimpleRatio(promote(num, denom)...)
ratio(num, denom) = num/denom
fast_trunc(::Type{Int}, x) = unsafe_trunc(Int, x)
fast_trunc(::Type{Int}, x::Rational) = x.num ÷ x.den
