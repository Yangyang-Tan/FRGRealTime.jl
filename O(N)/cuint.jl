using LinearAlgebra
using GFlops
using Memoize
using MKL
using BenchmarkTools
BLAS.get_config()


(typeof ∘ typeof ∘ typeof ∘ typeof ∘ typeof ∘ typeof)(1)

@gflops sum(rand(10, 10))


abstract type AbstractArray{T,N} end



function eltype(::Type{<:AbstractArray{T}}) where {T}
    T
end

AbstractArray{Float64,4}

Type{Array}

eltype(AbstractArray{Real,4})


typeof(Type{Real})

Type{Real}

using HCubature
using Plots
ENV["GKSwstype"] = "100"
plot(x -> x + 87, 5, 10)|>display
popdisplay()
scatter(rand(20000), rand(20000)) |> display


isa(1.0,Type{Float64})
begin
    return 1
    return 0
end

@memoize function fact(n)
           n < 0 && error("n must be non-negative")
           n <2  && return 1
           abs(fact(n-1)-fact(n-2))-n
end



function fact2(n)
        n < 0 && error("n must be non-negative")
        n <2  && return 2
        abs(fact2(n-1)-fact2(n-2))-n
end



@time fact2(20)


@time fact(50)




testv=rand(0:BigInt(10),1000)


@benchmark fact.($testv)

@benchmark fact2.($testv)


memoize_cache(fact)