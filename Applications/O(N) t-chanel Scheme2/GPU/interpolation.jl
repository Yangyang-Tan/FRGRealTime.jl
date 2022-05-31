using CUDA
using Interpolations
device!(1)
N = 1
T = Float32
using Test

@testset "interpolations" begin
    @testset "$interpolate $T" for T in (Float16, Float32,)
        @testset "$(N)D" for N in 1:2
            cpu_src = rand(T, fill(10, N)...)
            cpu_idx = [tuple(rand(1:0.1:10, N)...) for _ in 1:10]

            gpu_src = CuTextureArray(CuArray(cpu_src))
            gpu_idx = CuArray(cpu_idx)

            @testset "nearest neighbour" begin
                cpu_dst = similar(cpu_src, size(cpu_idx))
                cpu_int = interpolate(cpu_src, BSpline(Constant()))
                broadcast!(cpu_dst, cpu_idx, Ref(cpu_int)) do idx, int
                    int(idx...)
                end

                gpu_dst = CuArray{T}(undef, size(cpu_idx))
                gpu_tex = CuTexture(gpu_src; interpolation = CUDA.NearestNeighbour())
                broadcast!(gpu_dst, gpu_idx, Ref(gpu_tex)) do idx, tex
                    tex[idx...]
                end

                @test cpu_dst ≈ Array(gpu_dst)
            end

            @testset "linear interpolation" begin
                cpu_dst = similar(cpu_src, size(cpu_idx))
                cpu_int = interpolate(cpu_src, BSpline(Linear()))
                broadcast!(cpu_dst, cpu_idx, Ref(cpu_int)) do idx, int
                    int(idx...)
                end

                gpu_dst = CuArray{T}(undef, size(cpu_idx))
                gpu_tex = CuTexture(gpu_src; interpolation = CUDA.LinearInterpolation())
                broadcast!(gpu_dst, gpu_idx, Ref(gpu_tex)) do idx, tex
                    tex[idx...]
                end

                @test cpu_dst ≈ Array(gpu_dst) rtol = 0.01
            end

            N < 3 && @testset "cubic interpolation" begin
                cpu_dst = similar(cpu_src, size(cpu_idx))
                cpu_int = interpolate(cpu_src, BSpline(Cubic(Line(OnGrid()))))
                broadcast!(cpu_dst, cpu_idx, Ref(cpu_int)) do idx, int
                    int(idx...)
                end

                gpu_dst = CuArray{T}(undef, size(cpu_idx))
                gpu_tex = CuTexture(gpu_src; interpolation = CUDA.CubicInterpolation())
                broadcast!(gpu_dst, gpu_idx, Ref(gpu_tex)) do idx, tex
                    tex[idx...]
                end

                # FIXME: these results, although they look OK in an image,
                #        do not match the output from Interpolations.jl
                @test_skip cpu_dst ≈ Array(gpu_dst)
            end
        end
    end
end

src = rand(T, fill(10, N)...)

# indices we want to interpolate
idx = [tuple(rand(1:0.1:10, N)...) for _ in 1:10]

# upload to the GPU



gpu_src = CuArray(src)
gpu_idx = CuArray(idx)

# create a texture array for optimized fetching
# this is required for N=1, optional for N=2 and N=3
gpu_src = CuTextureArray(gpu_src)

# interpolate using a texture
gpu_dst = CuArray{T}(undef, size(gpu_idx))
gpu_tex = CuTexture(gpu_src; interpolation = CUDA.LinearInterpolation())
broadcast!((idx, tex) -> tex[idx...], gpu_dst, gpu_idx, Ref(gpu_tex))


# back to the CPU
dst = Array(gpu_dst)


myf2(x, y) = sin(x) + cos(y) + log(x * y) + cos(x / y)



myrangex = 0.2f0:0.2f0:800.0f0
myrangey = 0.2f0:0.2f0:800.0f0


myrange2x = 0.3f0:0.2f0:800.0f0
myrange2y = 0.3f0:0.2f0:800.0f0



cpu_src = myf2.(myrangex, myrangey')

surface(cpu_src)


cpu_idxx = ((myrange2x |> collect) ./ 0.2f0).+1.0f0
cpu_idxy = ((myrange2y |> collect) ./ 0.2f0).+1.0f0


gpu_src = CuArray(cpu_src)
gpu_idxx = CuArray(cpu_idxx)
gpu_idxy = CuArray(cpu_idxx)


gpu_dst = CuArray{T}(undef, (length(gpu_idxx), length(gpu_idxy)))
gpu_src_tex=CuTextureArray(gpu_src)
gpu_src_tex2=CuTextureArray(gpu_src)



@benchmark CUDA.@sync CuTexture(gpu_src_tex; interpolation = CUDA.NearestNeighbour())



gpu_src_tex2=CuTextureArray(gpu_src)

gpu_tex=CuTexture(gpu_src_tex; interpolation = CUDA.LinearInterpolation())

gpu_dst = CuArray{T}(undef, (length(gpu_idxx)))

broadcast!(gpu_dst, gpu_idxx, gpu_idxy', Ref(gpu_tex)) do idx, idy, tex
    tex[idx-1.0f0, idy-1.0f0]
end



broadcast!(gpu_dst, gpu_idxx, gpu_idxy',Ref(gpu_tex)) do idx, idy, tex
    tex[idx-1.0f0, idy-1.0f0]
end





cpu_int = interpolate(cpu_src, BSpline(Linear()))


surface(Array(gpu_dst) - myf2.(myrange2x, myrange2y'))

Array(gpu_dst) - myf2.(myrange2x, myrange2y')


surface(cpu_int.(cpu_idxx, cpu_idxy') - myf2.(myrange2x, myrange2y'))


using BenchmarkTools
gpuary=CUDA.randn(4000,4000)
texary1=CuTextureArray(gpuary)

begin
    texary=similar()
end