using FRGRealTime.MathFuns

function Imlambdaflow_Texture(
    p0,
    q0,
    p0min,
    q0min,
    dp0,
    dq0,
    k,
    m2,
    Ek,
    T,
    Imlambda,
    Relambda,
    ϵ = 0.1f0,
    Nf = 4.0f0,
)
    -2 * (
        TF.ImFb2(p0 - q0, k, m2, T, ϵ) * (
            Imlambda[(dp0+p0-p0min)/dp0, (dq0-q0-q0min)/dq0]^2 +
            Relambda[(dp0+p0-p0min)/dp0, (dq0+q0-q0min)/dq0]^2
        ) +
        TF.ImFb2(p0 + q0, k, m2, T, ϵ) * (
            Imlambda[(dp0+p0-p0min)/dp0, (dq0+q0-q0min)/dq0]^2 +
            Relambda[(dp0+p0-p0min)/dp0, (dq0+q0-q0min)/dq0]^2
        )
    ) / 3 +
    (
        TF.Fb2(k, m2, T) * (
            4 *
            Imlambda[(dp0-p0min+q0)/dp0, -((-dq0 + q0 + q0min) / dq0)] *
            Relambda[(dp0+p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
            (
                Imlambda[(dp0-p0min+q0)/dp0, (dq0+Ek-q0min)/dq0] +
                Imlambda[(dp0-p0min+q0)/dp0, (dq0-Ek-q0min)/dq0]
            ) * (
                Nf * Relambda[(dp0+p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
                2 * Relambda[(dp0+p0-p0min)/dp0, ((dq0+p0-q0min)/dq0)]
            ) +
            4 *
            Imlambda[(dp0+p0-p0min)/dp0, -((-dq0 + p0 + q0min) / dq0)] *
            Relambda[(dp0-p0min+q0)/dp0, (dq0+Ek-q0min)/dq0] +
            (
                Imlambda[(dp0+p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
                Imlambda[(dp0+p0-p0min)/dp0, -((-dq0 + Ek + q0min) / dq0)]
            ) * (
                Nf * Relambda[(dp0-p0min+q0)/dp0, (dq0+Ek-q0min)/dq0] +
                2 * Relambda[(dp0-p0min+q0)/dp0, (dq0+q0-q0min)/dq0]
            )
        )
    ) / 6
end

function Relambdaflow_Texture(
    p0,
    q0,
    p0min,
    q0min,
    dp0,
    dq0,
    k,
    m2,
    Ek,
    T,
    Imlambda,
    Relambda,
    ϵ = 0.1f0,
    Nf = 4.0f0,
)
    2 * (
        Imlambda[(dp0+p0-p0min)/dp0, (dq0+q0-q0min)/dq0]^2 *
        TF.ReFb2(p0 + q0, k, m2, T, ϵ) +
        Imlambda[(dp0+p0-p0min)/dp0, (dq0-q0-q0min)/dq0]^2 *
        TF.ReFb2(p0 - q0, k, m2, T, ϵ) +
        (TF.ReFb2(p0 + q0, k, m2, T, ϵ) + TF.ReFb2(p0 - q0, k, m2, T, ϵ)) *
        Relambda[(dp0+p0-p0min)/dp0, (dq0+q0-q0min)/dq0]^2
    ) / 3 -
    1 / 6 * (
        TF.Fb2(k, m2, T) * (
            2 *
            Imlambda[(dp0+p0-p0min)/dp0, -((-dq0 + p0 + q0min) / dq0)] *
            (
                Imlambda[(dp0-p0min+q0)/dp0, (dq0+Ek-q0min)/dq0] +
                Imlambda[(dp0-p0min+q0)/dp0, -((-dq0 + Ek + q0min) / dq0)]
            ) +
            Imlambda[(dp0+p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] * (
                Nf * Imlambda[(dp0-p0min+q0)/dp0, -((-dq0 + Ek + q0min) / dq0)] +
                2 * Imlambda[(dp0-p0min+q0)/dp0, -((-dq0 + q0 + q0min) / dq0)]
            ) +
            Imlambda[(dp0+p0-p0min)/dp0, -((-dq0 + Ek + q0min) / dq0)] * (
                Nf * Imlambda[(dp0-p0min+q0)/dp0, (dq0+Ek-q0min)/dq0] +
                2 * Imlambda[(dp0-p0min+q0)/dp0, -((-dq0 + q0 + q0min) / dq0)]
            ) -
            2 *
            (
                Nf * Relambda[(dp0+q0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
                2 * Relambda[(dp0+q0-p0min)/dp0, ((dq0+q0-q0min)/dq0)]
            ) *
            Relambda[(dp0-p0min+p0)/dp0, (dq0+Ek-q0min)/dq0] -
            4 *
            Relambda[(dp0+q0-p0min)/dp0, (dq0+Ek-q0min)/dq0] *
            Relambda[(dp0-p0min+p0)/dp0, (dq0+p0-q0min)/dq0] -
            0.5f0 *
            Nf *
            (
                Imlambda[(dp0+p0-p0min)/dp0, -((-dq0 + Ek + q0min) / dq0)] -
                Imlambda[(dp0+p0-p0min)/dp0, ((dq0+Ek-q0min)/dq0)]
            ) *
            (
                Imlambda[(dp0+q0-p0min)/dp0, -((-dq0 + Ek + q0min) / dq0)] -
                Imlambda[(dp0+q0-p0min)/dp0, ((dq0+Ek-q0min)/dq0)]
            )
        )
    )
end




function Twopointflow_Texture(p0, p0min, q0min, dp0, dq0, k, m2, Ek, T, lambda, Nf = 4.0f0)
    (1 / 12) *
    (
        Nf * (
            lambda[(dp0-p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
            lambda[(dp0-p0-p0min)/dp0, (dq0-Ek-q0min)/dq0]
        ) + 4 * lambda[(dp0-p0-p0min)/dp0, (dq0+p0-q0min)/dq0]
    ) *
    TF.Fb1(k, m2, T)
end


# function Twopointflow_Texture(p0, p0min, q0min, dp0, dq0, k, m2, Ek, T, lambda,lambda0, Nf = 4.0f0)
#     (1 / 12) *
#     (
#         Nf * (
#             lambda[(dp0-p0-p0min)/dp0, (dq0+Ek-q0min)/dq0] +
#             lambda[(dp0-p0-p0min)/dp0, (dq0-Ek-q0min)/dq0]
#         )+ 4 * lambda0
#     )*TF.Fb1(k, m2, T)
# end





function Imlambdaflow_GPU(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    (
        -2 * (
            (Imlambdap0q0^2 + Relambdap0q0^2) * TF.ImFb2(p0 + q0, k, m2, T, ϵ) +
            (Imlambdap0mq0^2 + Relambdap0q0^2) * TF.ImFb2(p0 - q0, k, m2, T, ϵ)
        )
    ) / 3 +
    (
        (
            4 * Imlambdaq0mq0 * Relambdap0Ek +
            Imlambdaq0Ek * Nf * Relambdap0Ek +
            Imlambdaq0mEk * Nf * Relambdap0Ek +
            2 * Imlambdaq0Ek * Relambdap0p0 +
            2 * Imlambdaq0mEk * Relambdap0p0 +
            4 * Imlambdap0mp0 * Relambdaq0Ek +
            Imlambdap0Ek * Nf * Relambdaq0Ek +
            Imlambdap0mEk * Nf * Relambdaq0Ek +
            2 * (Imlambdap0Ek + Imlambdap0mEk) * Relambdaq0q0
        ) * TF.Fb2(k, m2, T)
    ) / 6
end





function Relambdaflow_GPU(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    (
        2 * (
            (Imlambdap0mq0^2 + Relambdap0q0^2) * TF.ReFb2(p0 - q0, k, m2, T, ϵ) +
            (Imlambdap0q0^2 + Relambdap0q0^2) * TF.ReFb2(p0 + q0, k, m2, T, ϵ)
        )
    ) / 3 -
    1 / 12 * (
        (
            4 * Imlambdap0mp0 * (Imlambdaq0Ek + Imlambdaq0mEk) +
            4 * Imlambdap0Ek * Imlambdaq0mq0 +
            4 * Imlambdap0mEk * Imlambdaq0mq0 +
            Imlambdap0Ek * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf +
            Imlambdap0mEk * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf -
            4 * (Nf * Relambdap0Ek + 2 * Relambdap0p0) * Relambdaq0Ek -
            8 * Relambdap0Ek * Relambdaq0q0
        ) * TF.Fb2(k, m2, T)
    )
end


function Relambdaflow_GPU_A(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    (
        2 * (
            (Imlambdap0mq0^2 + Relambdap0q0^2) * TF.ReFb2A(p0 - q0, k, m2, T, ϵ) +
            (Imlambdap0q0^2 + Relambdap0q0^2) * TF.ReFb2A(p0 + q0, k, m2, T, ϵ)
        )
    ) / 3
end

function Relambdaflow_GPU_B(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    (
        2 * (
            (Imlambdap0mq0^2 + Relambdap0q0^2) * TF.ReFb2B(p0 - q0, k, m2, T, ϵ) +
            (Imlambdap0q0^2 + Relambdap0q0^2) * TF.ReFb2B(p0 + q0, k, m2, T, ϵ)
        )
    ) / 3 -
    1 / 12 * (
        (
            4 * Imlambdap0mp0 * (Imlambdaq0Ek + Imlambdaq0mEk) +
            4 * Imlambdap0Ek * Imlambdaq0mq0 +
            4 * Imlambdap0mEk * Imlambdaq0mq0 +
            Imlambdap0Ek * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf +
            Imlambdap0mEk * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf -
            4 * (Nf * Relambdap0Ek + 2 * Relambdap0p0) * Relambdaq0Ek -
            8 * Relambdap0Ek * Relambdaq0q0
        ) * TF.Fb2(k, m2, T)
    )
end



function ImRelambdaflow_Kernel(
    Imout,
    Reout,
    Imlambdaflow_fun,
    Relambdaflow_fun,
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    I = CartesianIndices(Imout)
    # idx0=Int32((size(I)[2]+1)/2)
    @inbounds if tid <= length(I)
        i, j = Tuple(I[tid])
        Imout[i, j] = Imlambdaflow_fun(
            p0[i],
            q0[j],
            Imlambdap0q0[i, j],
            Imlambdap0mq0[i, j],
            Imlambdap0mp0[i],
            Imlambdaq0mq0[j],
            Imlambdap0Ek[i],
            Imlambdap0mEk[i],
            Imlambdaq0Ek[j],
            Imlambdaq0mEk[j],
            Relambdap0q0[i, j],
            Relambdap0p0[i],
            Relambdaq0q0[j],
            Relambdap0Ek[i],
            Relambdaq0Ek[j],
            k,
            m2,
            T,
            ϵ,
            Nf,
        )
        Reout[i, j] = Relambdaflow_fun(
            p0[i],
            q0[j],
            Imlambdap0q0[i, j],
            Imlambdap0mq0[i, j],
            Imlambdap0mp0[i],
            Imlambdaq0mq0[j],
            Imlambdap0Ek[i],
            Imlambdap0mEk[i],
            Imlambdaq0Ek[j],
            Imlambdaq0mEk[j],
            Relambdap0q0[i, j],
            Relambdap0p0[i],
            Relambdaq0q0[j],
            Relambdap0Ek[i],
            Relambdaq0Ek[j],
            k,
            m2,
            T,
            ϵ,
            Nf,
        )
    end
    return
end







function Imlambdaflow_Kernelxy(
    Imout,
    Reout,
    Imlambdaflow_fun,
    Relambdaflow_fun,
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    j = threadIdx().y + (blockIdx().y - 1) * blockDim().y
    # I = CartesianIndices(Imout)
    # idx0=Int32((size(I)[2]+1)/2)
    if i<=size(Imout,1) && j<=size(Imout,2)
        Imout[i, j] = Imlambdaflow_fun(
            p0[i],
            q0[j],
            Imlambdap0q0[i, j],
            Imlambdap0mq0[i, j],
            Imlambdap0mp0[i],
            Imlambdaq0mq0[j],
            Imlambdap0Ek[i],
            Imlambdap0mEk[i],
            Imlambdaq0Ek[j],
            Imlambdaq0mEk[j],
            Relambdap0q0[i, j],
            Relambdap0p0[i],
            Relambdaq0q0[j],
            Relambdap0Ek[i],
            Relambdaq0Ek[j],
            k,
            m2,
            T,
            ϵ,
            Nf,
        )
    end
    # Reout[i, j] = Relambdaflow_fun(
    #         p0[i],
    #         q0[j],
    #         Imlambdap0q0[i, j],
    #         Imlambdap0mq0[i, j],
    #         Imlambdap0mp0[i],
    #         Imlambdaq0mq0[j],
    #         Imlambdap0Ek[i],
    #         Imlambdap0mEk[i],
    #         Imlambdaq0Ek[j],
    #         Imlambdaq0mEk[j],
    #         Relambdap0q0[i, j],
    #         Relambdap0p0[i],
    #         Relambdaq0q0[j],
    #         Relambdap0Ek[i],
    #         Relambdaq0Ek[j],
    #         k,
    #         m2,
    #         T,
    #         ϵ,
    #         Nf,
    #     )
    return
end


function Relambdaflow_Kernelxy(
    Imout,
    Reout,
    Imlambdaflow_fun,
    Relambdaflow_fun,
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    Relambdap0q0,
    Relambdap0p0,
    Relambdaq0q0,
    Relambdap0Ek,
    Relambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    j = threadIdx().y + (blockIdx().y - 1) * blockDim().y
    # I = CartesianIndices(Imout)
    # idx0=Int32((size(I)[2]+1)/2)
    if i<=size(Imout,1) && j<=size(Imout,2)
        Reout[i, j] = Relambdaflow_fun(
            p0[i],
            q0[j],
            Imlambdap0q0[i, j],
            Imlambdap0mq0[i, j],
            Imlambdap0mp0[i],
            Imlambdaq0mq0[j],
            Imlambdap0Ek[i],
            Imlambdap0mEk[i],
            Imlambdaq0Ek[j],
            Imlambdaq0mEk[j],
            Relambdap0q0[i, j],
            Relambdap0p0[i],
            Relambdaq0q0[j],
            Relambdap0Ek[i],
            Relambdaq0Ek[j],
            k,
            m2,
            T,
            ϵ,
            Nf,
        )
    end
    return
end


function ReCopy_Kernel(Reout)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    j = threadIdx().y + (blockIdx().y - 1) * blockDim().y
    lx=size(Reout,1)
    if i<=size(Reout,1) && j<=size(Reout,2)
        Reout[i, j] = Reout[min(i,lx-i+1), min(j,lx-j+1)]
    end
    return
end


function ImCopy_Kernel(Imout)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    j = threadIdx().y + (blockIdx().y - 1) * blockDim().y
    lx=size(Imout,1)
    if i<=size(Imout,1) && j<=(size(Imout,2)-1)÷2
        Imout[lx-i+1, lx-j+1] = -Imout[i, j]
    end
    return
end


# function symtrize(A,lengthy)
#     l=Int32((lengthy+1)/2)
#     tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
#     I = CartesianIndices(A)
#     @inbounds if tid <= length(I)
#         i, j = Tuple(I[tid])
#         A[i,j]=A[i,j]
#     end
# end

function ImRelambdaflow_Launch(
    Imout,
    Reout,
    p0,
    q0,
    A::LambdaGrid,
    k,
    m2,
    T,
    ϵ,
    Nf;
    threads = 1000,
    Imlambdaflow_fun = Imlambdaflow_GPU,
    Relambdaflow_fun = Relambdaflow_GPU,
)
    @cuda threads = threads blocks = cld(length(Imout), threads) ImRelambdaflow_Kernel(
        Imout,
        Reout,
        Imlambdaflow_fun,
        Relambdaflow_fun,
        p0,
        q0,
        A.Imlambdap0q0,
        A.Imlambdap0mq0,
        A.Imlambdap0mp0,
        A.Imlambdaq0mq0,
        A.Imlambdap0Ek,
        A.Imlambdap0mEk,
        A.Imlambdaq0Ek,
        A.Imlambdaq0mEk,
        A.Relambdap0q0,
        A.Relambdap0p0,
        A.Relambdaq0q0,
        A.Relambdap0Ek,
        A.Relambdaq0Ek,
        k,
        m2,
        T,
        ϵ,
        Nf,
    )

end




function ImRelambdaflow_Launchxy(
    Imout,
    Reout,
    p0,
    q0,
    A::LambdaGrid,
    k,
    m2,
    T,
    ϵ,
    Nf;
    threads = (32,32),
    Imlambdaflow_fun = Imlambdaflow_GPU,
    Relambdaflow_fun = Relambdaflow_GPU,
)
    @cuda threads = threads blocks = cld.(size(Imout), threads.*(1,2)) Imlambdaflow_Kernelxy(
        Imout,
        Reout,
        Imlambdaflow_fun,
        Relambdaflow_fun,
        p0,
        q0,
        A.Imlambdap0q0,
        A.Imlambdap0mq0,
        A.Imlambdap0mp0,
        A.Imlambdaq0mq0,
        A.Imlambdap0Ek,
        A.Imlambdap0mEk,
        A.Imlambdaq0Ek,
        A.Imlambdaq0mEk,
        A.Relambdap0q0,
        A.Relambdap0p0,
        A.Relambdaq0q0,
        A.Relambdap0Ek,
        A.Relambdaq0Ek,
        k,
        m2,
        T,
        ϵ,
        Nf,
    )
    @cuda threads = threads blocks = cld.(size(Imout), threads.*(2,2)) Relambdaflow_Kernelxy(
        Imout,
        Reout,
        Imlambdaflow_fun,
        Relambdaflow_fun,
        p0,
        q0,
        A.Imlambdap0q0,
        A.Imlambdap0mq0,
        A.Imlambdap0mp0,
        A.Imlambdaq0mq0,
        A.Imlambdap0Ek,
        A.Imlambdap0mEk,
        A.Imlambdaq0Ek,
        A.Imlambdaq0mEk,
        A.Relambdap0q0,
        A.Relambdap0p0,
        A.Relambdaq0q0,
        A.Relambdap0Ek,
        A.Relambdaq0Ek,
        k,
        m2,
        T,
        ϵ,
        Nf,
    )
    @cuda threads = threads blocks = cld.(size(Reout), threads) ReCopy_Kernel(Reout)
    @cuda threads = threads blocks = cld.(size(Imout), threads.*(1,2)) ImCopy_Kernel(Imout)
end



function ImTwopointflow_GPU(Imlambdap0Ek, Imlambdap0mEk, Imlambdap0mp0, k, m2, T, Nf)
    -(1 / 12) * (Nf * (Imlambdap0mEk + Imlambdap0Ek) + 4 * Imlambdap0mp0) * TF.Fb1(k, m2, T)
end
function ReTwopointflow_GPU(Relambdap0p0, Relambdap0Ek, k, m2, T, Nf)
    (1 / 12) * (Nf * (2 * Relambdap0Ek) + 4 * Relambdap0p0) * TF.Fb1(k, m2, T)
end


function ImReTwopointflow_Kernel(
    Imout,
    Reout,
    ImTwopointflow_fun,
    ReTwopointflow_fun,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdap0mp0,
    Relambdap0p0,
    Relambdap0Ek,
    k,
    m2,
    T,
    Nf,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    @inbounds if i <= length(Imout)
        Imout[i] = ImTwopointflow_fun(
            Imlambdap0Ek[i],
            Imlambdap0mEk[i],
            Imlambdap0mp0[i],
            k,
            m2,
            T,
            Nf,
        )
        Reout[i] = ReTwopointflow_fun(Relambdap0p0[i], Relambdap0Ek[i], k, m2, T, Nf)
    end
    return
end
function ImReTwopointflow_Launch(
    Imout,
    Reout,
    A::LambdaGrid,
    k,
    m2,
    T,
    Nf;
    ImTwopointflow_fun = ImTwopointflow_GPU,
    ReTwopointflow_fun = ReTwopointflow_GPU,
    threads = 1000,
)
    @cuda threads = threads blocks = cld(length(Imout), threads) ImReTwopointflow_Kernel(
        Imout,
        Reout,
        ImTwopointflow_fun,
        ReTwopointflow_fun,
        A.Imlambdap0Ek,
        A.Imlambdap0mEk,
        A.Imlambdap0mp0,
        A.Relambdap0p0,
        A.Relambdap0Ek,
        k,
        m2,
        T,
        Nf,
    )
end



function Imlambdaflow_GPU_InvRe(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    InvRelambdap0q0,
    InvRelambdap0p0,
    InvRelambdaq0q0,
    InvRelambdap0Ek,
    InvRelambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    (
        -2 * (
            (Imlambdap0q0^2 + 1 / InvRelambdap0q0^2) * TF.ImFb2(p0 + q0, k, m2, T, ϵ) +
            (Imlambdap0mq0^2 + 1 / InvRelambdap0q0^2) * TF.ImFb2(p0 - q0, k, m2, T, ϵ)
        )
    ) / 3 +
    (
        (
            4 * Imlambdaq0mq0 / InvRelambdap0Ek +
            Imlambdaq0Ek * Nf / InvRelambdap0Ek +
            Imlambdaq0mEk * Nf / InvRelambdap0Ek +
            2 * Imlambdaq0Ek / InvRelambdap0p0 +
            2 * Imlambdaq0mEk / InvRelambdap0p0 +
            4 * Imlambdap0mp0 / InvRelambdaq0Ek +
            Imlambdap0Ek * Nf / InvRelambdaq0Ek +
            Imlambdap0mEk * Nf / InvRelambdaq0Ek +
            2 * (Imlambdap0Ek + Imlambdap0mEk) / InvRelambdaq0q0
        ) * TF.Fb2(k, m2, T)
    ) / 6
end





function Relambdaflow_GPU_InvRe(
    p0,
    q0,
    Imlambdap0q0,
    Imlambdap0mq0,
    Imlambdap0mp0,
    Imlambdaq0mq0,
    Imlambdap0Ek,
    Imlambdap0mEk,
    Imlambdaq0Ek,
    Imlambdaq0mEk,
    InvRelambdap0q0,
    InvRelambdap0p0,
    InvRelambdaq0q0,
    InvRelambdap0Ek,
    InvRelambdaq0Ek,
    k,
    m2,
    T,
    ϵ,
    Nf,
)
    -(
        2 * (
            (InvRelambdap0q0^2 * Imlambdap0mq0^2 + 1) * TF.ReFb2(p0 - q0, k, m2, T, ϵ) +
            (InvRelambdap0q0^2 * Imlambdap0q0^2 + 1) * TF.ReFb2(p0 + q0, k, m2, T, ϵ)
        )
    ) / 3 +
    1 / 12 *
    InvRelambdap0q0^2 *
    (
        (
            4 * Imlambdap0mp0 * (Imlambdaq0Ek + Imlambdaq0mEk) +
            4 * Imlambdap0Ek * Imlambdaq0mq0 +
            4 * Imlambdap0mEk * Imlambdaq0mq0 +
            Imlambdap0Ek * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf +
            Imlambdap0mEk * (Imlambdaq0Ek + Imlambdaq0mEk) * Nf -
            4 * ((Nf / InvRelambdap0Ek) + (2 / InvRelambdap0p0)) / InvRelambdaq0Ek -
            8 / (InvRelambdap0Ek * InvRelambdaq0q0)
        ) * TF.Fb2(k, m2, T)
    )
end


function ReTwopointflow_GPU_InvRe(InvRelambdap0p0, InvRelambdap0Ek, k, m2, T, Nf)
    (1 / 12) * (Nf * (2 / InvRelambdap0Ek) + 4 / InvRelambdap0p0) * TF.Fb1(k, m2, T)
end
