function Imlambdaflow_Zero(p0,q0,k,m2,T,Relambda,Nf=4.0f0,ϵ=0.1f0)
    -2*(Relambda^2*(TF.ImFb2(p0-q0,k,m2,T,ϵ)+TF.ImFb2(p0+q0,k,m2,T,ϵ)))/3
end


function Relambdaflow_Zero(p0,q0,k,m2,T,Relambda,Nf=4.0f0,ϵ=0.1f0)
    (Relambda^2*((4+Nf)*TF.Fb2(k,m2,T)+2*(TF.ReFb2(p0-q0,k,m2,T,ϵ)+TF.ReFb2(p0+q0,k,m2,T,ϵ))))/3
end


function Twopointflow_Zero(lambdamp0Ek,lambdamp0mEk,lambdamp0p0, k, m2, T, Nf = 4.0f0)
    (1 / 12) *
    (
        Nf * (
            lambdamp0Ek +
            lambdamp0mEk
        )+ 4 * lambdamp0p0
    )*TF.Fb1(k, m2, T)
end

# function Twopointflow_Zero(lambdamp0Ek,lambdamp0mEk,lambdamp0p0, k, m2, T, Nf = 4.0f0)
#     (1 / 12) *
#     (
#         Nf * (
#             lambdamp0Ek +
#             lambdamp0mEk
#         )+ 2 * (lambdamp0Ek+lambdamp0mEk)
#     )*TF.Fb1(k, m2, T)
# end
