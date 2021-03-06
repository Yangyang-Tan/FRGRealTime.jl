## Check F1 & F2

load package
```julia
using FRGRealTime,Plots,HCubature
```





Define integrate $\tilde{\partial_k}F1$ function:
```julia
testloopfunpp(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpp(p0,ps,x,m,T),0.0,k,initdiv=2000)[1]+FRGRealTime.loopfunpp(p0,ps,0.0,m,T)
```

```
testloopfunpp (generic function with 1 method)
```





Define integrate $\tilde{\partial_k}F2$ function:
```julia
testloopfunpm(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpm(p0,ps,x,m,T),0.0,k)[1]+FRGRealTime.loopfunpm(p0,ps,0.0,m,T)
```

```
testloopfunpm (generic function with 1 method)
```







Define integrate $\tilde{\partial_k}F1$ function for m<0:
```julia
testloopfunpp2(p0,ps,k,m,T)=-hquadrature(x->FRGRealTime.flowpp(p0,ps,x,m,T),k,400.0,rtol=1e-14,atol=1e-14,initdiv=10000)[1]+FRGRealTime.loopfunpp(p0,ps,400.0,m,T)
```

```
testloopfunpp2 (generic function with 1 method)
```





Define integrate $\tilde{\partial_k}F2$ function:
```julia
testloopfunpm2(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpm(p0,ps,x,m,T),2.0,k,initdiv=1000)[1]+FRGRealTime.loopfunpm(p0,ps,2.0,m,T)
```

```
testloopfunpm2 (generic function with 1 method)
```









### F1 Checking Results

 $k>ps,m>0,p_0\in (1.0,400.0)$
```julia
plot(p0->testloopfunpp(p0,10.0,20.0,1.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,1.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_6_1.png)




  $k>ps,m>0,p_0\in (1.0,10.0)$
```julia
plot(p0->testloopfunpp(p0,10.0,20.0,1.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,1.0,145.0),1.0,20.0,label="F1")
```

![](figures/Check_7_1.png)





 $ps/2<k<ps,m>0$
```julia
plot(p0->testloopfunpp(p0,10.0,8.0,1.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,8.0,1.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_8_1.png)



  $k<ps/2,m>0$
```julia
plot(p0->testloopfunpp(p0,10.0,2.0,1.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,2.0,1.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_9_1.png)



 $k>ps,m<0,p_0\in (1.0,400.0)$
```julia
plot(p0->testloopfunpp2(p0,10.0,20.0,-3.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,-3.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_10_1.png)



 $ps/2<k<ps,m<0$
```julia
plot(p0->testloopfunpp2(p0,10.0,8.0,-3.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,8.0,-3.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_11_1.png)




 $k<ps/2,m>0$
```julia
plot(p0->testloopfunpp2(p0,10.0,2.0,-3.0,145.0),1.0,400.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,2.0,-3.0,145.0),1.0,400.0,label="F1")
```

![](figures/Check_12_1.png)





 $k<ps/2,m<0$
```julia
plot(p0->testloopfunpp2(p0,10.0,4.0,-3.0,145.0),9.25,9.5,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpp(p0,10.0,4.0,-3.0,145.0),9.25,9.5,label="F1")
```

![](figures/Check_13_1.png)




### F2 Checking Results

 $k>ps,m>0$
```julia
plot(p0->testloopfunpm(p0,10.0,20.0,1.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,20.0,1.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_14_1.png)




 $ps/2<k<ps,m>0$
```julia
plot(p0->testloopfunpm(p0,10.0,8.0,1.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,8.0,1.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_15_1.png)



 $k<ps/2,m>0$
```julia
plot(p0->testloopfunpm(p0,10.0,2.0,1.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,2.0,1.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_16_1.png)




 $k>ps,m<0$
```julia
plot(p0->testloopfunpm2(p0,10.0,20.0,-3.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,20.0,-3.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_17_1.png)



 $ps/2<k<ps,m<0$
```julia
plot(p0->testloopfunpm2(p0,10.0,8.0,-3.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,8.0,-3.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_18_1.png)



 $k<ps/2,m<0$
```julia
plot(p0->testloopfunpm2(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label="Integrate flow",xaxis="p0")
plot!(p0->FRGRealTime.loopfunpm(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label="F2")
```

![](figures/Check_19_1.png)




## Check Vertex

The flow equation of vertex has the following form:
```math
\begin{aligned}
\partial_kV_k=&\lambda_{4\pi}(k)^2 \tilde{\partial_k}\Bigg\{3\pi\Bigg.
(2 + N) \Big[\Big.
F_1(p_0-q_0,|\mathbf{p}_s+\mathbf{q}_s|,k,m(k))+F_2(p_0-q_0,|\mathbf{p}_s+\mathbf{q}_s|,k,m(k))\\
&\Bigg.\Big.+F_1(p_0+q_0,|\mathbf{p}_s+\mathbf{q}_s|,k,m(k))+F_2(p_0+q_0,|\mathbf{p}_s+\mathbf{q}_s|,k,m(k))\Big]\Bigg\}
\end{aligned}
```

Define Vertex testing function, we set $q_0=E_\pi(k,m)$:
```julia
testVIm(p0, ps, k, m, T, lam4pik)=-hquadrature(x->FRGRealTime.dkVIm(p0, ps, Epi(k,m), x, m, T, 4.0, lam4pik),k,400.0,atol=1e-6,rtol=1e-6,initdiv=1000)[1]
```

```
testVIm (generic function with 1 method)
```





Define a check plot function with `p0` is x-axis:
```julia
function checkplotVIm(ps,qs,costh, k, m, T, lam4pik)
  plot(p0->testVIm(p0, sqrt(ps^2+qs^2+2*ps*qs*costh), k, m, T, lam4pik),1.0,100.0,label="Integrate Flow",xaxis="p0")
  plot!(p0->FRGRealTime.VImSimple(p0, sqrt(ps^2+qs^2+2*ps*qs*costh), Epi(k,m), k, m, T,4.0, lam4pik,400.0),1.0,100.0,label="F1&F2")
end
```

```
checkplotVIm (generic function with 1 method)
```





$m<0,k>ps/2,\cos\theta=1.0,ps=0.0$
```julia
checkplotVIm(0.0,5.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_22_1.png)



$m<0,k>ps/2,\cos\theta=1.0,ps=0.0,qs=8,k=10$
```julia
checkplotVIm(0.0,8.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_23_1.png)



$m<0,k>ps/2,\cos\theta=1.0,ps=0.0,qs=9,k=10$
```julia
checkplotVIm(0.0,9.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_24_1.png)



$m<0,k>ps/2,\cos\theta=1.0,ps=0.0,qs=ps$
```julia
checkplotVIm(0.0,10.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_25_1.png)







$m<0,k>ps/2,\cos\theta=-1.0,ps=4.0,k=10.0$
```julia
checkplotVIm(4.0,5.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_26_1.png)



$m<0,k>ps/2,\cos\theta=0.0,ps=4.0,k=10.0$
```julia
checkplotVIm(4.0,5.0,0.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_27_1.png)




$m<0,k>ps/2,\cos\theta=1.0,ps=4.0,k=10.0$
```julia
checkplotVIm(4.0,5.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_28_1.png)



$m<0,k>ps/2,\cos\theta=-1.0,ps=10.0,k=10.0$
```julia
checkplotVIm(10.0,5.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_29_1.png)



$m<0,k>ps/2,\cos\theta=0.0,ps=10.0,k=10.0$
```julia
checkplotVIm(10.0,5.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_30_1.png)



$m<0,k>ps/2,\cos\theta=1.0,ps=10.0,k=10.0$
```julia
checkplotVIm(10.0,5.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_31_1.png)









$m<0,k=ps/2,\cos\theta=-1.0$
```julia
checkplotVIm(20.0,5.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_32_1.png)



$m<0,k=ps/2,\cos\theta=0.0$
```julia
checkplotVIm(20.0,5.0,0.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_33_1.png)



$m<0,k=ps/2,\cos\theta=1.0$
```julia
checkplotVIm(20.0,5.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_34_1.png)







$m<0,k<ps/2,\cos\theta=-1.0$
```julia
checkplotVIm(40.0,5.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_35_1.png)



$m<0,k<ps/2,\cos\theta=0.0$
```julia
checkplotVIm(40.0,5.0,0.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_36_1.png)



$m<0,k<ps/2,\cos\theta=1.0$
```julia
checkplotVIm(40.0,5.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_37_1.png)






$m<0,k<ps/2,\cos\theta=-1.0,qs=k$
```julia
checkplotVIm(40.0,10.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_38_1.png)



$m<0,k<ps/2,\cos\theta=0.0,qs=k$
```julia
checkplotVIm(40.0,10.0,0.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_39_1.png)



$m<0,k<ps/2,\cos\theta=1.0,qs=k$
```julia
checkplotVIm(40.0,10.0,1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_40_1.png)





$m<0,k=ps/2,\cos\theta=-1.0,qs=k$
```julia
checkplotVIm(20.0,10.0,-1.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_41_1.png)



$m<0,k=ps/2,\cos\theta=0.0,qs=k$
```julia
checkplotVIm(20.0,10.0,0.0,10.0,-2.0,145.0,1.0)
```

![](figures/Check_42_1.png)











## Check Two Point function


Define mass function and $\lambda_4\pi$ function:
```julia
using Dierckx
testmfun(x)=-2.0
testlamfun(x)=1.0
testmfun2=Spline1D(collect(0.1:0.5:800.0),fill(-2.0,length(collect(0.1:0.5:800.0))))
testlamfun2=Spline1D(collect(0.1:0.5:800.0),fill(1.0,length(collect(0.1:0.5:800.0))))
```

```
Spline1D(knots=[0.1,1.1 ??? 798.6,799.6] (1598 elements), k=3, extrapolation=
"nearest", residual=0.0)
```






$p_s=0.5, m=-2,T=145,p_0\in (1,100)$
```julia
v1=collect(1.0:2.0:100.0)
propout1=(p0->FRGRealTime.propImintqs(p0, 0.5, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propdeltaout1=-(p0->FRGRealTime.propImintqs_delta1(p0, 0.5, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propsimpleout1=FRGRealTime.propImSimple.(v1, 0.5, 145.0, 5.0,400.0, 4.0, -2.0, 1.0)
plot(v1,propout1+propdeltaout1,lable="Integrate Flow",xaxis="p0")
plot!(v1,propsimpleout1,lable="F1 & F2")
```

![](figures/Check_44_1.png)




$p_s=0.5, m=-2,T=145,p_0\in (1,400)$
```julia
v1=collect(1.0:5.0:400.0)
propout1=(p0->FRGRealTime.propImintqs(p0, 0.5, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propdeltaout1=-(p0->FRGRealTime.propImintqs_delta1(p0, 0.5, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propsimpleout1=FRGRealTime.propImSimple.(v1, 0.5, 145.0, 5.0,400.0, 4.0, -2.0, 1.0)
plot(v1,propout1+propdeltaout1,lable="Integrate Flow",xaxis="p0")
plot!(v1,propsimpleout1,lable="F1 & F2")
```

![](figures/Check_45_1.png)




$p_s=100.0, m=-2,T=145,p_0\in (1,100)$
```julia
v1=collect(1.0:2.0:100.0)
propout1=(p0->FRGRealTime.propImintqs(p0, 100.0, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propdeltaout1=-(p0->FRGRealTime.propImintqs_delta1(p0, 100.0, 145.0,5.0,400.0, 4.0, testmfun2, testlamfun2)).(v1)
propsimpleout1=FRGRealTime.propImSimple.(v1, 100.0, 145.0, 5.0,400.0, 4.0, -2.0, 1.0)
plot(v1,propout1+propdeltaout1,lable="Integrate Flow",xaxis="p0")
plot!(v1,propsimpleout1,lable="F1 & F2")
```

![](figures/Check_46_1.png)
