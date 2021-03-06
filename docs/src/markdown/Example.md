# Theory Preset

## Two point function flow

The flow of two point function is:

```math
\begin{aligned}
&\partial_{\tau} \Gamma_{k, \phi_{q} \phi_{c}}^{(2)}\left(p_{0},p_s\right) \\
&=-\frac{1}{24} \frac{(N+2)}{(2 \pi)^{2}}\left[-\frac{\coth\left(\frac{E_{\pi, k}(k)}{2 T}\right)}{\left(E_{\pi, k}(k)\right)^{3}}-\frac{\operatorname{csch}^{2}\left(\frac{E_{\pi, k}(k)}{2 T}\right)}{2 T\left(E_{\pi, k}(k)\right)^{2}}\right] \\
&\times\left(2 k^{2}\right) \int_{0}^{k} dq  q_s^{2} \int_{-1}^{1} d \cos \theta
\end{aligned}
```



## Load Package
```julia
using FRGRealTime, Plots
```






### $F_1$ & $F_2$ at zero momentum


 $F_1$ will goes to zero when the width of $\delta$ function is narrow.
```julia
plot(δ->FRGRealTime.loopfunppfix(1e-6,1e-5,1.0,0.5,2.0,δ),1e-10,1e-2)
```

![](figures/Example_2_1.png)



 At zero momentum, $F2$ will goes to zero when $p_0$ is far small than $p_s$, or
$p_0>p_s$.


 $p_s=1e-5$ case:
```julia
plot(p0->FRGRealTime.loopfunpm(p0,1e-5,1.0,0.5,2.0),1e-8,1.5*1e-5)
```

![](figures/Example_3_1.png)




 $p_s=1e-10$ case:
```julia
plot(p0->FRGRealTime.loopfunpm(p0,1e-10,1.0,0.5,2.0),1e-12,1.5*1e-10)
```

![](figures/Example_4_1.png)
