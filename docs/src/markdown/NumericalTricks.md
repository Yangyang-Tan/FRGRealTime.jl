# Delta Function

## Integrate Dirac Delta in $F_1$


### Im Part
There are two types of dirac delta function we need to handle. One is due to the existence of
$\theta\left(p_0-2E_\pi\right)$, the other is due to the existence of $\delta\left(p_0-2E_\pi\right)$.

For the first one, we can see the plot around $p_0=2E_\pi$
```julia
using FRGRealTime,Plots
#k>ps
plot(p0->FRGRealTime.loopfunppfix(p0,1.0,2.0,1.0,2.0),2*Epi(2.0,1.0)-0.01,2*Epi(2.0,1.0)+0.01)
```

![](figures/NumericalTricks_1_1.png)



Considering that our function has the following form:

```math
    F(k)=f(k)\theta\Big(g(k)-p_0\Big)
```
and this relation is only true when $p_0$ is in a infinitesimal neighborhood around $g(k)$.
$F'(k)$ is:

```math
    F'(k)=f(k) g'(k) \delta \Big(g(k)-p_0\Big)+f'(k) \theta \Big(g(k)-p_0\Big)
```

Now we want to integrate back to $F(k)$, any finite part will disappear due to the infinitesimal measure, so:

```math
    F'(k)=f(k) g'(k) \delta \Big(g(k)-p_0\Big)=f(k)\partial_k \theta\Big(g(k)-p_0\Big)
```