Load package and data
```julia
using Distributed
using SharedArrays
addprocs(8)
nprocs()
@everywhere using FRGRealTime, DelimitedFiles,Dierckx,Plots
@everywhere kdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/k_pde_break.dat")[:,1]|>reverse
@everywhere lamdata=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/lamdak_pde_break.dat")[:,1]|>reverse
@everywhere m2data=readdlm("/home/tyy/Documents/CTP-fRG-Test/PDE/m2k_pde_break.dat")[:,1]|>reverse
@everywhere lamfun=Spline1D(kdata,lamdata)
@everywhere m2fun=Spline1D(kdata,m2data)
```




prerun & `p0` range
```julia
@everywhere FRGRealTime.propImsimpleintqs(4.0, 10.0, 145.0,1.0,400.0, 4.0, m2fun, lamfun)
v1=SharedArray(collect(0.1:3.0:400.0))
```

```
134-element SharedArrays.SharedVector{Float64}:
   0.1
   3.1
   6.1
   9.1
  12.1
  15.1
  18.1
  21.1
  24.1
  27.1
   ⋮
 375.1
 378.1
 381.1
 384.1
 387.1
 390.1
 393.1
 396.1
 399.1
```





two point function without type1-delta:
```julia
outv1=pmap(
    p0 -> FRGRealTime.propImsimpleintqs(
        p0,
        0.05,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),v1
)
```

```
Error: On worker 2:
UndefVarError: ##WeaveSandBox#261 not defined
Stacktrace:
  [1] deserialize_module
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:965
  [2] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:864
  [3] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782
  [4] deserialize_datatype
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:1287
  [5] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:835
  [6] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782
  [7] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:842
  [8] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782 [inlined]
  [9] deserialize_msg
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/messages.jl:87
 [10] #invokelatest#2
    @ ./essentials.jl:708 [inlined]
 [11] invokelatest
    @ ./essentials.jl:706 [inlined]
 [12] message_handler_loop
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/process_messages.jl:169
 [13] process_tcp_streams
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/process_messages.jl:126
 [14] #99
    @ ./task.jl:411
```





two function from type1-delta:
```julia
outv2=pmap(
    p0 -> FRGRealTime.propImsimpleintqs_delta1(
        p0,
        0.05,
        145.0,
        1.0,
        400.0,
        4.0,
        m2fun,
        lamfun,
    ),v1
)
```

```
Error: On worker 6:
UndefVarError: ##WeaveSandBox#261 not defined
Stacktrace:
  [1] deserialize_module
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:965
  [2] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:864
  [3] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782
  [4] deserialize_datatype
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:1287
  [5] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:835
  [6] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782
  [7] handle_deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:842
  [8] deserialize
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Serialization/src/Serialization.jl:782 [inlined]
  [9] deserialize_msg
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/messages.jl:87
 [10] #invokelatest#2
    @ ./essentials.jl:708 [inlined]
 [11] invokelatest
    @ ./essentials.jl:706 [inlined]
 [12] message_handler_loop
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/process_messages.jl:169
 [13] process_tcp_streams
    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6
/Distributed/src/process_messages.jl:126
 [14] #99
    @ ./task.jl:411
```





two point function total plot:
```julia
plot(v1,outv1.-outv2)
```

```
Error: UndefVarError: outv1 not defined
```





two point function without type1-delta plot:
```julia
plot(v1,outv1)
```

```
Error: UndefVarError: outv1 not defined
```





two function from type1-delta plot:
```julia
plot(v1,outv2)
```

```
Error: UndefVarError: outv2 not defined
```


