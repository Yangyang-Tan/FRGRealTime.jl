var documenterSearchIndex = {"docs":
[{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Load package and data","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"using Distributed\nusing SharedArrays\naddprocs(8)\nnprocs()\n@everywhere using FRGRealTime, DelimitedFiles,Dierckx,Plots\n@everywhere kdata=readdlm(\"/home/tyy/Documents/CTP-fRG-Test/PDE/k_pde_break.dat\")[:,1]|>reverse\n@everywhere lamdata=readdlm(\"/home/tyy/Documents/CTP-fRG-Test/PDE/lamdak_pde_break.dat\")[:,1]|>reverse\n@everywhere m2data=readdlm(\"/home/tyy/Documents/CTP-fRG-Test/PDE/m2k_pde_break.dat\")[:,1]|>reverse\n@everywhere lamfun=Spline1D(kdata,lamdata)\n@everywhere m2fun=Spline1D(kdata,m2data)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"prerun & p0 range","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"@everywhere FRGRealTime.propImsimpleintqs(4.0, 10.0, 145.0,1.0,400.0, 4.0, m2fun, lamfun)\nv1=SharedArray(collect(0.1:3.0:400.0))","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"134-element SharedArrays.SharedVector{Float64}:\n   0.1\n   3.1\n   6.1\n   9.1\n  12.1\n  15.1\n  18.1\n  21.1\n  24.1\n  27.1\n   ⋮\n 375.1\n 378.1\n 381.1\n 384.1\n 387.1\n 390.1\n 393.1\n 396.1\n 399.1","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"two point function without type1-delta:","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"outv1=pmap(\n    p0 -> FRGRealTime.propImsimpleintqs(\n        p0,\n        0.05,\n        145.0,\n        1.0,\n        400.0,\n        4.0,\n        m2fun,\n        lamfun,\n    ),v1\n)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Error: On worker 2:\nUndefVarError: ##WeaveSandBox#261 not defined\nStacktrace:\n  [1] deserialize_module\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:965\n  [2] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:864\n  [3] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782\n  [4] deserialize_datatype\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:1287\n  [5] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:835\n  [6] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782\n  [7] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:842\n  [8] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782 [inlined]\n  [9] deserialize_msg\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/messages.jl:87\n [10] #invokelatest#2\n    @ ./essentials.jl:708 [inlined]\n [11] invokelatest\n    @ ./essentials.jl:706 [inlined]\n [12] message_handler_loop\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/process_messages.jl:169\n [13] process_tcp_streams\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/process_messages.jl:126\n [14] #99\n    @ ./task.jl:411","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"two function from type1-delta:","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"outv2=pmap(\n    p0 -> FRGRealTime.propImsimpleintqs_delta1(\n        p0,\n        0.05,\n        145.0,\n        1.0,\n        400.0,\n        4.0,\n        m2fun,\n        lamfun,\n    ),v1\n)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Error: On worker 6:\nUndefVarError: ##WeaveSandBox#261 not defined\nStacktrace:\n  [1] deserialize_module\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:965\n  [2] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:864\n  [3] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782\n  [4] deserialize_datatype\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:1287\n  [5] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:835\n  [6] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782\n  [7] handle_deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:842\n  [8] deserialize\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Serialization/src/Serialization.jl:782 [inlined]\n  [9] deserialize_msg\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/messages.jl:87\n [10] #invokelatest#2\n    @ ./essentials.jl:708 [inlined]\n [11] invokelatest\n    @ ./essentials.jl:706 [inlined]\n [12] message_handler_loop\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/process_messages.jl:169\n [13] process_tcp_streams\n    @ /buildworker/worker/package_linux64/build/usr/share/julia/stdlib/v1.6\n/Distributed/src/process_messages.jl:126\n [14] #99\n    @ ./task.jl:411","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"two point function total plot:","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"plot(v1,outv1.-outv2)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Error: UndefVarError: outv1 not defined","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"two point function without type1-delta plot:","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"plot(v1,outv1)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Error: UndefVarError: outv1 not defined","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"two function from type1-delta plot:","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"plot(v1,outv2)","category":"page"},{"location":"markdown/ImPartTwopointfunction/","page":"-","title":"-","text":"Error: UndefVarError: outv2 not defined","category":"page"},{"location":"markdown/Check/#Check-F1-and-F2","page":"Check","title":"Check F1 & F2","text":"","category":"section"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"load package","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"using FRGRealTime,Plots,HCubature","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"Define integrate tildepartial_kF1 function:","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpp(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpp(p0,ps,x,m,T),0.0,k,initdiv=1000)[1]+FRGRealTime.loopfunpp(p0,ps,0.0,m,T)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpp (generic function with 1 method)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"Define integrate tildepartial_kF2 function:","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpm(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpm(p0,ps,x,m,T),0.0,k)[1]+FRGRealTime.loopfunpm(p0,ps,0.0,m,T)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpm (generic function with 1 method)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"Define integrate tildepartial_kF1 function for m<0:","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpp2(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpp(p0,ps,x,m,T),2.0,k,rtol=1e-10,atol=1e-10,initdiv=1000)[1]+FRGRealTime.loopfunpp(p0,ps,2.0,m,T)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpp2 (generic function with 1 method)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"Define integrate tildepartial_kF2 function:","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpm2(p0,ps,k,m,T)=hquadrature(x->FRGRealTime.flowpm(p0,ps,x,m,T),2.0,k,initdiv=1000)[1]+FRGRealTime.loopfunpm(p0,ps,2.0,m,T)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"testloopfunpm2 (generic function with 1 method)","category":"page"},{"location":"markdown/Check/#F1-Checking-Results","page":"Check","title":"F1 Checking Results","text":"","category":"section"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kpsm0p_0in (104000)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp(p0,10.0,20.0,1.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,1.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kpsm0p_0in (10100)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp(p0,10.0,20.0,1.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,1.0,145.0),1.0,20.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"ps2kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp(p0,10.0,8.0,1.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,8.0,1.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kps2m0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp(p0,10.0,2.0,1.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,2.0,1.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kpsm0p_0in (104000)","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp2(p0,10.0,20.0,-3.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,20.0,-3.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"ps2kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp2(p0,10.0,8.0,-3.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,8.0,-3.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kps2m0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp2(p0,10.0,2.0,-3.0,145.0),1.0,400.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,2.0,-3.0,145.0),1.0,400.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kps2m0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpp2(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpp(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label=\"F1\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/#F2-Checking-Results","page":"Check","title":"F2 Checking Results","text":"","category":"section"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm(p0,10.0,20.0,1.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,20.0,1.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"ps2kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm(p0,10.0,8.0,1.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,8.0,1.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kps2m0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm(p0,10.0,2.0,1.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,2.0,1.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm2(p0,10.0,20.0,-3.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,20.0,-3.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"ps2kpsm0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm2(p0,10.0,8.0,-3.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,8.0,-3.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"kps2m0","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"plot(p0->testloopfunpm2(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label=\"Integrate flow\",xaxis=\"p0\")\nplot!(p0->FRGRealTime.loopfunpm(p0,10.0,2.0,-3.0,145.0),1.0,20.0,label=\"F2\")","category":"page"},{"location":"markdown/Check/","page":"Check","title":"Check","text":"(Image: )","category":"page"},{"location":"markdown/Example/#Theory-Preset","page":"Example","title":"Theory Preset","text":"","category":"section"},{"location":"markdown/Example/#Two-point-function-flow","page":"Example","title":"Two point function flow","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"The flow of two point function is:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"beginaligned\npartial_tau Gamma_k phi_q phi_c^(2)left(p_0p_sright) \n=-frac124 frac(N+2)(2 pi)^2left-fraccothleft(fracE_pi k(k)2 Tright)left(E_pi k(k)right)^3-fracoperatornamecsch^2left(fracE_pi k(k)2 Tright)2 Tleft(E_pi k(k)right)^2right \ntimesleft(2 k^2right) int_0^k dq  q_s^2 int_-1^1 d cos theta\nendaligned","category":"page"},{"location":"markdown/Example/#Load-Package","page":"Example","title":"Load Package","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"using FRGRealTime, Plots","category":"page"},{"location":"markdown/Example/#F_1-and-F_2-at-zero-momentum","page":"Example","title":"F_1 & F_2 at zero momentum","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"F_1 will goes to zero when the width of delta function is narrow.","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(δ->FRGRealTime.loopfunppfix(1e-6,1e-5,1.0,0.5,2.0,δ),1e-10,1e-2)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"At zero momentum, F2 will goes to zero when p_0 is far small than p_s, or p_0p_s.","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"p_s=1e-5 case:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(p0->FRGRealTime.loopfunpm(p0,1e-5,1.0,0.5,2.0),1e-8,1.5*1e-5)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"p_s=1e-10 case:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(p0->FRGRealTime.loopfunpm(p0,1e-10,1.0,0.5,2.0),1e-12,1.5*1e-10)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FRGRealTime","category":"page"},{"location":"#FRGRealTime","page":"Home","title":"FRGRealTime","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FRGRealTime.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FRGRealTime]","category":"page"},{"location":"#FRGRealTime.Epi-Tuple{Any, Any}","page":"Home","title":"FRGRealTime.Epi","text":"Epi(k,m)\n\ncompute sqrt(x^2+m)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.VImintqs-NTuple{9, Any}","page":"Home","title":"FRGRealTime.VImintqs","text":"VImintqs(p0, ps, k, T, Npi,UVScale,mfun::Function,lamfun::Function)\n\ncompute int_0^kdq_s qs^2int_-1^1dcostheta mathrmImV(q_0k). In our code, we perform integration over kprim, q0 & qs does not involved, so qs=k, q0=Epi(k, mfun(k)).\n\nArguments\n\nmfun::Function: m^2(k), input from zero momentum result\nlampifun::Function: lambda_4pi(k), input from zero momentum result.\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.deltafun","page":"Home","title":"FRGRealTime.deltafun","text":"deltafun(x,dϵ=0.02)\n\ncompute delta function with the approximation fracepsilonpi(epsilon^2+x^2) and epsilon is set to be 002 by default\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.dkF1Allintqs-NTuple{6, Any}","page":"Home","title":"FRGRealTime.dkF1Allintqs","text":"dkF1Allintqs(p0, ps, qsmax, k, m, T)\n\nF_1 is an odd function about p_0\n\nwhen p0>0,dkF1Allintqs=flowpp_intcostheqs\n\nwhen p0<0,dkF1Allintqs=-flowpp_intcostheqs(-p0)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.dkF2Allintqs-NTuple{6, Any}","page":"Home","title":"FRGRealTime.dkF2Allintqs","text":"flowpm(p0, ps, k, m, T)\n\nF_2 is an odd function about p_0\n\nwhen p0>0,dkF2Allintqs=flowpm_intcostheqs\n\nwhen p0<0,dkF2Allintqs=-flowpm_intcostheqs(-p0)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.dkVImintqs-NTuple{9, Any}","page":"Home","title":"FRGRealTime.dkVImintqs","text":"dkVImintqs(p0, ps, q0, qsmax, k, m, T, Npi, lam4pik)\n\ncompute int_0^qsmaxdq_s qs^2int_-1^1dcostheta tildepartial_kmathrmImV(q_0).\n\ndkV4piImintqs only contains V(q_0), for -q_0, we have int dcostheta V(q_0)=int dcostheta V(-q_0), so we need an extra 2 at somewhere.\n\nArguments\n\nqsmax: we integrate q_s from 0 to k, qsmax will set to k when we do the integration dk, it should be distinguished from k\nm: mass square, it will be m(k) when we do the integration dk.\nlam4pik: lambda_4pi, it will be lambda_4pi(k) when we do the integration dk .\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.flowpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.flowpm","text":"flowpm(p0, ps, k, m, T)\n\ncompute -frac1pitildepartial_kIm I_2 k(p)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.flowpm_intcostheqs-NTuple{6, Any}","page":"Home","title":"FRGRealTime.flowpm_intcostheqs","text":"flowpm_intcostheqs(p0, ps, qsmax, k, m, T)\n\ncompute\n\nbeginaligned\nint_0^qsmaxdq_s q_s^2int_-1^1dcostheta tildepartial_kF_2left(sqrtp_s^2+q_s^2+2p_sq_scosthetaright)\r\n=int_0^qsmaxdq_s q_s^2int_-1^1dcostheta tildepartial_kF_2left(sqrtp_s^2+q_s^2-2*p_s*q_scosthetaright)\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.flowpp","page":"Home","title":"FRGRealTime.flowpp","text":"flowpp(p0, ps, k, m, T,δ=0.02)\n\ncompute -frac1pitildepartial_kIm I_1 k(p)\n\nTo be noticed that, flowpp doesn't constains tildepartial_kmathcalF_4, we will consider it seperately.\n\nAt p_0=2E_pik, flowpp has a delta function contribution, we use a rectangle function with width delta to approximate.\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.flowpp_intcostheqs-NTuple{6, Any}","page":"Home","title":"FRGRealTime.flowpp_intcostheqs","text":"flowpp_intcostheqs(p0, ps, qsmax, k, m, T)\n\ncompute\n\nbeginaligned\nint_0^qsmaxdq_s q_s^2int_-1^1dcostheta tildepartial_kF_1left(sqrtp_s^2+q_s^2+2*p_s*q_scosthetaright)\r\n=int_0^qsmaxdq_s q_s^2int_-1^1dcostheta tildepartial_kF_1left(sqrtp_s^2+q_s^2-2*p_s*q_scosthetaright)\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.loopfunpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.loopfunpm","text":"loopfunpm(p0, ps, k, m, T)\n\ncompute -frac1piIm I_2 k(p)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.loopfunpp-NTuple{5, Any}","page":"Home","title":"FRGRealTime.loopfunpp","text":"loopfunpp(p0, ps, k, m, T)\n\ncompute -frac1piIm I_2 k(p), but without mathcalF_4 contribution\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.loopfunppfix","page":"Home","title":"FRGRealTime.loopfunppfix","text":"loopfunppfix(p0, ps, k, m, T)\n\ncompute -frac1piIm I_1 k(p)\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.star1fun-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star1fun","text":"star1fun(qp, qm, ps, m, T)\n\ncompute -frac1pimathcalF_1left(q_+ q_- p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star1funpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star1funpm","text":"star1funpm(qp, qm, ps, m, T)\n\ncompute -frac1pimathcalF_1^primeleft(q_+ q_- p barm_pi k^2right)\n\nCompared with mathcalF, the mathcalF has subtracted vacuum contributions. it will be used in mathrmIm I_2\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star2fun-NTuple{6, Any}","page":"Home","title":"FRGRealTime.star2fun","text":"star2fun(qp, qm, ps, k, m, T)\n\ncompute -frac1pimathcalF_2left(q_+ q_- k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star2funpm-NTuple{6, Any}","page":"Home","title":"FRGRealTime.star2funpm","text":"star2funpm(qp, qm, ps, k, m, T)\n\ncompute -frac1pimathcalF^prime_2left(q_+ q_- k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star3fun-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star3fun","text":"star3fun(qp, ps, k, m, T)\n\ncompute -frac1pimathcalF_3left(q_+k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star3funpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star3funpm","text":"star3funpm(qp, ps, k, m, T)\n\ncompute -frac1pimathcalF^prime_3left(q_+ k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star4fun","page":"Home","title":"FRGRealTime.star4fun","text":"star4fun(p0, p, k, m, T)\n\ncompute -frac1pimathcalF_4left(p_0 k p barm_pi k^2right)\n\n\n\n\n\n","category":"function"},{"location":"markdown/NumericalTricks/#Delta-Function","page":"Numerical Tricks","title":"Delta Function","text":"","category":"section"},{"location":"markdown/NumericalTricks/#Integrate-Dirac-Delta-in-F_1","page":"Numerical Tricks","title":"Integrate Dirac Delta in F_1","text":"","category":"section"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"There are two types of dirac delta function we need to handle. One is due to the existence of thetaleft(p_0-2E_piright), the other is due to the existence of deltaleft(p_0-2E_piright).","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"For the first one, we can see the plot around p_0=2E_pi","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"using FRGRealTime,Plots\n#k>ps/2\nplot(p0->FRGRealTime.loopfunppfix(p0,1.0,2.0,1.0,2.0),2*Epi(2.0,1.0)-0.01,2*Epi(2.0,1.0)+0.01)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"(Image: )","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"Considering that our function has the following form:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    F(k)=f(k)thetaBig(g(k)-p_0Big)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"and this relation is only true when p_0 is in a infinitesimal neighborhood around g(k). F(k) is:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    F(k)=f(k) g(k) delta Big(g(k)-p_0Big)+f(k) theta Big(g(k)-p_0Big)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"Now we want to integrate back to F(k), any finite part will disappear due to the infinitesimal measure, so:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    F(k)=f(k) g(k) delta Big(g(k)-p_0Big)=f(k)partial_k thetaBig(g(k)-p_0Big)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"F_1 is theta function when kfracp_s2 and p_0 is around 2E_pi(km). It has a drop=fracleft(2 k-p_sright) coth left(fracsqrtk^2+m^22 Tright)16 pi ^2 sqrtk^2+m^2 drop in p_0=2E_pi(km). F_1 can be written as:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    F_1=dropthetaBig(p_0-2E_pi(km)Big)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"when p_0=2E_pi(km) :","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    partial_k F_1=-dropfrac2kE_pi(km)deltaBig(p_0-2E_pi(km)Big)","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"alternatively:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"    partial_k F_1=-drop deltaleft(k-frac12 sqrtp_0^2-4 m^2right)","category":"page"},{"location":"markdown/NumericalTricks/#Im-Part","page":"Numerical Tricks","title":"Im Part","text":"","category":"section"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"The only chance that we can get rid of these annoying delta function in the Im part is doing the k or k integral:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"mathrmdeltaintkqsequivint_Lambda^kdkint_0^qsmaxdq_sq_s^2int_-1^1dcosthetatildepartial_kF_1","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"the q_s and costheta integration will act on the drop, the k integration will act on the delta function.","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"We also noticed that this delta function only appears in kp2, so the internel momentum integration can't performed on the whole range. Assuming that the solution of  p_0=2E_pi  is k_0. So we only do the k integration around k_0, then p2k_0.","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"Integrate[((2 k - Sqrt[ps^2 + qs^2 + 2*ps*qs*x])*\n   Coth[Sqrt[k^2 + m]/(2 T)])/(16 \\[Pi]^2 Sqrt[k^2 + m]) qs^2, {x, -1,\n   1}, Assumptions -> qsmax > 0 && 1 > x > -1 && ps > 0 && qs > 0]","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"we get:","category":"page"},{"location":"markdown/NumericalTricks/","page":"Numerical Tricks","title":"Numerical Tricks","text":"beginaligned\np_sk_0 int_0^k_0\np_s3k_0 0\np_s3k_0 0\np_s3k_0 0\n\nendaligned","category":"page"}]
}
