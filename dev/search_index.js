var documenterSearchIndex = {"docs":
[{"location":"markdown/Example/#Theory-Preset","page":"Example","title":"Theory Preset","text":"","category":"section"},{"location":"markdown/Example/#Two-point-function-flow","page":"Example","title":"Two point function flow","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"The flow of two point function is:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"beginaligned\npartial_tau Gamma_k phi_q phi_c^(2)left(p_0p_sright) \n=-frac124 frac(N+2)(2 pi)^2left-fraccothleft(fracE_pi k(k)2 Tright)left(E_pi k(k)right)^3-fracoperatornamecsch^2left(fracE_pi k(k)2 Tright)2 Tleft(E_pi k(k)right)^2right \ntimesleft(2 k^2right) int_0^k dq  q_s^2 int_-1^1 d cos theta\nendaligned","category":"page"},{"location":"markdown/Example/#Load-Package","page":"Example","title":"Load Package","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"using FRGRealTime, Plots","category":"page"},{"location":"markdown/Example/#F_1-and-F_2-at-zero-momentum","page":"Example","title":"F_1 & F_2 at zero momentum","text":"","category":"section"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"F_1 will goes to zero when the width of delta function is narrow.","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(δ->FRGRealTime.loopfunppfix(1e-6,1e-5,1.0,0.5,2.0,δ),1e-10,1e-2,dpi=250)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"At zero momentum, F2 will goes to zero when p_0 is far small than p_s, or p_0p_s.","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"p_s=1e-5","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"case:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(p0->FRGRealTime.loopfunpm(p0,1e-5,1.0,0.5,2.0),1e-8,1.5*1e-5,dpi=250)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"p_s=1e-10","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"case:","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"plot(p0->FRGRealTime.loopfunpm(p0,1e-10,1.0,0.5,2.0),1e-12,1.5*1e-10,dpi=250)","category":"page"},{"location":"markdown/Example/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FRGRealTime","category":"page"},{"location":"#FRGRealTime","page":"Home","title":"FRGRealTime","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FRGRealTime.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FRGRealTime]","category":"page"},{"location":"#FRGRealTime.Epi-Tuple{Any, Any}","page":"Home","title":"FRGRealTime.Epi","text":"Epi(k,m)\n\ncompute sqrt(x^2+m)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.deltafun","page":"Home","title":"FRGRealTime.deltafun","text":"deltafun(x,dϵ=0.02)\n\ncompute delta function with the approximation fracepsilonpi(epsilon^2+x^2) and epsilon is set to be 002 by default\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.flowpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.flowpm","text":"flowpm(p0, ps, k, m, T)\n\ncompute -frac1pitildepartial_kIm I_2 k(p)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.flowpp","page":"Home","title":"FRGRealTime.flowpp","text":"flowpp(p0, ps, k, m, T,δ=0.02)\n\ncompute -frac1pitildepartial_kIm I_1 k(p)\n\nTo be noticed that, flowpp doesn't constains tildepartial_kmathcalF_4, we will consider it seperately.\n\nAt p_0=2E_pik, flowpp has a delta function contribution, we use a rectangle function with width delta to approximate.\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.loopfunpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.loopfunpm","text":"loopfunpm(p0, ps, k, m, T)\n\ncompute -frac1piIm I_2 k(p)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.loopfunpp-NTuple{5, Any}","page":"Home","title":"FRGRealTime.loopfunpp","text":"loopfunpp(p0, ps, k, m, T)\n\ncompute -frac1piIm I_2 k(p), but without mathcalF_4 contribution\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.loopfunppfix","page":"Home","title":"FRGRealTime.loopfunppfix","text":"loopfunppfix(p0, ps, k, m, T)\n\ncompute -frac1piIm I_1 k(p)\n\n\n\n\n\n","category":"function"},{"location":"#FRGRealTime.star1fun-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star1fun","text":"star1fun(qp, qm, ps, m, T)\n\ncompute -frac1pimathcalF_1left(q_+ q_- p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star1funpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star1funpm","text":"star1funpm(qp, qm, ps, m, T)\n\ncompute -frac1pimathcalF_1^primeleft(q_+ q_- p barm_pi k^2right)\n\nCompared with mathcalF, the mathcalF has subtracted vacuum contributions. it will be used in mathrmIm I_2\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star2fun-NTuple{6, Any}","page":"Home","title":"FRGRealTime.star2fun","text":"star2fun(qp, qm, ps, k, m, T)\n\ncompute -frac1pimathcalF_2left(q_+ q_- k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star2funpm-NTuple{6, Any}","page":"Home","title":"FRGRealTime.star2funpm","text":"star2funpm(qp, qm, ps, k, m, T)\n\ncompute -frac1pimathcalF^prime_2left(q_+ q_- k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star3fun-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star3fun","text":"star3fun(qp, ps, k, m, T)\n\ncompute -frac1pimathcalF_3left(q_+k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star3funpm-NTuple{5, Any}","page":"Home","title":"FRGRealTime.star3funpm","text":"star3funpm(qp, ps, k, m, T)\n\ncompute -frac1pimathcalF^prime_3left(q_+ k p barm_pi k^2right)\n\n\n\n\n\n","category":"method"},{"location":"#FRGRealTime.star4fun","page":"Home","title":"FRGRealTime.star4fun","text":"star4fun(p0, p, k, m, T)\n\ncompute -frac1pimathcalF_4left(p_0 k p barm_pi k^2right)\n\n\n\n\n\n","category":"function"}]
}
