ϵ5=0.05
lpakrang = RGscale(800.0, 1.0)
const q0grid = collect(0.0:0.2:lpakrang.UV)
Tlpa = 145.0*lpakrang.UV/800
Tself=Tlpa-8.07
mass2_UV=-0.18216546875*lpakrang.UV^2
lambda_UV=6*8.0

p0min_Tc=1e-4
p0max_Tc=1.0
N_p0grid_Tc=400
p0step_Tc = (log(p0max_Tc) - log(p0min_Tc)) / (N_p0grid_Tc - 1)
p0grid_Tc = SharedArray(exp.(collect(log(p0min_Tc):p0step_Tc:log(p0max_Tc))))
