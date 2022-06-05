include("electrical_synapse_adaptive.jl")
include("spike_train2firing_rate.jl")

using Plots
gr()
l = @layout [a b; c d]
# fig a and b are firing rate of neuron1 and neuron2 respectively
# fig c is the correlation of neuron1 and neuron2
# fid d is g_es change with time
tspan = (0,500)
# electrical synapse adaptive
u0a=[-60;-60;0.025]
ha(p,t)=[-60,-60,0.025]
pa=[0.5, 0.025, -70, -50, 0.3, 10000, 50, 0.05, 0.002, 0.1, 20, 0.5, 0]
# C, g_L, V_rest,V_th, d_es, τ_f, τ_l, u_f, g_0, g_max, a, σ, c = p
# τ_f=10s=10000ms, τ_l=50ms, u_f=50us=0.05ms
# g_0=0.002, g_max=0.1 are fitted from experiment data
# c is the same with fig1, which is 0 for electrical synapse
proba = DDEProblem(electric_synapse_adaptive!,u0a,ha,tspan,pa)
# sol = solve(prob)
sola = solve(proba,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)

fr1,fr2=spike_train2firing_rate(sola,5)

a=bar(fr1)
b=bar(fr2)

ges=[ui[3] for ui in sola.u]
d=plot(ges)

plot(a, b, c, d, layout = l)

