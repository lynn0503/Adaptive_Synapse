include("electrical_synapse_adaptive.jl")
include("spike_train2firing_rate.jl")
include("ccf.jl")

using DifferentialEquations
using Plots

gr()
l = @layout [a b;c d]
# fig a and b are firing rate of neuron1 and neuron2 respectively
# fig c is the correlation of neuron1 and neuron2
# fid d is g_es change with time
tspan = (0,10000)
# 100s
# electrical synapse adaptive
u0a=[-60;-60;0.025]
ha(p,t)=[-60,-60,0.025]
pa=[0.5, 0.025, -70, -50, 0.3, 1000, 5, 0.5, 0.022, 0.1, 5, 0.5]
#  C, g_L, V_rest,V_th, d_es, τ_f, τ_l, u_f, g_0, g_max, σ, c = p

# τ_f=10s=10000ms, τ_l=50ms, u_f=50us=0.05ms
# according to fig4, g_0=0.022, g_max=0.03
# according to code, g_0=0.002, g_max=0.1
# c is the same with fig1, which is 0 for electrical synapse
proba = DDEProblem(electrical_synapse_adaptive!,u0a,ha,tspan,pa)
# sola = solve(proba)
sola = solve(proba,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)
# p1=plot(sola)

fr1,fr2=spike_train2firing_rate(sola,100)
p3=bar(fr1,label="neuron1",ylabel="Firing Rate",xlims=(0,50))
p4=bar(fr2,label="neuron2",ylabel="Firing Rate",xlims=(0,50))

ccf_by_time=ccf_sampling(sola,10,100,20)
p1=plot(ccf_by_time,label="correlation strength")
ges=[ui[3] for ui in sola.u]
p2=plot(ges,label="g_es")

plot(p3, p4, p1, p2, layout = l)



