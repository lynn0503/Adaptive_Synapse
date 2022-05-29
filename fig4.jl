include("electrical_synapse_adaptive.jl")
include("ccf.jl")

using Plots
gr(leg=false)
l = @layout [a b]

tspan = (0,500)
# electrical synapse adaptive
u0a=[-60;-60;0.025]
ha(p,t)=[-60,-60,0.025]
pa=[0.5, 0.025, -70, -50, 0.3, 10000, 50, 0.05, 0.022, 0.03, 20, 0.5, 0]
# C, g_L, V_rest,V_th, d_es, τ_f, τ_l, u_f, g_0, g_max, a, σ, c = p
# τ_f=10s=10000ms, τ_l=50ms, u_f=50us=0.05ms
# g_0=0.022, g_max=0.03, g_init=0.025
proba = DDEProblem(electric_synapse_adaptive!,u0a,ha,tspan,pa)
# sol = solve(prob)
sola = solve(proba,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)

plot(sola)

# u2=[ui[3] for ui in sola.u]
# plot(u2)