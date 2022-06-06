include("electrical_synapse.jl")
include("chemical_synapse.jl")
include("ccf.jl")

using Plots
using DifferentialEquations
gr(leg=false)
l = @layout [a b; c d]

tspan = (0,200)
# electrical synapse
u0=[-60;-60]
h(p,t)=[-60,-60]
p=[0.5, 0.025, -70, -50, 0.1, 0.3, 0.62, 0.5, 0.5]
# C, g_L, V_rest,V_th, g_es, d_es, μ, σ, c = p
prob = DDEProblem(electric_synapse!,u0,h,tspan,p)
# sol = solve(prob)
sol = solve(prob,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)

p1=plot(sol,ylims=(-80,50))
p2=bar(ccf(sol,2),ylims=(0,1))


# chemical synapse
u0c=[-60;-60;0.025]
hc(p,t)=[-60,-60,0.025]
pc=[0.5, 0.025, -70, -50, 3, 0, 5, 0.05, 0.62, 0.5, 0.5]
# C, g_L, V_rest,V_th, d_cs, V_rev, τ_s, uu, μ, σ, c = p
probc = DDEProblem(chemical_synapse!,u0c,hc,tspan,pc)
# solc = solve(probc)
solc = solve(probc,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)

p3=plot(solc,ylims=(-80,50))
p4=bar(ccf(solc,2),ylims=(0,1))

plot(p1, p2, p3, p4, layout = l)


# u2=[ui[3] for ui in solc.u]
# plot(u2)
# xlims!(0,1000)