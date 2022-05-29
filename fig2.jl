# constant g_es and uu 
include("electrical_synapse.jl")
include("chemical_synapse.jl")
include("ccf.jl")

using Statistics
using Plots

gr()
l = @layout [a b]

cs=0:0.1:1
tspan = (0.0,100.0)

# electrical synapse
corrs=Array{Float64}(undef,length(cs),10)
u0=[-60,-60]
h(p,t)=[-60,-60]

for (idx,c) in enumerate(cs)
    for i in 1:10
        @show "electrical c=", c, "i=",i
        p=[0.5, 0.025, -70, -50, 0.025, 0.3, 0.62, 0.5, c]
        # C, g_L, V_rest,V_th, g_es, d_es, μ, σ, c = p
        prob = DDEProblem(electric_synapse!,u0,h,tspan,p)
        sol = solve(prob,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)
        corrs[idx,i]=maximum(ccf(sol,2))
    end
end

p1=plot(mean(corrs,dims=2))

# chemical synapse
u0c=[-60;-60;0.025]
hc(p,t)=[-60,-60,0.025]
corrsc=Array{Float64}(undef,length(cs),10)

for (idx,c) in enumerate(cs)
    for i in 1:10
        @show "chemical c=", c,  i
        pc=[0.5, 0.025, -70, -50, 3, 0, 5, 0.05, 0.62, 0.5, c]
        # C, g_L, V_rest,V_th, d_cs, V_rev, τ_s, uu, μ, σ, c = p
        probc = DDEProblem(chemical_synapse!,u0c,hc,tspan,pc)
        solc = solve(probc,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)
        corrsc[idx,i]=maximum(ccf(solc,2))
    end
end

p2=plot(mean(corrsc,dims=2))

# plot all
plot(p1, p2, layout = l)