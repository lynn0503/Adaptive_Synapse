
include("electrical_synapse.jl")
include("chemical_synapse.jl")
include("ccf.jl")

using Statistics
using Plots
using DifferentialEquations

gr()
l = @layout [a b]

cs=0:0.1:1
tspan = (0.0,200.0)
trials=10

# electrical synapse

corrs=Array{Float64}(undef,4,length(cs),trials)
u0=[-60,-60]
h(p,t)=[-60,-60]

for (i,g_es) in enumerate([0.025, 0.05, 0.1, 0.2])
    for (j,c) in enumerate(cs)
        for k in 1:trials
            @show "electrical", g_es, c, k
            p=[0.5, 0.025, -70, -50, g_es, 0.3, 0.62, 0.5, c]
            # C, g_L, V_rest,V_th, g_es, d_es, μ, σ, c = p
            prob = DDEProblem(electric_synapse!,u0,h,tspan,p)
            sol = solve(prob,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)
            corrs[i,j,k]=maximum(ccf(sol,2))
        end
    end
end

corrs_mean=mean(corrs,dims=3)[:,:,1]

p1=plot(
    cs,
    corrs_mean',
    ylims=(0,1),
    xlabel="Input Correlation",
    ylabel="Correlation Strength",
    label=["g_es=0.001" "g_es=0.025" "g_es=0.05" "g_es=0.1"],
    legend=:bottomright,
    lw=3
    )


# chemical synapse
u0c=[-60;-60;0.025]
hc(p,t)=[-60,-60,0.025]
corrsc=Array{Float64}(undef,4,length(cs),trials)
tspan = (0.0,100.0)
for (i,uu) in enumerate([0,0.05,0.08,0.1])
    for (j,c) in enumerate(cs)
        for k in 1:trials
            @show "chemical", uu, c, k
            pc=[0.5, 0.025, -70, -50, 3, 0, 5, uu, 0.62, 0.5, c]
            # C, g_L, V_rest,V_th, d_cs, V_rev, τ_s, uu, μ, σ, c = p
            probc = DDEProblem(chemical_synapse!,u0c,hc,tspan,pc)
            solc = solve(probc,dt=0.1,dtmax=0.1, dtmin=0.1,force_dtmin=true)
            corrsc[i,j,k]=maximum(ccf(solc,2))
        end
    end
end
corrsc_mean=mean(corrsc,dims=3)[:,:,1]

p2=plot(
    cs,
    corrsc_mean',
    ylims=(0,1),
    xlabel="Input Correlation",
    ylabel="Correlation Strength",
    label=["u=0" "u=0.05" "u=0.08" "u=0.1"],
    legend=:bottomright,
    lw=3
    )

# plot all
plot(p1, p2, layout = l)