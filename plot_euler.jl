include("electrical_synapse_adaptive_euler.jl")
include("spike_train2firing_rate.jl")
using Plots
gr()
l = @layout [a b;c d]
tspan = 10000
u1, u2, ges, st1, st2 = electrical_synapse_adaptive_euler(20,tspan)
plot(ges[2,:])
bar_fire_rate(st1[20,:],100)