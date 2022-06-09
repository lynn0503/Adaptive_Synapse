include("electrical_synapse_adaptive_euler.jl")
include("spike_train2firing_rate.jl")
include("ccf.jl")

using Plots
using Statistics

gr()
l = @layout [a b; c d]
tspan = 10000
trials=10
u1, u2, ges, st1, st2 = electrical_synapse_adaptive_euler(trials,tspan)
# ges=mean(ges,dims=1)
# u1=mean(u1,dims=1)
# u2=mean(u2,dims=1)
# st1=mean(st1,dims=1)
# st2=mean(st2,dims=1)
trial=8
a=plot(ges[trial,:],ylabel="g_es")
b=bar_fire_rate(st1[trial,:],50)
c=bar_fire_rate(st2[trial,:],25)

# c=plot(u1_mean',label="u1",xlims=(0,1000))
# d=plot(u2_mean',label="u2",xlims=(0,1000))
binsize=10
corr_len=100
samples= floor(Int,tspan/(binsize*corr_len))
ccf_by_time=ccf_sampling(u1[trial,:],u2[trial,:],binsize,corr_len,samples)
d=plot(ccf_by_time[1:5],ylabel="correlation strength")

plot(a,b,c,d,layout=l)


