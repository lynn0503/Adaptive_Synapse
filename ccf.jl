using NaNStatistics
using StatsBase

function ccf(sol,binsize)
    u1=[ui[1] for ui in sol.u]
    u2=[ui[2] for ui in sol.u]
    n=length(u1)
#     @show n
    idx=1:n
    fire1_idx=idx[u1.>0]
    fire1_binned,bins=histcountindices(fire1_idx,0:binsize/0.1:n)
    fire2_idx=idx[u2.>0]
    fire2_binned,bins=histcountindices(fire2_idx,0:binsize/0.1:n)
    # @show sum(fire1_binned),sum(fire2_binned)
    N=length(fire1_binned)
    @show n,N
    # ccf=xcorr(fire1_binned,fire2_binned)
    # positive t
    # ccf=[N*sum(fire1_binned[t+1:N-t].*fire2_binned[1+2*t:N])/((N-t)*sqrt(sum(fire1_binned)*sum(fire2_binned))) for t in 0:N-1]
    # # negative t
    # ccf1=[N*sum(fire1_binned[1-t:N+t].*fire2_binned[1:N+2*t])/((N+t)*sqrt(sum(fire1_binned)*sum(fire2_binned))) for t in -100:0]
    # ccf=vcat(ccf1,ccf2)
    # all 
    # ccf=[N*sum(fire1_binned[abs(t)+1:N-abs(t)].*fire2_binned[1+2*abs(t):N])/((N-2*abs(t))*sqrt(sum(fire1_binned)*sum(fire2_binned))) for t in 1-N:N-1]
    ccf=abs.(crosscor(fire1_binned,fire2_binned))
    return ccf
    
end