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

function ccf_sampling(u1,u2,binsize,corr_len,samples)
    # u1=[ui[1] for ui in sol.u]
    # u2=[ui[2] for ui in sol.u]
    # u1=u1[1,:]
    # u2=u2[1,:]
    
    idx=1:length(u1)
    fire1_idx=idx[u1.>0]
    s1,bins=histcountindices(fire1_idx,0:binsize:length(u1))
    fire2_idx=idx[u2.>0]
    s2,bins=histcountindices(fire2_idx,0:binsize:length(u1))
    # samples=10
    N=length(s1)
    n=N/samples
    ccf_by_time=Array{Float64}(undef,samples,1)
    for i in 1:samples-1
        start_idx=floor(Int,(i-1)*n+1)
        s1_temp=s1[start_idx:start_idx+corr_len-1]
        s2_temp=s2[start_idx:start_idx+corr_len-1]
        corr=crosscor(s1_temp,s2_temp)
        ccf_by_time[i]=maximum(abs.(corr))
    end
    return ccf_by_time

end
