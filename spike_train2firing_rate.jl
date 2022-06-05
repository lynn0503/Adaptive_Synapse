function spike_train2firing_rate(sol,binsize)
    # unit of binsize is ms
    # stepsize is 0.1ms by default
    u1=[ui[1] for ui in sol.u]
    u2=[ui[2] for ui in sol.u]

    n=length(u1)
#     @show n
    idx=1:n
    fire1_idx=idx[u1.>0]
    fire1_binned,bins=histcountindices(fire1_idx,0:binsize/0.1:n)
    fire2_idx=idx[u2.>0]
    fire2_binned,bins=histcountindices(fire2_idx,0:binsize/0.1:n)
    
    return fire1_binned,fire2_binned
    
end