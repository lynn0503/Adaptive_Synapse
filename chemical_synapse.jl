function chemical_synapse!(du,u,h,p,t)
    # parameter
    C, g_L, V_rest,V_th, d_cs, V_rev, τ_s, uu, μ, σ, c = p
    
    # if spiked, reset
    if u[1]>0
        u[1]=V_rest
    end
    if u[2]>0
        u[2]=V_rest
    end
    
    # dynamical equation
    I_ext=μ+σ*(sqrt(c)*randn())
    I_ext1=σ*(sqrt(1-c)*randn())
    I_ext2=σ*(sqrt(1-c)*randn())
    du[1]=1/C*(-g_L*(u[1]-V_rest)+I_ext+I_ext1)
    du[2]=1/C*(-g_L*(u[2]-V_rest)-h(p,t-d_cs)[3]*(u[2]-V_rev)+I_ext+I_ext2)
    du[3]=1/τ_s*(-u[3]+uu*(h(p,t)[1]>0))
    # u[3] is g_21, connection from neuron1 to neuron2
    
    # if above threshold, then spike
    if u[1]>V_th
        u[1]=30
    end
    if u[2]>V_th
        u[2]=30
    end
    
end

