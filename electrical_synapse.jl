using DifferentialEquations
# using Logging
# Logging.disable_logging(Logging.Info)

function electric_synapse!(du,u,h,p,t)
    # parameter
    C, g_L, V_rest,V_th, g_es, d_es, μ, σ, c = p
    
    # if spiked, reset
    if u[1]>0
        u[1]=V_rest
    end
    if u[2]>0
        u[2]=V_rest
    end
    
    # dynamical equation
    I_ext=μ+σ*(sqrt(c)*randn())
    I_ext1=σ*sqrt(1-c)*randn()
    I_ext2=σ*sqrt(1-c)*randn()
    du[1]=1/C*(-g_L*(u[1]-V_rest)+g_es*(h(p,t-d_es)[2]-u[1])+I_ext+I_ext1)
    du[2]=1/C*(-g_L*(u[2]-V_rest)+g_es*(h(p,t-d_es)[1]-u[2])+I_ext+I_ext2)
    
    # if above threshold, spike
    if u[1]>V_th
        u[1]=30
    end
    if u[2]>V_th
        u[2]=30
    end
    
end

