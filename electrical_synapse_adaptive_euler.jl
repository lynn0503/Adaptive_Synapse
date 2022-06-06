function electrical_synapse_adaptive_euler(trial_number,tspan)
    # C = 0.5
    g_L = 0.025 # leaky conductance
    V_rest = -70
    V_th = -50
    V_peak = 40
    d_es = 0.3 # time dedaly
    τ_f = 1000 # time constant of short-term facilitation
    τ_l = 50 # time window of plasticity
    τ_u = 10 # τ_u = RC
    u_f = 0.05 #rate of facilitation
    g_0 =0.022
    g_max =0.1
    σ = 5 # 
    c = 0.5 # common input current
    trials =  trial_number
    tspan = tspan 
    τ_ref = 2 # refractory period 2ms
    dt = 1 # step size 1ms

    # init global 
    u1 = Array{Float64}(undef,trials,tspan)
    u2 = Array{Float64}(undef,trials,tspan)
    ges = Array{Float64}(undef,trials,tspan)
    spike_train1 = fill(0,trials,tspan)
    spike_train2 = fill(0,trials,tspan)

    g_temp = g_0
    # step in all trials
    for k in 1:trials
        if k==1
            ges[k,1]=g_temp
        else
            ges[k,1]=ges[k-1,tspan]
            g_temp = ges[k-1,tspan]
        end
        # init local
        u1[k,1]= V_rest + randn()
        u2[k,1]= V_rest + randn()
        neu_t1=0
        spike_t1=0
        neu_t2=0
        spike_t2=0
        t_ref1 = τ_ref/dt+1
        t_ref2 = τ_ref/dt+1
        I_c = randn() #common input 
        I_input = [20*exp(-i/20)+σ*(sqrt(c)*I_c+sqrt(1-c)*randn()) for i in 1:tspan]
        I_input2 = [20*exp(-i/20)+σ*(sqrt(c)*I_c+sqrt(1-c)*randn()) for i in 1:tspan]

        # step in one trial
        for i in 2:tspan
            ges[k,i] = g_temp
            # u1
            if t_ref1 > τ_ref/dt
                if u2[k,i-1]>u1[k,i-1]
                    delta_u1 = - (u1[k,i-1] - V_rest)/τ_u + I_input[i-1]/(τ_u*g_L) + (u2[k,i-1]-u1[k,i-1])*g_temp/(τ_u*g_L)
                else
                    delta_u1 = - (u1[k,i-1] - V_rest)/τ_u + I_input[i-1]/(τ_u*g_L) 
                end
                u1[k,i] = u1[k,i-1] + delta_u1
            else
                delta_u1 = 0
                u1[k,i] = V_rest
                t_ref1 = t_ref1 + 1
            end
            # u2
            if t_ref2 > τ_ref/dt
                if u1[k,i-1]>u2[k,i-1]
                    delta_u2 = - (u2[k,i-1] - V_rest)/τ_u + I_input2[i-1]/(τ_u*g_L) + (u1[k,i-1]-u2[k,i-1])*g_temp/(τ_u*g_L)
                else
                    delta_u2 = - (u2[k,i-1] - V_rest)/τ_u + I_input2[i-1]/(τ_u*g_L) 
                end
                u2[k,i] = u2[k,i-1] + delta_u2
            else
                delta_u2 = 0
                u2[k,i] = V_rest
                t_ref2 = t_ref2 + 1
            end
            # ges and u1
            if u1[k,i] >= V_th
                # fire leads to ltp
                u1[k,i] = V_peak
                spike_train1[k,i] = 1
                neu_t1 = i 
                spike_t1 = i 
                t_ref1 = 1 
                if spike_t2 > 0
                    delta_ges = (-g_temp + g_0)/τ_f + u_f/τ_f*(g_max-g_temp)*exp(-abs(i-spike_t2)*dt/τ_l)
                    g_temp = g_temp + delta_ges
                end
            else
                # no fire leads to ltd
                delta_ges = (-g_temp + g_0)/τ_f
                g_temp = g_temp + delta_ges

            end
            # ges and u2
            if u2[k,i] >= V_th
                # fire leads to ltp
                u2[k,i] = V_peak
                spike_train2[k,i] = 1
                neu_t2 = i 
                spike_t2 = i 
                t_ref2 = 1 
                if spike_t1 > 0
                    delta_ges = (-g_temp + g_0)/τ_f + u_f/τ_f*(g_max-g_temp)*exp(-abs(i-spike_t1)*dt/τ_l)
                    g_temp = g_temp + delta_ges
                end
            else
                # no fire leads to ltd
                delta_ges = (-g_temp + g_0)/τ_f
                g_temp = g_temp + delta_ges

            end
        end
    end
    return u1, u2, ges, spike_train1, spike_train2
end