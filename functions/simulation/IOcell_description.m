================================================]
          INFERIOR OLIVE MODEL DESCRIPTION
================================================]


 Our network model currently does the computation in two steps:
 1. Updates the cross gap junction currents
 2. Updates the three compartments

 Note: In step 1 the matlab implementation uses gpuArrays    

 gap junctions between all connected cells computed according to this function:

 step 1 : for all gaps
    fgap = @(DeltaV) (0.8 .* exp(-1.*DeltaV^2/100) + 0.2);
     Where DeltaV is the difference of potential between two cells

 Inputs are currently like this

     1 ampa conductance per cell
     1 gaba conductance per soma
     1 gaba conductance per dendrite
     1 Current input to the soma (we compute the ohrstein-ulenbeck at all dt)

================================================]
          SINGLE IO CELL
================================================]
 Single cells have three compartments:

 CURRENTS OF THE SOMA
    dVs_dt = (-(I_CaL   + I_ds  + I_as + I_Na_s + I_ls  + I_h_s   + I_Kdr_s +I_K_s  + I_amp + I_gab_soma) ) / C_m

 CURRENTS OF DENDRITE
    dVd_dt = (-(I_CaH   + I_sd  + I_ld + I_K_Ca + I_cx36 + I_h     + I_gab_dend)  + I_app ) / C_m

 CURRENTS OF THE AXON
    dVa_dt = (-(I_K_a   + I_sa  + I_la + I_Na_a)                                   ) / C_m
    V_soma =  dVs_dt + V_soma
    V_dend =  dVd_dt + V_dend
    V_axon =  dVa_dt + V_axon

================================================]
          SOMA
================================================]

 1 somatic calcium Low Threshold current

    I_CaL = g_CaL * k * k * k * l * (V_soma - V_Ca)

    k_inf = (1 / (1 + exp(-1 * (V_soma + 61)   / 4.2)))
    l_inf = (1 / (1 + exp((     V_soma + 85.5) / 8.5)))
    tau_k = 1
    tau_l = ((20 * exp((V_soma + 160) / 30) / (1 + exp((V_soma + 84) / 7.3))) +35)
    dk_dt = (k_inf - Calcium_k) / tau_k
    dl_dt = (l_inf - Calcium_l) / tau_l
        k = dk_dt + Calcium_k
        l = dl_dt + Calcium_l


 1.2 somatic Sodium current

    I_Na_s  = g_Na_s * m * m * m * h * (V_soma - V_Na)

    m_inf = 1 / (1 + (exp((-30 - V_soma)/ 5.5)))
    h_inf = 1 / (1 + (exp((-70 - V_soma)/-5.8)))
    tau_h =       3 * exp((-40 - V_soma)/33)

    dh_dt = (h_inf - Sodium_h)/tau_h
     
        m     = m_inf
        h     = Sodium_h +  dh_dt

 1.3 somatic Potassium rectifier current

    I_Kdr_s = g_Kdr_s * n * n * n * n * (V_soma - V_K)

    n_inf = 1 / (1 + exp( ( -3 - V_soma) /  10))
    tau_n =   5 + (  47 * exp( -(-50 - V_soma) /  900))
    dn_dt = (n_inf - Potassium_n) / tau_n
    n     =  dn_dt + Potassium_n

 1.4 Potassium 

    I_K_s   = g_K_s * (x_s ^ 4) * (V_soma - V_K)

    alpha_x_s = 0.13 * (V_soma + 25) / (1 - exp(-(V_soma + 25) / 10))
    beta_x_s  = 1.69 * exp(-0.0125 * (V_soma + 35))

    x_inf_s   = alpha_x_s / (alpha_x_s + beta_x_s)
    tau_x_s   =         1 / (alpha_x_s + beta_x_s)
    dx_dt_s   = (x_inf_s - Potassium_x_s) / tau_x_s

    x_s       = dx_dt_s + Potassium_x_s

================================================]
 2 Dendritic components
=================================================]

 2.1 Soma-dendrite interaction current I_sd
    I_sd   = (g_int / (1 - p1)) * (V_dend - V_soma)

    
 2.2 H current (Inward anomalous rectifier)
    I_h    =  g_h * q * (V_dend
        q_inf = 1 /(1 + exp((V_dend + 80) / 4))
        tau_q = 1 /(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87))
        dq_dt = (q_inf - Hcurrent_q) / tau_q
        q =  dq_dt + Hcurrent_q  original

    
    
    
 2.4 Outward Ca-dependent K current I_K_Ca (Potassium BK channel)
    I_K_Ca =  g_K_Ca * s * (V_dend - V_K)
    
      alpha_s = min([0.00002 * Ca2Plus, 0.01])
       alternatively: alpha_s = (0.00002 * Ca2Plus) * (0.00002 * Ca2Plus < 0.01) + 0.01*((0.00002 * Ca2Plus)> 0.01) 
      
       beta_s = 0.015
        s_inf = alpha_s / (alpha_s + beta_s)
        tau_s = 1 / (alpha_s + beta_s)
        
        ds_dt = (s_inf - Potassium_s) / tau_s
        s =  * ds_dt + Potassium_s
      

 2.5 Inward high-threshold Ca current I_CaH
    I_CaH  =  g_CaH * r * r * (V_dend - V_Ca)
    
        alpha_r = 1.7 / (1 + exp( -(V_dend - 5) / 13.9))
        beta_r = 0.02 * (V_dend + 8.5) / (exp((V_dend + 8.5) / 5) - 1)

            r_inf = alpha_r / (alpha_r + beta_r)
            tau_r = 5 / (alpha_r + beta_r)
        
                dr_dt = (r_inf - Calcium_r) / tau_r 
         
                    r = dr_dt + Calcium_r
    
 2.6 !!!! Calcium concentration !!!!
    dCa_dt = -3 * I_CaH - 0.075 * Ca2Plus * arbitrary
    Ca2Plus =  dCa_dt + Ca2Plus  uM 
    
 2.7 Leakage current I_ld
    I_ld   =  g_Ld * (V_dend - V_l)

 2.8 ***** GABA A current *****
    I_gab_dend   = gbar_gaba_dend * g_gaba_dend * (V_dend - V_gaba_dend)



 [================================================]
   3 axon hillock components
 [================================================]


                               3.1 Inward Sodium

    I_Na_a  = g_Na_a  * m_a * m_a * m_a * h_a * (V_axon - V_Na)

        m_inf_a   = 1 / (1 + (exp((-30 - V_axon)/ 5.5)))
        h_inf_a   = 1 / (1 + (exp((-60 - V_axon)/-5.8)))
        tau_h_a   =     1.5 * exp((-40 - V_axon)/33)

        dh_dt_a   = (h_inf_a - Sodium_h_a)/tau_h_a
         
        m_a       = m_inf_a
        h_a       = Sodium_h_a + dh_dt_a

 3.4 Potassium (transient)
    
    I_K_a   = g_K_a * (x_a ^ 4) * (V_axon - V_K)
    
        alpha_x_a = 1e-3*0.13 * (V_axon + 25) / (1 - exp(-(V_axon + 25) / 10))
        beta_x_a  = 1e-3*1.69 * exp(-0.0125 * (V_axon + 35))

        x_inf_a   = alpha_x_a / (alpha_x_a + beta_x_a)
        tau_x_a   =         1 / (alpha_x_a + beta_x_a)

        dx_dt_a   = (x_inf_a - Potassium_x_a) / tau_x_a 

        x_a       =dx_dt_a + Potassium_x_a


 3.2 Leak
    I_la    = g_la    * (V_axon - V_l)

  3.3 Soma-axon interaction current I_sa
    I_sa    = (g_int / p2) * (V_axon - V_soma)

    
  

  
  