function ...
    [ V_soma,...            
 I_K_s,...                  
 I_Kdr_s,...                
 I_CaL,...                  
 I_ds,...                   
 m,...                      
 h,...                      
 n,...                      
 x_s,...                    
  k,...                     
  l,...                     
  V_dend,...                
  r,...                     
  s,...                     
  q,...                     
  Ca2Plus,...               
  I_CaH,...                 
  I_K_Ca,...                
  I_c,...                   
  V_axon,...                
  m_a,...                   
  h_a,...                   
  x_a] = IOcellFunctionCa_gpu(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,S26,S27,S28)


  V_soma             = S1 ; %V_soma;    
  PotassiumSlow_soma = S2; %I_K_s;     
  PotassiumFast_soma = S3; %I_Kdr_s;
  Calcium_soma       = S4; %I_CaL;     
  Leak_soma_dend     = S5; %I_ds;      
  Sodium_m           = S6; %m;         
  Sodium_h           = S7; %h;         
  Potassium_n        = S8; %n;         
  Potassium_x_s      = S9; %x_s;       
  Calcium_k          = S10; %k;        
  Calcium_l          = S11; %l;        
  V_dend             = S12; %V_dend;   
  Calcium_r          = S13; %r;        
  Potassium_s        = S14; %s;        
  Hcurrent_q         = S15; %q;        
  Ca2Plus            = S16; %Ca2Plus;  
  I_CaH              = S17; %I_CaH;    
  I_K_Ca             = S18; %I_K_Ca;   
  I_cx36             = S19; %I_c;      
  V_axon             = S20; %V_axon;   
  Sodium_m_a         = S21; %m_a;      
  Sodium_h_a         = S22; %h_a;      
  Potassium_x_a      = S23; %x_a;      


delta = S24;
I_c   = S25;
I_app = S26;
V_app = S27;
g_CaL = S28;



% >> A = rand(5000);         
% >> B = gpuArray(rand(5000));
% >> 
% >> tic,arrayfun(@exp,B);toc
% Elapsed time is 2.299962 seconds.
% >> tic,arrayfun(@exp,B);toc
% Elapsed time is 0.002737 seconds.
% >> tic,arrayfun(@exp,A);toc
% Elapsed time is 40.080442 seconds.
% >> 
% >> 



%% Cell properties

% Capacitance
C_m    =   1;
    
% Somatic conductances (mS/cm2)
g_Na_s   =  150;      % Sodium
g_Kdr_s  =    9.0;    % Potassium
g_K_s    =    5;      % Potassium
g_ls     =    0.016;  % Leak
    
% Dendritic conductances (mS/cm2)
g_K_Ca   =  35;       % Potassium
g_CaH    =   4.5;     % High-threshold calcium
g_ld     =   0.016;   % Leak
g_h      =   0.12;    % H current

% Axon hillock conductances (mS/cm2)
g_Na_a   =  240;      % Sodium
g_K_a    =   20;      % Potassium
g_la     =    0.016;  % Leak
    
% Cell morphology
p1     = 0.25;        % Cell surface ratio soma/dendrite
p2     = 0.15;        % Cell surface ratio axon(hillock)/soma

g_int  = 0.13;        % Cell internal conductance

%% Reversal potentials
V_Na =  55;       % Sodium
V_K  = -75;       % Potassium
V_Ca = 120;       % Calcium
V_h  = -43;       % H current
V_l  =  10;       % Leak

%% New state...

    %% update somatic components
    
    k_inf = (1 / (1 + exp(-1 * (V_soma + 61)   / 4.2)));
    l_inf = (1 / (1 + exp((     V_soma + 85.5) / 8.5)));

    tau_k = 1;
    tau_l = ((20 * exp((V_soma + 160) / 30) / (1 + exp((V_soma + 84) / 7.3))) +35);
        
    dk_dt = (k_inf - Calcium_k) / tau_k;
    dl_dt = (l_inf - Calcium_l) / tau_l;

        k = delta * dk_dt + Calcium_k;
        l = delta * dl_dt + Calcium_l;

    m_inf = 1 / (1 + (exp((-30 - V_soma)/ 5.5)));
    h_inf = 1 / (1 + (exp((-70 - V_soma)/-5.8)));
    tau_h =       3 * exp((-40 - V_soma)/33);

    dh_dt = (h_inf - Sodium_h)/tau_h;
     
    m     = m_inf;
    h     = Sodium_h + delta * dh_dt;

     n_inf = 1 / (1 + exp( ( -3 - V_soma) /  10));
     tau_n =   5 + (  47 * exp( -(-50 - V_soma) /  900));
     dn_dt = (n_inf - Potassium_n) / tau_n;
     n     = delta * dn_dt + Potassium_n;
          
     alpha_x_s = 0.13 * (V_soma + 25) / (1 - exp(-(V_soma + 25) / 10));
     beta_x_s  = 1.69 * exp(-0.0125 * (V_soma + 35));
    
     x_inf_s   = alpha_x_s / (alpha_x_s + beta_x_s);
     tau_x_s   =         1 / (alpha_x_s + beta_x_s);
    
     dx_dt_s   = (x_inf_s - Potassium_x_s) / tau_x_s;
    
     x_s       = 0.05 * dx_dt_s + Potassium_x_s;
          
          
    
    %% update dendritic components
    
    % Update dendritic H current component
    q_inf = 1 /(1 + exp((V_dend + 80) / 4));
    tau_q = 1 /(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87));
        
    dq_dt = (q_inf - Hcurrent_q) / tau_q;
        q = delta * dq_dt + Hcurrent_q;
        
     % Update dendritic high-threshold Ca current component
     alpha_r = 1.7 / (1 + exp( -(V_dend - 5) / 13.9));
      beta_r = 0.02 * (V_dend + 8.5) / (exp((V_dend + 8.5) / 5) - 1);
       r_inf = alpha_r / (alpha_r + beta_r);
       tau_r = 5 / (alpha_r + beta_r);
        
       dr_dt = (r_inf - Calcium_r) / tau_r;
           r = delta * dr_dt + Calcium_r;
          
      % Update dendritic Ca-dependent K current component
      % alpha_s = min([0.00002 * Ca2Plus, 0.01]);
      
      alpha_s = 0.00002 * Ca2Plus;
       beta_s = 0.015;
        s_inf = alpha_s / (alpha_s + beta_s);
        tau_s = 1 / (alpha_s + beta_s);
        
        ds_dt = (s_inf - Potassium_s) / tau_s;
        
            s = delta * ds_dt + Potassium_s;
        
       % update Calcium concentration
         dCa_dt = -3 * I_CaH - 0.075 * Ca2Plus;
        Ca2Plus = delta * dCa_dt + Ca2Plus;
       
       %% update axon hillock components
       
       % Update axonal Na components
       % NOTE: current has shortened inactivation to account for high
       % firing frequencies in axon hillock
     m_inf_a   = 1 / (1 + (exp((-30 - V_axon)/ 5.5)));
     h_inf_a   = 1 / (1 + (exp((-60 - V_axon)/-5.8)));
     tau_h_a   =     1.5 * exp((-40 - V_axon)/33);

     dh_dt_a   = (h_inf_a - Sodium_h_a)/tau_h_a;
     
     m_a       = m_inf_a;
     h_a       = Sodium_h_a + delta * dh_dt_a;
     
     % Update potassium components
     alpha_x_a = 0.13 * (V_axon + 25) / (1 - exp(-(V_axon + 25) / 10));
     beta_x_a  = 1.69 * exp(-0.0125 * (V_axon + 35));
    
     x_inf_a   = alpha_x_a / (alpha_x_a + beta_x_a);
     tau_x_a   =         1 / (alpha_x_a + beta_x_a);
    
     dx_dt_a   = (x_inf_a - Potassium_x_a) / tau_x_a;
    
     x_a       = 0.05 * dx_dt_a + Potassium_x_a;
     
     
        %% update currents
    
    % SOMATIC CURRENTS
        
    % Dendrite-soma interaction current
    I_ds  = (g_int / p1) * (V_soma - V_dend);
    % Inward low-threshold Ca current
    I_CaL = g_CaL * k * k * k * l * (V_soma - V_Ca);
    % Inward Na current
    I_Na_s  = g_Na_s * m * m * m * h * (V_soma - V_Na);
    % Leak current
    I_ls  = g_ls * (V_soma - V_l);
    % Potassium current
    I_Kdr_s = g_Kdr_s * n * n * n * n * (V_soma - V_K);
    I_K_s   = g_K_s * (x_s ^ 4) * (V_soma - V_K);
    % Axon-soma interaction current
    I_as    = (g_int / (1 - p2)) * (V_soma - V_axon);
    
    
    % DENDRITIC CURRENTS
    
    % Soma-dendrite interaction current I_sd
    I_sd   = (g_int / (1 - p1)) * (V_dend - V_soma);
    % Inward high-threshold Ca current I_CaH
    I_CaH  =  g_CaH * r * r * (V_dend - V_Ca);
    % Outward Ca-dependent K current I_K_Ca
    I_K_Ca =  g_K_Ca * s * (V_dend - V_K);
    % Leakage current I_ld
    I_ld   =  g_ld * (V_dend - V_l);
    % Inward anomalous rectifier I_h
    I_h    =  g_h * q * (V_dend - V_h);
    
    
    % AXONAL CURRENTS
    % Sodium
    I_Na_a  = g_Na_a  * m_a * m_a * m_a * h_a * (V_axon - V_Na);
    % Leak
    I_la    = g_la    * (V_axon - V_l);
    % Soma-axon interaction current I_sa
    I_sa    = (g_int / p2) * (V_axon - V_soma);
    % Potassium (transient)
    I_K_a   = g_K_a * (x_a ^ 4) * (V_axon - V_K);
    
    %% update voltages
    
    dVs_dt = (-(I_CaL   + I_ds  + I_as + I_Na_s + I_ls   + I_Kdr_s + I_K_s ) ) / C_m;
    dVd_dt = (-(I_CaH   + I_sd  + I_ld + I_K_Ca + I_c    + I_h)     + I_app) / C_m;
    dVa_dt = (-(I_K_a   + I_sa  + I_la + I_Na_a)                           ) / C_m;
        
    if isnan(V_app)
        V_soma = delta * dVs_dt + V_soma;
    else
        V_soma = V_app;
    end
    V_dend = delta * dVd_dt + V_dend;
    V_axon = delta * dVa_dt + V_axon;
    
