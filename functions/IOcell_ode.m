function [DYdt] =  IOcell_ode(t, Y)

V_soma               = Y( 1 );
Sodium_h             = Y( 2 );
Potassium_n          = Y( 3 );
Potassium_x_s        = Y( 4 );
Calcium_k            = Y( 5 );
Calcium_l            = Y( 6 );
V_dend               = Y( 7 );
Calcium_r            = Y( 8 );
Potassium_s          = Y( 9 );
Hcurrent_q           = Y( 10);
Ca2Plus              = Y( 11);
I_CaH                = Y( 12);
V_axon               = Y( 13);
Sodium_h_a           = Y( 14);
Potassium_x_a        = Y( 15);
I_cx36               = Y( 16);
curr_noise           = Y( 17);
vclamp               = Y( 18);
g_CaL                = Y( 19);
g_int                = Y( 20);
g_K_Ca               = Y( 21);
g_ld                 = Y( 22);
C_m                  = Y( 23);
g_Na_s               = Y( 24);
g_Kdr_s              = Y( 25);
g_K_s                = Y( 26);
g_ls                 = Y( 27);
g_CaH                = Y( 28);
g_h                  = Y( 29);
g_Na_a               = Y( 30);
g_K_a                = Y( 31);
g_la                 = Y( 32);
p1                   = Y( 33);
p2                   = Y( 34);
V_Na                 = Y( 35);
V_K                  = Y( 36);
V_Ca                 = Y( 37);
V_h                  = Y( 38);
V_l                  = Y( 39);



% for defaults see createDefaultNeurons.m
% [================================================]
%    Change log
% [================================================]
%  see #change hashtag
%  



%% New state...

    %%================================================]
    %          update somatic components
    %=================================================]
    
    k_inf = (1 / (1 + exp(-1 * (V_soma + 61)   / 4.2)));
    l_inf = (1 / (1 + exp((     V_soma + 85.5) / 8.5)));

    tau_k = 1;
    tau_l = ((20 * exp((V_soma + 160) / 30) / (1 + exp((V_soma + 84) / 7.3))) +35);
        
    dk_dt = (k_inf - Calcium_k) / tau_k;
    dl_dt = (l_inf - Calcium_l) / tau_l;

        
    m_inf = 1 / (1 + (exp((-30 - V_soma)/ 5.5)));
    h_inf = 1 / (1 + (exp((-70 - V_soma)/-5.8)));
    tau_h =       3 * exp((-40 - V_soma)/33);

    dh_dt = (h_inf - Sodium_h)/tau_h;
    m     = m_inf;
    

     n_inf = 1 / (1 + exp( ( -3 - V_soma) /  10));
     tau_n =   5 + (  47 * exp( -(-50 - V_soma) /  900));
     dn_dt = (n_inf - Potassium_n) / tau_n;
          
     alpha_x_s = 0.13 * (V_soma + 25) / (1 - exp(-(V_soma + 25) / 10));
     beta_x_s  = 1.69 * exp(-0.0125 * (V_soma + 35));
    
     x_inf_s   = alpha_x_s / (alpha_x_s + beta_x_s);
     tau_x_s   =         1 / (alpha_x_s + beta_x_s);
    
    
     dx_dt_s   = (x_inf_s - Potassium_x_s) / tau_x_s;
          
   
    %%================================================]
    % update dendritic components
    %=================================================]
    
    % Update dendritic H current component
    q_inf = 1 /(1 + exp((V_dend + 80) / 4));
    tau_q = .1 /(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87));

    dq_dt = (q_inf - Hcurrent_q) / tau_q;

    % Update dendritic high-threshold Ca current component
     alpha_r = 1.7 / (1 + exp( -(V_dend - 5) / 13.9));
      beta_r = 0.02 * (V_dend + 8.5) / (exp((V_dend + 8.5) / 5) - 1);

       r_inf = alpha_r / (alpha_r + beta_r);
       tau_r = 5 / (alpha_r + beta_r);
        
        dr_dt = (r_inf - Calcium_r) / tau_r; 
         


      % Update dendritic Ca-dependent K current component

  % original:
       alpha_s = (0.00002 * Ca2Plus) * (0.00002 * Ca2Plus < 0.01) + 0.01*((0.00002 * Ca2Plus)> 0.01); % nefarious calcium creates instabilities?
       beta_s = 0.015;
    
    % NEW:
        % alpha_s =  .9e-3*log(Ca2Plus)*exp(V_dend/24);
        % alpha_s =  .9e-4*log(Ca2Plus)*exp(V_dend/24);
        % beta_s  =  .75e-4*exp(-V_dend/24);

        s_inf = alpha_s / (alpha_s + beta_s);
        tau_s = .1 / (alpha_s + beta_s);
        
        ds_dt = (s_inf - Potassium_s) / tau_s;
        
        
    % update Calcium concentration
    dCa_dt = -3 * I_CaH - 0.075 * Ca2Plus;
    
       % [================================================]
       %%   update axon hillock components
       % [================================================]
       
       % Update axonal Na components
       % NOTE: current has shortened inactivation to account for high
       % firing frequencies in axon hillock
     m_inf_a   = 1 / (1 + (exp((-30 - V_axon)/ 5.5)));
     h_inf_a   = 1 / (1 + (exp((-60 - V_axon)/-5.8)));
     tau_h_a   =     1.5 * exp((-40 - V_axon)/33);

     dh_dt_a   = (h_inf_a - Sodium_h_a)/tau_h_a;
     
     m_a       = m_inf_a;

     % Update potassium components
     alpha_x_a = 0.13 * (V_axon + 25) / (1 - exp(-(V_axon + 25) / 10));
     beta_x_a  = 1.69 * exp(-0.0125 * (V_axon + 35));
    
     x_inf_a   = alpha_x_a / (alpha_x_a + beta_x_a);
     tau_x_a   =         1 / (alpha_x_a + beta_x_a);
    
     dx_dt_a   = (x_inf_a - Potassium_x_a) / tau_x_a; 
    

     % axonal spike creates huge potassium 


DYdt(1     ,1) = dh_dt;                          % Sodium_h;
DYdt(2     ,1) = dn_dt;                        % Potassium_n;
DYdt(3     ,1) = dx_dt_s;                      % Potassium_x_s;
DYdt(4     ,1) = dk_dt;                         % Calcium_k ;
DYdt(5     ,1) = dl_dt;                          % Calcium_l;
DYdt(6     ,1) = m;
DYdt(7     ,1) = dr_dt; ;                         % Calcium_r ;
DYdt(8     ,1) = ds_dt;;                        % Potassium_s;
DYdt(9     ,1) = dq_dt;;                         % Hcurrent_q;
DYdt(10    ,1) = dCa_dt;                            % Ca2Plus;
DYdt(11    ,1) = m_a;
DYdt(12    ,1) = dh_dt_a;                         % Sodium_h_a;
DYdt(13    ,1) = dx_dt_a;                      % Potassium_x_a;
DYdt(14:39 ,1) = 0;



