function out = ode_update_cells(state, DYdt)


dh_dt   =       DYdt( 1 ,:) ; % Sodium_h
dn_dt   =       DYdt( 2 ,:) ; % Potassium_n
dx_dt_s =       DYdt( 3 ,:) ; % Potassium_x_s
dk_dt   =       DYdt( 4 ,:) ; % Calcium_k 
dl_dt   =       DYdt( 5 ,:) ; % Calcium_l
m       =       DYdt( 6 ,:)  
dr_dt   =       DYdt( 7 ,:) ; % Calcium_r 
ds_dt   =       DYdt( 8 ,:) ; % Potassium_s
dq_dt   =       DYdt( 9 ,:) ; % Hcurrent_q
dCa_dt  =       DYdt( 10,:) ; % Ca2Plus
m_a     =       Dydt( 11,:) ; % axon: m_a = m_a_inf 
dh_dt_a =       DYdt( 12,:) ; % Sodium_h_a
dx_dt_a =       DYdt( 13,:) ; % Potassium_x_a
DYdt(14:49,:) = 0

V_gaba_dend         = state.V_gaba_dend     ;
V_gaba_soma         = state.V_gaba_soma     ;
V_ampa_soma         = state.V_ampa_soma     ;
gbar_gaba_soma      = state.gbar_gaba_soma  ;
g_gaba_soma         = state.g_gaba_soma     ;
g_ampa_soma         = state.g_ampa_soma     ;
gbar_ampa_soma      = state.gbar_ampa_soma  ;
gbar_gaba_dend      = state.gbar_gaba_dend  ;
g_gaba_dend         = state.g_gaba_dend     ;
I_cx36 = state.I_cx36(:);
curr_noise  state.curr_noise(:);
vclamp state.vclamp(:)



out.V_soma          =  state.V_soma        + DYdt( 1 , : );; % V_soma               = DYdt( 1 );
out.Sodium_h        =  state.Sodium_h      + DYdt( 2 , : );; % Sodium_h             = DYdt( 2 );
out.Potassium_n     =  state.Potassium_n   + DYdt( 3 , : );; % Potassium_n          = DYdt( 3 );
out.Potassium_x_s   =  state.Potassium_x_s + DYdt( 4 , : );; % Potassium_x_s        = DYdt( 4 );
out.Calcium_k       =  state.Calcium_k     + DYdt( 5 , : );; % Calcium_k            = DYdt( 5 );
out.Calcium_l       =  state.Calcium_l     + DYdt( 6 , : );; % Calcium_l            = DYdt( 6 );
out.V_dend          =  state.V_dend        + DYdt( 7 , : );; % V_dend               = DYdt( 7 );
out.Calcium_r       =  state.Calcium_r     + DYdt( 8 , : );; % Calcium_r            = DYdt( 8 );
out.Potassium_s     =  state.Potassium_s   + DYdt( 9 , : );; % Potassium_s          = DYdt( 9 );
out.Hcurrent_q      =  state.Hcurrent_q    + DYdt( 10, : );; % Hcurrent_q           = DYdt( 10);
out.Ca2Plus         =  state.Ca2Plus       + DYdt( 21, : );; % Ca2Plus              = DYdt( 21);
out.V_axon          =  state.V_axon        + DYdt( 22, : );; % I_CaH                = DYdt( 22);
out.Sodium_m_a      =  state.Sodium_m_a    + DYdt( 23, : );; % V_axon               = DYdt( 23);
out.Sodium_h_a      =  state.Sodium_h_a    + DYdt( 24, : );; % Sodium_h_a           = DYdt( 24);
out.Potassium_x_a   =  state.Potassium_x_a + DYdt( 25, : );; % Potassium_x_a        = DYdt( 25);
out.I_cx36          = DYdt( 26);
curr_noise          = DYdt( 27);
% vclamp            = DYdt( 28);


g_CaL           =     cell_parameters.g_CaL         ;
g_int           =     cell_parameters.g_int         ;
g_K_Ca          =     cell_parameters.g_K_Ca        ;
g_ld            =     cell_parameters.g_ld          ;
C_m             =     cell_parameters.C_m           ;
g_Na_s          =     cell_parameters.g_Na_s        ;
g_Kdr_s         =     cell_parameters.g_Kdr_s       ;
g_K_s           =     cell_parameters.g_K_s         ;
g_ls            =     cell_parameters.g_ls          ;
g_CaH           =     cell_parameters.g_CaH         ;
g_h             =     cell_parameters.g_h           ;
g_Na_a          =     cell_parameters.g_Na_a        ;
g_K_a           =     cell_parameters.g_K_a         ;
g_la            =     cell_parameters.g_la          ;
p1              =     cell_parameters.p1            ;
p2              =     cell_parameters.p2            ;
V_Na            =     cell_parameters.V_Na          ;
V_K             =     cell_parameters.V_K           ;
V_Ca            =     cell_parameters.V_Ca          ;
V_h             =     cell_parameters.V_h           ;
V_l             =     cell_parameters.V_l           ;
V_gaba_dend     =     cell_parameters.V_gaba_dend   ;
V_gaba_soma     =     cell_parameters.V_gaba_soma   ;
V_ampa_soma     =     cell_parameters.V_ampa_soma   ;
gbar_gaba_soma   =    cell_parameters.gbar_gaba_soma;
gbar_ampa_soma   =    cell_parameters.gbar_ampa_soma;
gbar_gaba_dend   =    cell_parameters.gbar_gaba_dend;





    % SOMATIC CURRENTS
        
    % Dendrite-soma interaction current
    I_ds  = (g_int ./ p1) .* (V_soma - V_dend) 
    % Inward low-threshold Ca current
    I_CaL = g_CaL .* k .* k .* k .* l .* (V_soma - V_Ca)
    % Inward Na current
    I_Na_s  = g_Na_s .* m .* m .* m .* h .* (V_soma - V_Na)
    % Leak current
    I_ls  = g_ls .* (V_soma - V_l)
    % Potassium current
    I_Kdr_s = g_Kdr_s .* n .* n .* n .* n .* (V_soma - V_K)
    I_K_s   = g_K_s .* (x_s .^ 4) .* (V_soma - V_K)
    % Axon-soma interaction current
    I_as    = (g_int ./ (1 - p2)) .* (V_soma - V_axon)
    %.*.*.*.*.* AMPA current .*.*.*.*.*
    I_amp   = gbar_ampa_soma .* g_ampa_soma .* (V_dend - V_ampa_soma)
    %.*.*.*.*.* GABA A current .*.*.*.*.*
    I_gab_soma   = gbar_gaba_soma .* g_gaba_soma .* (V_soma - V_gaba_soma)
        

    
    % DENDRITIC CURRENTS
    
    % Soma-dendrite interaction current I_sd
    I_sd   = (g_int ./ (1 - p1)) .* (V_dend - V_soma)
    % Inward high-threshold Ca current I_CaH
    I_CaH  =  g_CaH .* r .* r .* (V_dend - V_Ca)
    % Outward Ca-dependent K current I_K_Ca
    I_K_Ca =  g_K_Ca .* s .* (V_dend - V_K)
    % Leakage current I_ld
    I_ld   =  g_ld .* (V_dend - V_l)
    % Inward anomalous rectifier I_h
    I_h    =  g_h .* q .* (V_dend - V_h)
    %.*.*.*.*.* GABA A current .*.*.*.*.*
    I_gab_dend   = gbar_gaba_dend .* g_gaba_dend .* (V_dend - V_gaba_dend)
    

    
    % AXONAL CURRENTS

    % Sodium
    I_Na_a  = g_Na_a  .* m_a .* m_a .* m_a .* h_a .* (V_axon - V_Na)
    % Leak
    I_la    = g_la    .* (V_axon - V_l)
    % Soma-axon interaction current I_sa
    I_sa    = (g_int ./ p2) .* (V_axon - V_soma)
    % Potassium (transient)
    I_K_a   = g_K_a .* (x_a .^ 4) .* (V_axon - V_K)
    
    
    %% update voltages

    dVs_dt = (-(I_CaL   + I_ds  + I_as + I_Na_s + I_ls   + I_Kdr_s + I_K_s  + I_amp + I_gab_soma) ) ./ C_m
    dVd_dt = (-(I_CaH   + I_sd  + I_ld + I_K_Ca + I_cx36 + I_h     + I_gab_dend) + curr_noise ) ./ C_m
    dVa_dt = (-(I_K_a   + I_sa  + I_la + I_Na_a)                                   ) ./ C_m


    V_soma =  V_soma + dVs_dt 
    V_dend =  V_dend + dVd_dt 
    V_axon =  V_axon + dVa_dt 


out.I_CaL           = I_CaL   ;
out.I_ds            = I_ds    ;
out.I_as            = I_as    ;
out.I_Na_s          = I_Na_s  ;
out.I_ls            = I_ls    ;
out.I_Kdr_s         = I_Kdr_s ;
out.I_K_s           = I_K_s   ;
out.I_CaH           = I_CaH   ;
out.I_sd            = I_sd    ;
out.I_ld            = I_ld    ;
out.I_K_Ca          = I_K_Ca  ;
out.I_cx36          = I_cx36  ;
out.I_h             = I_h     ;
out.I_K_a           = I_K_a   ;
out.I_sa            = I_sa    ;
out.I_la            = I_la    ;
out.I_Na_a          = I_Na_a  ;

out.V_soma          =  V_soma        ;
out.Sodium_h        =  Sodium_h      ;
out.Potassium_n     =  Potassium_n   ;
out.Potassium_x_s   =  Potassium_x_s ;
out.Calcium_k       =  Calcium_k     ;
out.Calcium_l       =  Calcium_l     ;
out.V_dend          =  V_dend        ;
out.Calcium_r       =  Calcium_r     ;
out.Potassium_s     =  Potassium_s   ;
out.Hcurrent_q      =  Hcurrent_q    ;
out.Ca2Plus         =  Ca2Plus       ;
out.V_axon          =  V_axon        ;
out.Sodium_m_a      =  Sodium_m_a    ;
out.Sodium_h_a      =  Sodium_h_a    ;
out.Potassium_x_a   =  Potassium_x_a ;
                    

