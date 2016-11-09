function out = ode_update_cells(state, cell_parameters, DYdt)

h    =      DYdt(1  ,:);
n    =      DYdt(2  ,:);
x_s  =      DYdt(3  ,:);
k    =      DYdt(4  ,:);
l    =      DYdt(5  ,:);
m    =      DYdt(6  ,:);
r    =      DYdt(7  ,:);
s    =      DYdt(8  ,:);
q    =      DYdt(9  ,:);
Ca   =      DYdt(10 ,:);
m_a  =      DYdt(11 ,:);
h_a     =   DYdt(12 ,:);
x_a     =   DYdt(13 ,:);

% x_dt_s     =   Dydt(11 ,1);

V_soma          =  state.V_soma;
V_dend          =  state.V_dend;
V_axon          =  state.V_axon;

g_gaba_soma         = state.g_gaba_soma;
g_ampa_soma         = state.g_ampa_soma;
g_gaba_dend         = state.g_gaba_dend;
I_cx36              = state.I_cx36(:);
curr_noise          = state.curr_noise(:);

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


    keyboard
    % SOMATIC CURRENTS
        
    % Dendrite-soma interaction current
    I_ds  = (g_int ./ p1) .* (V_soma - V_dend);
    % Inward low-threshold Ca current
    I_CaL = g_CaL .* k .* k .* k .* l .* (V_soma - V_Ca);
    % Inward Na current
    I_Na_s  = g_Na_s .* m .* m .* m .* h .* (V_soma - V_Na);
    % Leak current
    I_ls  = g_ls .* (V_soma - V_l);
    % Potassium current
    I_Kdr_s = g_Kdr_s .* n .* n .* n .* n .* (V_soma - V_K);
    I_K_s   = g_K_s .* (x_s .^ 4) .* (V_soma - V_K);
    % Axon-soma interaction current
    I_as    = (g_int ./ (1 - p2)) .* (V_soma - V_axon);
    %.*.*.*.*.* AMPA current .*.*.*.*.*
    I_amp   = gbar_ampa_soma .* g_ampa_soma .* (V_dend - V_ampa_soma);
    %.*.*.*.*.* GABA A current .*.*.*.*.*
    I_gab_soma   = gbar_gaba_soma .* g_gaba_soma .* (V_soma - V_gaba_soma);
        

    
    % DENDRITIC CURRENTS
    
    % Soma-dendrite interaction current I_sd
    I_sd   = (g_int ./ (1 - p1)) .* (V_dend - V_soma);
    % Inward high-threshold Ca current I_CaH
    I_CaH  =  g_CaH .* r .* r .* (V_dend - V_Ca);
    % Outward Ca-dependent K current I_K_Ca
    I_K_Ca =  g_K_Ca .* s .* (V_dend - V_K);
    % Leakage current I_ld
    I_ld   =  g_ld .* (V_dend - V_l);
    % Inward anomalous rectifier I_h
    I_h    =  g_h .* q .* (V_dend - V_h);
    %.*.*.*.*.* GABA A current .*.*.*.*.*
    I_gab_dend   = gbar_gaba_dend .* g_gaba_dend .* (V_dend - V_gaba_dend);
    

    
    % AXONAL CURRENTS

    % Sodium
    I_Na_a  = g_Na_a  .* m_a .* m_a .* m_a .* h_a .* (V_axon - V_Na);
    % Leak
    I_la    = g_la    .* (V_axon - V_l);
    % Soma-axon interaction current I_sa
    I_sa    = (g_int ./ p2) .* (V_axon - V_soma);
    % Potassium (transient)
    I_K_a   = g_K_a .* (x_a .^ 4) .* (V_axon - V_K);
    
    
    %% update voltages

    dVs_dt = (-(I_CaL   + I_ds  + I_as + I_Na_s + I_ls   + I_Kdr_s + I_K_s  + I_amp + I_gab_soma) ) ./ C_m;
    dVd_dt = (-(I_CaH   + I_sd  + I_ld + I_K_Ca + I_cx36 + I_h     + I_gab_dend) + curr_noise ) ./ C_m;
    dVa_dt = (-(I_K_a   + I_sa  + I_la + I_Na_a)                                   ) ./ C_m;


    out.V_soma =  V_soma + dVs_dt ;
    out.V_dend =  V_dend + dVd_dt ;
    out.V_axon =  V_axon + dVa_dt ;


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
out.Sodium_h        =  h      ;
out.Potassium_n     =  n   ;
out.Potassium_x_s   =  x_s ;
out.Calcium_k       =  k     ;
out.Calcium_l       =  l     ;
out.V_dend          =  V_dend        ;
out.Calcium_r       =  r     ;
out.Potassium_s     =  s   ;
out.Hcurrent_q      =  q    ;
out.Ca2Plus         =  Ca       ;
out.V_axon          =  V_axon        ;
out.Sodium_m_a      =  m_a    ;
out.Sodium_h_a      =  h_a    ;
out.Potassium_x_a   =  x_a ;
                    

