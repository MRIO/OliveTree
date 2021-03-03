% leak_compensation.m
% for a given set of cells calculate necessary leak compensation
 % to balance excitability of cells without gap junction
 % 
 % gl(mean(V(t))-Vl) ~ sum(gc_i(mean(Vi(t) - Vj(t) ) )


% two strategies:
% 1. decrease leak in networks without gap
% 2. increase leak in networks with gaps

% without gaps -> decrease leak -> smaller depo current -> lower Vm -> less excitability
% with gaps -> increase leak -> smaller depo current -> lower Vm -> less excitability



given:
1. total gap junction conductance per cell
2. average membrane potential of cells in disconnected network

return:
1. vector of values of leaks to compensate for gap junction leak difference

% Somatic conductances (mS/cm2)
cell_parameters.g_CaL    =  1.1   .* O; 
cell_parameters.g_Na_s   =  150   .* O;      % Sodium
cell_parameters.g_Kdr_s  =  9.0   .* O;      % Potassium
cell_parameters.g_K_s    =  5     .* O;      % Potassium
cell_parameters.g_ls     =  0.016 .* O;      % Leaks
    
% Dendritic conductances (mS/cm2)
cell_parameters.g_K_Ca   =  35      .* O;       % Potassium: not voltage dependent 
cell_parameters.g_CaH    =  4.5     .* O;     % High-threshold calcium
cell_parameters.g_ld     =  0.016   .* O;   % Leak
cell_parameters.g_h      =  .12     .* O;    % H current .12
cell_parameters.g_h_s    =  .12     .* O;    % H current, somatic

cell_parameters.V_Na =  55 .* O;       % Sodium
cell_parameters.V_K  = -75 .* O;       % Potassium
cell_parameters.V_Ca = 120 .* O;       % Calcium
cell_parameters.V_h  = -43 .* O;       % H current
cell_parameters.V_l  =  10 .* O;       % Leak

cell_parameters.V_gaba_dend = -70 .*O; % from Devor and Yarom, 2002
cell_parameters.V_gaba_soma = -63 .*O; % from Devor and Yarom, 2002
cell_parameters.V_ampa_soma = 0 	  .*O; % from Cian McDonnel et al 2012
cell_parameters.V_ampa_dend = 0 	  .*O; % from Cian McDonnel et al 2012
