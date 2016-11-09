% calc_steady_state.m
%
% Sanity checks for continuating simulations from saved last state.

gpu = 1;
netsize = [1 2 1];
simtime  = 800;
W = [0 1; 1 0];
dt = .05

cell_parameters = createDefaultNeurons(2);
cell_parameters.g_CaL = [0.4 1.1]; % overwrite defaults

toreport =  {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r', ...
 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'I_CaH', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a' , 'Ca2_soma'};

simtime = 200;
rndState = initNetState(prod(netsize),1, 0);
rndState.V_soma = [-100 -80]';
st_st = IOnet_new('delta', dt, 'networksize', netsize, 'time',simtime ,'W', W*0 ,'to_report',toreport,'gpu', gpu, 'tempState', rndState,'cell_parameters', cell_parameters);

% to continue the network and see if it stitches well
    simtime = 200;
    cont  = IOnet_new('tempState', st_st.lastState ,'cell_parameters', cell_parameters, 'networksize', netsize,'time',simtime ,'W', W*0  ,'to_report',toreport,'gpu', gpu);
    plot([st_st.networkHistory.V_soma cont.networkHistory.V_soma]')
