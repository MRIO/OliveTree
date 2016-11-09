% hyperpol_tests.m

 % Note that 1 􏰈A/cm2 corresponds to 0.1 nA for a total cell surface of 10,000 􏰈m2./

clear

gpu = 1;
netsize = [1 6 1];
noneurons = prod(netsize);
simtime  = 1000;
W = zeros(noneurons);
% subset1 =  {'V_soma','V_dend','V_axon','Calcium_r', 'I_CaH','Ca2Plus', 'Potassium_s', 'Hcurrent_q'};
subset1 =  {'V_soma'};
dt = .05;
% cell_function = 'vanilla';
cell_function = 'devel';

def_neurons = createDefaultNeurons(6);
% def_neurons.g_CaL = [.4 .5 .6 1.1 1.2 1.3];
def_neurons.g_CaL = [.4 .4 .4 .4 .4 .4];
def_neurons.g_h   = linspace(.12, 1.2,6)';
% def_neurons.g_la   = linspace(.01, 0.02,6)';
% def_neurons.g_ls   = linspace(.01, 0.02,6)';
% def_neurons.g_ld   = linspace(.01, 0.02,6)';
% def_neurons.g_int  = linspace(.01, 0.02,6)';
% def_neurons.arbitrary = [1e3 1e3 1e3 1e0 1e0 1e0 ];

% calc steady state
st_st = IOnet_new('W', W, 'cell_function', cell_function, 'networksize', netsize, 'cell_parameters', def_neurons, 'time', 1000 ,'gpu', gpu,'to_report', {'V_soma'},'delta',0.05);

hyperpolsteps = 5;%; linspace(1,20, 5);


simcount= 0;
for hstep = hyperpolsteps
		simcount = simcount+1;

		
		I_app = zeros(noneurons, simtime*(1/dt));
		I_app(:,(100*(1/dt):500*(1/dt))) = -hstep; % nAmpere 20/dt [nA/s.cm^2] 
		% I_app(1,(351:360)*20) =  21/20; % nAmpere 20/dt [nA/s.cm^2] 

		gnoise = [0 0 0 0];



	[transients] = IOnet_new('tempState', st_st.lastState ,'cell_function', cell_function, 'cell_parameters', def_neurons, 'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', subset1,'gpu', gpu ,'cell_parameters', def_neurons, 'delta',dt);
	
		V_soma_cell_1(simcount,:) = transients.networkHistory.V_soma(1,:);
		V_soma_cell_2(simcount,:) = transients.networkHistory.V_soma(2,:);
		V_soma_cell_3(simcount,:) = transients.networkHistory.V_soma(3,:);
		V_soma_cell_4(simcount,:) = transients.networkHistory.V_soma(4,:);
		V_soma_cell_5(simcount,:) = transients.networkHistory.V_soma(5,:);
		V_soma_cell_6(simcount,:) = transients.networkHistory.V_soma(6,:);

end

% availablefieldstosave =  {'V_soma', 'I_cx36', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r',...
%  'I_CaH', 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a', 'Ca2_soma'};
% subset1 =  {'V_soma','V_dend','V_axon', 'I_cx36','Calcium_r', 'I_CaH','Ca2Plus', 'Potassium_s', 'Hcurrent_q'};


figure
subplot(2,3,1), plot(V_soma_cell_1'); ylim([-200,10])
subplot(2,3,2), plot(V_soma_cell_2'); ylim([-200,10])
subplot(2,3,3), plot(V_soma_cell_3'); ylim([-200,10])
subplot(2,3,4), plot(V_soma_cell_4'); ylim([-200,10])
subplot(2,3,5), plot(V_soma_cell_5'); ylim([-200,10])
subplot(2,3,6), plot(V_soma_cell_6'); ylim([-200,10])





