% asymetry tests.m

% This script shows how coupling coefficient varies as a function of CaT expression

% Is =gsð(Vm - VRevÞ)
% Vs = gs Rm Vm + gs Rm V Rev:

% From Lefler et al. 2014
% 
% Thus, plotting Vs as a function of Vm will yield a linear relation where the slope is equal to gs * Rm.
%  Normalizing the slope by the input resistance will thus give the synaptic conductance.

% The coupling coefficient (CC) was calculated as the ratio between steady- state voltage change 
% of the postjunctional cell (DVpost) and that of the prejunctional cell (DVpre) in response 
% to current injection in the prejunctional cell. 
% 
% The DV was calculated as the average voltage drop during the last 35–50 ms of the injected pulse.
 % CCs smaller than 0.002 were treated as zero.
 % 
% The values for CC, their symmetry, and the effect of light activation were measured nearly 
% simultaneously, using a protocol that obtained these values within a single sweep. A light pulse 
% train was followed by two hyperpolarizing current injections to one cell with the second pulse given
% simultaneously with a light pulse train. This was followed by same protocol given to the other cell.
% This protocol was repeated 20–30 times and the recorded signals were averaged. The postjunctional 
% voltage response during light stimulation was obtained by subtracting the light-evoked response 
% from the combined response to light train and current injection in the prejunctional cell.


% CC = DVpost/DVpre;


clear

gpu = 0;


% some defaults
to_report =  {'V_soma', 'I_cx36', 'V_dend'};
pert = []; I_app = [];
dt = 0.025;
cell_function = 'vanilla';
% cell_function = 'devel'

% default noise
stnoise = [0 0 0 0];
ounoise = [3 5 0 1];

Ca_set = [.3:.2:.5];


% available variables for reporting:
% to_report =  {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r', ...
%  'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a' , 'g_GABA_soma','g_GABA_dend', 'g_AMPA', 'Ca2_soma'};


experiment = 1;
switch experiment
	
	case 1 % stimulate first cell then second cell to compute asymmetries

		% to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend'};
	
	 % to_report = {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', ...
	 % 'Potassium_s', 'Hcurrent_q', 'Ca2Plus',  'Sodium_m_a', 'curr_noise'};
		netsize = [1 2 1];
		simtime  = 2000;
		gap = .1;
		W = [0 1; 1 0]*gap;
		dt = 0.01;
		I_app = zeros(2, simtime*(1/dt));
		I_app(1,300*(1/dt):1000*(1/dt))   = -5;
		
		cell_parameters = createDefaultNeurons(2);

		ounoise = [0 0 0 0];



		
	
	case 2 % compute coupling coefficient with different calcium T types
	

		% to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend'};	
	 
		netsize = [1 2 1];
		simtime  = 2000;
		gap = .5;
		W = [0 1; 1 0]*gap;
		I_app = zeros(2, simtime*(1/dt));
		I_app(2,300*(1/dt):800*(1/dt))   = -10;
		
		cell_parameters = createDefaultNeurons(2);
		cell_parameters.g_CaL = [.4 1.1]';


	% case 3 % the cell surrounded 3D
	% 	% to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend','curr_noise'};
	% 	simtime = 1000;
	% 	netsize = [3 3 3];
	% 	gap = .05;
	% 	centercell = 14;
	% 	W = zeros(prod(netsize)); 
	% 	allcells = 1:prod(netsize);
	% 	othercells = setdiff(allcells, centercell);
	% 	W(centercell, othercells) = 1*gap;
	% 	W = W+W';

	% 	cell_parameters = createDefaultNeurons(27);
	% 		cell_parameters.g_CaL = ones(1,prod(netsize))*1.;
	% 		cell_parameters.g_CaL(centercell) = 1.1;
	% 	I_app = zeros(prod(netsize), simtime*(1/dt));
	% 	I_app(centercell,300*(1/dt):400*(1/dt))   = 10;


	% case 4 % the cell surrounded 2D #jochen #xcorrs
	% 	to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend'};
	% 	simtime = 1000;
	% 	netsize = [1 3 3];
	% 	gap = .05;
	% 	centercell = 5;
	% 	W = zeros(prod(netsize)); 
	% 	allcells = 1:prod(netsize);
	% 	othercells = setdiff(allcells, centercell);
	% 	W(centercell, othercells) = 1*gap;
	% 	W = W+W';

	% 	cell_parameters = createDefaultNeurons(9);
	% 		cell_parameters.g_CaL = ones(1,prod(netsize))*.3;
	% 		cell_parameters.g_CaL(centercell) = 1.;
			
	% 		cell_parameters.g_int = ones(1,prod(netsize))*.15;
	% 		cell_parameters.g_int(centercell) = .15;

	% 	I_app = zeros(prod(netsize), simtime*(1/dt));
	% 	I_app(centercell,300*(1/dt):400*(1/dt))   = 10;


	% case 5 % the cell surrounded 2D
	% 	% to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend','curr_noise'};
	% 	simtime = 1000;
	% 	netsize = [1 3 3];
	% 	gap = .05;
	% 	centercell = 5;
	% 	W = zeros(prod(netsize));
	% 	allcells = 1:prod(netsize);
	% 	othercells = setdiff(allcells, centercell);
	% 	W(centercell, othercells) = 1*gap;
	% 	W = W+W';

	% 	cell_parameters = createDefaultNeurons(9);
	% 		cell_parameters.g_CaL = ones(1,prod(netsize))*.4;
	% 		cell_parameters.g_CaL(centercell) = 1.;
			
	% 		cell_parameters.g_int = ones(1,prod(netsize))*.05;
	% 		cell_parameters.g_int(centercell) = .15;

	% 	I_app = zeros(prod(netsize), simtime*(1/dt));
	% 	I_app(centercell,300*(1/dt):400*(1/dt))   = 10;



	% case 6 % a chain of cells
	% 	% to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend','curr_noise'};
	% 	simtime = 1000;
	% 	netsize = [1 6 1];
	% 	gap = .1;
	% 	stimcell = 1;
	% 	W = [0 1 0 0 0 0;
	% 		 0 0 1 0 0 0;
	% 		 0 0 0 1 0 0;
	% 		 0 0 0 0 1 0;
	% 		 0 0 0 0 0 1;
	% 		 0 0 0 0 0 0];
		
	% 	% W = W+W';
	% 	W = W*gap;

	% 	cell_parameters = createDefaultNeurons(6);
	% 		% cell_parameters.g_CaL = ones(linspace(.4,1,6);
	% 	I_app = zeros(prod(netsize), simtime*(1/dt));
	% 	I_app(1,300*(1/dt):302*(1/dt))   = 10;
	% 	% I_app(6,800*(1/dt):820*(1/dt))   = 30;


	% case 7 % rebound spike and sync after depo
	% 	dt = 0.01;
	% 	simtime = 1000;
	% 	netsize = [1 3 3];
	% 	gap = .05;
	% 	W = zeros(9);
	% 	W(5,[1 2 3 4 6 7 8 9]) = 1*gap;
	% 	W = W+W';
	% 	cell_parameters = createDefaultNeurons(9);
	% 	cell_parameters.g_CaL = ones(1,9)*.9;
	% 	% cell_parameters.g_ih = ones(1,9)*1.2;
	% 	% g_CaL(5) = 1.;
	% 	I_app = zeros(9, simtime*(1/dt));
	% 	I_app(5,(100*(1/dt):105*(1/dt))) = 10;% uA/cm^2 
	% 	stnoise = [3 5 0 0]; % 
	% 	ounoise = [3 1 0 0]; % 
	% 	% cell_function = 'devel'; 
	% 	cell_function = 'vanilla'; 

	% case 8 % rebound spike and sync after depo - cx36 test
	% 	dt = 0.01;
	% 	simtime = 1000;
	% 	netsize = [1 3 3];
	% 	gap = .0;
	% 	W = zeros(9);
	% 	W(5,[1 2 3 4 6 7 8 9]) = 1*gap;
	% 	W = W+W';
	% 	cell_parameters = createDefaultNeurons(9);
	% 	cell_parameters.g_CaL = ones(1,9)*.9;
	% 	% g_CaL(5) = 1.;
	% 	I_app = zeros(9, simtime*(1/dt));
	% 	I_app(5,(100*(1/dt):120*(1/dt))) = 10;% uA/cm^2 
	% 	stnoise = [3 5 0 0]; % 
	% 	ounoise = [3 1 0 0]; % 
	% 	% cell_function = 'devel'; 
	% 	cell_function = 'vanilla'; 


end	


st_st = IOnet_new('delta', 0.05, 'networksize', netsize, 'time',2000 ,'W', W*0 ,'to_report',{'V_soma'},'gpu', gpu,'cell_parameters', cell_parameters , 'ou_noise', stnoise);

	% to continue the network and see if it stitches well


f = 0; g =0;
for ca1 = Ca_set
	 f = f +1;
	for ca2 = Ca_set
		g = g + 1
	
		cell_parameters.g_CaL = [ca1 ca2]';

		if f >= g

			[transients] = ...
				  IOnet_new('tempState', st_st.lastState ,'delta', dt, 'cell_parameters', cell_parameters,  'cell_function', cell_function, 'networksize', netsize,'perturbation', pert ,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', ounoise,'to_report', to_report,'gpu', gpu);


			% replayResults(transients, [],0, 1)


			Vr1 = mean(transients.networkHistory.V_soma(1,200:299));
			Vr2 = mean(transients.networkHistory.V_soma(2,200:299));

			Vst1_stim1 = mean(transients.networkHistory.V_soma(1,800:1000)) ;
			Vst2_stim1 = mean(transients.networkHistory.V_soma(2,800:1000)) ;
			
			DVpre_stim1  = Vr1- Vst1_stim1;
			DVpost_stim1 = Vr2- Vst2_stim1;

			CC1(f,g) = DVpost_stim1/DVpre_stim1;
		end

	end
	g = 0;

end
