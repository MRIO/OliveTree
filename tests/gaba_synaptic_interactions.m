% sanity_checks.m

clear

gpu = 0;

% fieldstosave =  {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r', ...
 % 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a' , 'g_GABA_soma','g_GABA_dend', 'g_AMPA', 'Ca2_soma'};
pert = []; I_app = [];

experiment = 2;
switch experiment
	case 0 % gaba soma in 2D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend', 'g_GABA_soma'};


		simtime =500;
		netsize = [1 3 3];s
		gap = .1;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		
		cell_parameters = createDefaultNeurons(9);
		cell_parameters.g_CaL = ones(1,prod(netsize))*.4;
		cell_parameters.g_CaL(centercell) = .6;	

		g_CaL = ones(1,prod(netsize))*.4;
		g_CaL(centercell) = .6;
		
		pert.mask  	   {1} = [0 0 0 0 1 0 0 0 0 ]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [100];
		pert.duration  {1} = 2;
		% pert.type	   {1} = 'gaba_dend'; % gaba_soma ; ampa_noise ; ampa_soma
		pert.type	   {1} = 'gaba_soma'; % gaba_soma ; ampa_noise ; ampa_soma


		pert.mask  	   {2} = [0 0 0 0 1 0 0 0 0 ]';
		pert.amplitude {2} = 1;
		pert.triggers  {2} = [50 120 200];
		pert.duration  {2} = 2;
		pert.type	   {2} = 'ampa_soma'; % gaba_soma ; ampa_noise ; ampa_soma

	
	case 1 % gaba soma in 2D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend'};


		simtime =1000;
		netsize = [1 3 3];
		% gap = .1;
		gap = 0;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		
		cell_parameters = createDefaultNeurons(9);
		cell_parameters.g_CaL = ones(1,prod(netsize))*.4;
		cell_parameters.g_CaL(centercell) = .6;	

		pert.mask  	   {1} = [0 0 0 0 1 0 0 0 0 ]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [100 :10: 200 ];
		pert.duration  {1} = 1;
		% pert.type	   {1} = 'gaba_dend'; % gaba_soma ; ampa_noise ; ampa_soma
		pert.type	   {1} = 'gaba_soma'; % gaba_soma ; ampa_noise ; ampa_soma
		

	case 2 % gaba soma and dendrite in 2D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend'};


		simtime = 2000;
		netsize = [1 3 3];
		% gap = .05;
		gap = .0;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		
		cell_parameters = createDefaultNeurons(9);
		cell_parameters.g_CaL = ones(1,prod(netsize))*.6;
		cell_parameters.g_CaL(centercell) = 1.1;	

		cell_parameters.gbar_gaba_dend = ones(1,prod(netsize))*2;
		cell_parameters.gbar_gaba_soma = ones(1,prod(netsize))*1;

		cell_parameters.V_gaba_soma = ones(1,prod(netsize))*-70;



		pert.mask  	   {1} = [1 1 1 1 1 1 1 1 1]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [200:5:250 800 1500:5:1530];
		pert.duration  {1} = 5;
		pert.type	   {1} = 'gaba_soma';
		
		pert.mask  	   {2} = [1 1 1 1 1 1 1 1 1 ]';
		pert.amplitude {2} = 1;
		pert.triggers  {2} = [800:805 1505:5:1555];
		pert.duration  {2} = 5;
		pert.type	   {2} = 'gaba_dend'; 
		


end


[transients] = IOnet('networksize', netsize,'perturbation', pert ,'appCurrent',I_app,'time',simtime,'cell_parameters', cell_parameters ,'W', W ,'ou_noise', [0 0 0 0],'to_report', to_report,'gpu', gpu);

plot(transients.networkHistory.V_soma')









