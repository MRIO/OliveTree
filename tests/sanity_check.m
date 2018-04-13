 % sanity_checks.m


gpu = 1;

cell_function = 'vanilla'; % unless specified below

availablefieldstosave =  {'V_soma', 'I_cx36', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r',...
 'I_CaH', 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a', 'Ca2_soma'};
to_report =  {'V_soma','V_dend', 'Ca2Plus', 'I_cx36'};

pert = []; I_app = [];

dt = .02;

iapp = false;

experiment = 'tycho'
experiment = 'devel against vanilla'

switch experiment

	case 'one_cell'

		dt = 0.01;
		simtime = 1000;
		netsize = [1 1 1];
		gap = 0;

		W = 0;
		cell_parameters = createDefaultNeurons(1);
		
		I_app = zeros(1, simtime*(1/dt));
		I_app(1,(500*(1/dt):700*(1/dt))) = -5;% uA/cm^2 
		gnoise = [0 0 0 0]; % 



	case 'tycho' % two cells 
		dt = 0.025;
		simtime = 200;
		netsize = [1 5 5];
		gap = .1;
		noneurons = prod(netsize);

		W = rand(noneurons)*gap;
		W(find(eye(noneurons))) = 0;
		cell_parameters = createDefaultNeurons(prod(netsize));

		cell_parameters.g_CaL = .4 + rand(1,noneurons)'*.4;
		
		I_app = zeros(noneurons, simtime*(1/dt));
		I_app(5,(500*(1/dt):700*(1/dt))) = -5;% uA/cm^2 
		gnoise = [3 5 0 0]; % 

	
	case 'three cells' % three cells
		dt = 0.01;
		netsize = [1 3 1];
		simtime  = 300;
		gap = 0.2;
		W = (ones(3)-eye(3))'*gap;
		W(3,:) = 0;W(:,3) = 0;

		cell_parameters = createDefaultNeurons(3);
		cell_parameters.g_CaL = [.4 .4 .4]; 

		I_app = zeros(3, simtime*(1/dt));
		I_app(1,(200*(1/dt):220*(1/dt))) = -5;% uA/cm^2 
		gnoise = [0 0 0 0];
		cell_function = 'devel'; % experimental!!!!!

	
	case '2Dcenter' % the cell surrounded 2D
		dt = 0.01;
		simtime = 1000;
		netsize = [1 3 3];
		gap = .05;
		W = zeros(9);
		W(5,[1 2 3 4 6 7 8 9]) = 1*gap;
		W = W+W';
		cell_parameters = createDefaultNeurons(9);
		cell_parameters.g_CaL = ones(1,9)*.9;
		% g_CaL(5) = 1.;
		I_app = zeros(9, simtime*(1/dt));
		I_app(5,(500*(1/dt):520*(1/dt))) = 10;% uA/cm^2 
		gnoise = [3 1 0 0]; % 
		cell_function = 'devel'; 


	case '3Dcenter' % the cell surrounded 3D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend','curr_noise'};
		simtime = 1000;
		netsize = [3 3 3];
		gap = .5;
		centercell = 14;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		cell_parameters = createDefaultNeurons(prod(netsize));
		cell_parameters.g_CaL = ones(1,prod(netsize))*1.1;
		cell_parameters.g_CaL(centercell) = .4;
		
		I_app = zeros(prod(netsize), simtime*20);
		I_app(centercell,(500:600)*20) = 20;
		gnoise = [0 0 0 0];


	case 'gaba' % gaba soma in 2D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend', 'g_GABA_soma'};

		simtime = 1000;
		netsize = [1 3 3];
		cell_parameters = createDefaultNeurons(prod(netsize));
		cell_parameters.gbar_gaba_soma = .5*ones(prod(netsize),1);
		cell_parameters.gbar_gaba_dend = .5*ones(prod(netsize),1);
		gap = .01;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		g_CaL = ones(1,prod(netsize))*.4;
		g_CaL(centercell) = 1.1;
		
		pert.mask  	   {1} = [0 0 1 1 1 1 1 0 0 ]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [100:10:150];
		pert.duration  {1} = 5;
		pert.type	   {1} = 'gaba_soma'; % gaba_soma ; ampa_noise ; ampa_soma

		pert.mask  	   {2} = [0 0 0 0 1 0 0 0 0 ]';
		pert.amplitude {2} = 1;
		pert.triggers  {2} = [500:10:600];
		pert.duration  {2} = 3;
		pert.type	   {2} = 'gaba_dend'; % gaba_soma ; ampa_noise ; ampa_soma

		gnoise = [0 0 0 0];

	case 'gaba_soma_dend' % gaba soma and dendrite in 2D
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend', 'g_GABA_soma','g_GABA_dend'};


		simtime = 2000;
		netsize = [1 3 3];
		gap = .05;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';

		cell_parameters = createDefaultNeurons(prod(netsize));
		cell_parameters.g_CaL = ones(1,prod(netsize))*.4;
		cell_parameters.g_CaL(centercell) = 1.1;
		
		pert.mask  	   {1} = [1 1 1 1 1 1 1 1 1]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [500 1500];
		pert.duration  {1} = 2;
		pert.type	   {1} = 'gaba_soma'; % gaba_soma ; ampa_noise ; ampa_soma
		
		pert.mask  	   {2} = [1 1 1 1 1 1 1 1 1 ]';
		pert.amplitude {2} = 1;
		pert.triggers  {2} = [500 1000 ];
		pert.duration  {2} = 2;
		pert.type	   {2} = 'gaba_dend'; % gaba_soma ; ampa_noise ; ampa_soma
		
		gnoise = [0 0 0 0];

	case 'two_cells' % two cells - current input
		to_report = availablefieldstosave;
		% to_report = {'V_soma','Ca2Plus', 'I_CaH', 'V_axon', 'Calcium_r','Hcurrent_q', 'Potassium_s'};
		netsize = [1 2 1];
		simtime  = 500;
		gap = .0;
		W = [0 1; 1 0]*gap;
		I_app = zeros(2, simtime*20);
		I_app(2,(100:350)*20) = -10; % nAmpere 20/dt [nA/cm^2] 
		I_app(2,(351:360)*20) =  21; % nAmpere 20/dt [nA/cm^2] 
		g_CaL = [.4 1.2];

		gnoise = [0 0 0 0];

		cell_parameters = createDefaultNeurons(2);
		iapp = true;


	case 'continuations'
		 % Sanity checks for continuating simulations from saved last state.

		netsize = [1 2 1];
		simtime  = 800;
		W = [0 1; 1 0];
		dt = .05;

		% overwrite defaults
		cell_parameters = createDefaultNeurons(2);
		cell_parameters.g_CaL = [0.4 1.1]; 

		% create initial condition
		rndState = initNetState(prod(netsize),gpu, 0);
		rndState.V_soma = [-100 -80]';
		st_st = IOnet('delta', dt, 'networksize', netsize, 'time',simtime ,'W', W*0 ,'to_report',{'V_soma'},'gpu', gpu, 'tempState', rndState,'cell_parameters', cell_parameters);

	% to continue the network and see if it stitches well
	    simtime = 200;
	    cont  = IOnet('tempState', st_st.lastState ,'cell_parameters', cell_parameters, 'networksize', netsize,'time',simtime ,'W', W*0  ,'to_report',to_report,'gpu', gpu);
	    plot([st_st.networkHistory.V_soma cont.networkHistory.V_soma]')

	    % break


	case 'somatic_ampa' % ampa in soma
		to_report = {'V_soma','I_cx36', 'Ca2Plus', 'Ca2_soma' 'V_dend', 'g_GABA_soma'};

		simtime = 300;
		netsize = [1 3 3];
		cell_parameters = createDefaultNeurons(prod(netsize));
		cell_parameters.g_CaL = ones(1,prod(netsize))*.5;
		cell_parameters.gbar_ampa_soma = linspace(.01, .1, prod(netsize));
		gap = .01;
		centercell = 5;
		W = zeros(prod(netsize)); 
		allcells = 1:prod(netsize);
		othercells = setdiff(allcells, centercell);
		W(centercell, othercells) = 1*gap;
		W = W+W';
		
		
		pert.mask  	   {1} = [1 1 1 1 1 1 1 1 1]';
		pert.amplitude {1} = 1;
		pert.triggers  {1} = [100];
		pert.duration  {1} = 10;
		pert.type	   {1} = 'ampa'; % gaba_soma ; ampa_noise ; ampa_soma

		gnoise = [0 0 0 0];

	case 'devel against vanilla'

		to_report = availablefieldstosave;
		% to_report = {'V_soma','Ca2Plus', 'I_CaH', 'V_axon', 'Calcium_r','Hcurrent_q', 'Potassium_s'};
		netsize = [1 2 1];
		simtime  = 200;
		gap = 0;
		W = [0 1; 1 0]*gap;
		g_CaL = [.4 1.2];
		dt = 0.05;

		gnoise = [0 0 0 0];

		cell_parameters = createDefaultNeurons(2);
		cell_parameters.g_CaL = [.4 1.2];
		

	
		[tr_vanilla] = IOnet('cell_function', 'vanilla', 'networksize', netsize, 'cell_parameters', cell_parameters,  ...
	'perturbation', pert ,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', to_report ,'gpu', gpu, 'delta', dt);

		[tr_devel]   = IOnet('cell_function', 'devel', 'networksize', netsize, 'cell_parameters', cell_parameters,  ...
	'perturbation', pert ,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', to_report ,'gpu', gpu, 'delta', dt);

		replayResults(tr_vanilla, 'plotallfields',1)
		replayResults(tr_devel, 'plotallfields',1)


		plot(tr_devel.networkHistory.V_soma', 'b'); hold on
		plot(tr_vanilla.networkHistory.V_soma', 'r');
		plot([tr_devel.networkHistory.V_soma-tr_vanilla.networkHistory.V_soma]', 'g');
		
		


end


% [transients] = IOnet('networksize', netsize,'perturbation', pert ,'appCurrent',I_app,'time',simtime,'g_CaL', g_CaL ,'W', W ,'ou_noise', gnoise ,'to_report', to_report,'gpu', gpu);

if iapp
[transients] = IOnet('cell_function', cell_function, 'networksize', netsize, 'cell_parameters', cell_parameters,  ...
	'perturbation', pert ,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', to_report ,'gpu', gpu, 'delta', dt);
else
	
[transients] = IOnet('cell_function', cell_function, 'networksize', netsize, 'cell_parameters', cell_parameters,  ...
	'perturbation', pert ,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', to_report ,'gpu', gpu, 'delta', dt);
end

replayResults_3(transients)

replayResults(transients, 'plotallfields',1)






