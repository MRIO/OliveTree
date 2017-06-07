% Input correlations and center peak sharpnesss?

% pspace tests for distributions of ampa input over a wee time

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

cell_function = 'vanilla'; % 'devel'
steadystate_time = 10; %ms
simtime  = 1000; %ms
delta = .025;
gpu = 1;

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
selection = {'V_soma' 'I_cx36' 'I_CaH' 'Calcium_r'};

% variables to store
to_report = selection ;


% [================================================]
% 		 create default neurons
% [================================================]
noneurons = 20;
rng(0,'twister')
gap_neurons = createDefaultNeurons(noneurons,'celltypes', 'randomized','gapcompensation', 1);
gap_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; %linspace(0,.3,noneurons);
gap_neurons.gbar_gaba_soma = ones(noneurons,1)*.5; %ones(noneurons,1)*.3;

rng(0,'twister')
def_neurons = createDefaultNeurons(noneurons,'celltypes', 'randomized','gapcompensation', 0);
def_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; %linspace(0,.3,noneurons);
def_neurons.gbar_gaba_soma = ones(noneurons,1)*.5; %ones(noneurons,1)*.3;

netsize = [1 noneurons 1];

% [================================================]
% 		gap connections
% [================================================]
 W = zeros(noneurons);
 % W = createW(noneurons);

gaps = [0.05];
noconnections = 10;
symmetrize = 1;
randomize = 0;
W = createW('all to all', netsize, [],   1,       randomize,         0,        0,       noconnections, [],symmetrize, [0 0 0 0], 1);
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)


% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% gnoise = [0 0 0 0];
gnoise = [0 0 0 4]; sametoall = 0;
% [================================================]
%  Ornstein Uhlenbeck Perturbation over masks
% [================================================]

% create overlapping ou-masks

numberofmasks = 20;

onset_of_stim = 100;
stim_dur      = 2000;
offset_stim   = 100;
% synapseprobability = .3;

th =	 1/5 ; % decay time parameter
mu = 	 0 ; % pA
sig = 	 2 ; % pA
mix =    1;


	for nm = 1:numberofmasks

		pert.mask  	  {nm} = zeros(noneurons,1);
		pert.mask  	  {nm}(nm) = 1;

		% pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', synapseprobability);
		pert.amplitude{nm} = 1;
		pert.triggers {nm} = onset_of_stim + (nm-1)*offset_stim;
		pert.duration {nm} = stim_dur;
		pert.type	  {nm} = 'ou_noise';

		pert.param{nm}(1)  = th  ;
		pert.param{nm}(2)  = mu  ;
		pert.param{nm}(3)  = sig ;
		pert.param{nm}(4)  = mix;


	end 


% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
end


% [===========================================================================================================]
 ns = 0;
 for g = gaps
	 	ns = ns +1;


		if g>0
			neurons = gap_neurons;
		else
			neurons = def_neurons;
		end


		   [simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


end
% [===========================================================================================================]

% save simresults_ounoisemask simresults

spks{1} = spikedetect(simresults{1}, 0,0);


% netsize = [1 112 1];
% numberofmasks = 15; 

% onset_of_stim = 100;
% stim_dur      = 300;
% offset_stim   = 1;

% th =	 3  ; % decay time parameter
% mu = 	 0  ; % pA
% sig = 	 10 ; % pA
% mix = 	 .2 ;

% for nm = 1:numberofmasks

% 			pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', .25);
% 			pert.amplitude{nm} = 1;
% 			pert.triggers {nm} = onset_of_stim + (nm-1)*offset_stim;
% 			pert.duration {nm} = stim_dur;
% 			pert.type	  {nm} = 'ou_noise';

% 			pert.param{nm}(1)  = th  ;
% 			pert.param{nm}(2)  = mu  ;
% 			pert.param{nm}(3)  = sig ;
% 			pert.param{nm}(4)  = mix ;

% 			pert.rng{nm} = rng(nm);

% 		end 

% tic
% curr_noise = zeros(112,50000);

% for t = 1:5000
% [g_ampa_soma g_gaba_dend g_gaba_soma Ca2_soma curr_noise(:,t)] = ...
%             apply_perturbation(pert,t,1000, zeros(112,1), 0, 0, 0);
% end
% toc


