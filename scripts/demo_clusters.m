% demo_clusters.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear
rng(0,'twister') % random seed

steadystate_time = 1000; %ms
simtime  = 10000; %ms
delta = .025;
gpu = 1;


activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

% variables to store
to_report = gapcur;

% [=================================================================]
%  create cells
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
netsize = [5 10 10];
noneurons = prod(netsize);

def_neurons = createDefaultNeurons(noneurons,'celltypes','param_sweep');
% neurons.g_CaL = linspace(.5, 1, noneurons);


% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)
W = createW('all to all', netsize, [], 1, 1, 1, [], 15, [], 0, [1 20 1 .001]);
gaps = .05;

% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
gnoise = [.2 1 0 0];
% gnoise = [0 0 0 0];
sametoall = 0.2;


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]

% create overlapping ampa masks

numberofmasks = 10; 
stim_interval = [25 50 100 250 500];
onset_of_stim = 8000;
n_of_pulses = 10;

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
	% st_st.Plist = Plist;
end



% [=================================================================]
%  do it 
% [=================================================================]

% 'tempState', simresults.lastState ,
if exist('sim','var')
	sim = sim+1;
else
	sim = 2; V = []; simR{sim-1}.lastState = st_st.lastState;
end

 simR{sim} = IOnet('tempState', simR{sim-1}.lastState, 'cell_parameters', neurons, ...
	   	'networksize', netsize ,'time',simtime ,'W', W.W*gap ,'ou_noise', gnoise , ...
	   	'to_report', to_report ,'gpu', gpu , ...
	   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

simR{sim}.clusters = W.stats.clusters;

[bla sortedclusters] = sort(W.stats.clusters);


V = [V simR{sim}.networkHistory.V_soma(sortedclusters,:)];

clf, imagesc(V);




% [===========================================================================================================]
 ns = 0;
 for g = gaps
 	for si = stim_interval 
	 	ns = ns +1;

 	
 		for nm = 1:numberofmasks

			pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', .25);
			pert.amplitude{nm} = 1;
			pert.triggers {nm} = onset_of_stim + round(si*rand) + cumsum(poissrnd(si,n_of_pulses,1)) ;
			pert.duration {nm} = 1;
			pert.type	  {nm} = 'gaba_soma';
			
		end 

		if g > 0
			gaussnoise(2) = 10;
		end


	   [simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, ...
	   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W*g ,'ou_noise', gnoise , ...
	   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
	   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

	end
end
% [===========================================================================================================]
