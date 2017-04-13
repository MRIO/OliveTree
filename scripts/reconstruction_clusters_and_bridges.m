% demo_clusters.m


%  4 simulations
% : stimulate group 41 and one of its neighbors
% with and without bridges
% looking at kuramoto and mystical neighbor


% [================================================]
% 		 simulation parameters
% [================================================]
% clear
rng(0,'twister') % random seed

steadystate_time = 1000; %ms
simtime  = 10000; %ms
delta = .025;
gpu = 1;

% [================================================]
% simulations to perform
% [================================================]

frombrick_to_clusters = 0;
bridge_conductance_pspace = 0;
conjuctive_stimulation = 1
	cluster1 = 41; cluster2 = 34;
bridges_and_curlies_with_gaba = 0;
maskstimulation = 1;
nostimulation = 1;

% [================================================]
% variables to report
% [================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

% variables to store
to_report = vsoma;


% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)

% nconns_curlies = 8;
% nconns_bridges = 6;

% gap_curlies = .1;
% gap_bridges = .1;
% plotconn = 1;
% normalize = 1;

nconns_curlies = 5;
nconns_bridges = 5;
cells_in_cluster = 20;
gap_curlies = .05;
gap_bridges = .05;
plotconn = 0;
normalize = 1;



load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);

if not(exist('curlies'))
	% curlies:
	% create a network with distance based connectivity for close by connections
	% this network is clusterized with about 20cells per cluster, according to a k-means algo.
	curlies = createW('3d_reconstruction', [], 4*40, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[1 cells_in_cluster 1 0]);

	% create a network with distance based connectivity for further apart cells: bridges
	% these cells are not bound to specific clusters.
	bridges = createW('3d_reconstruction', [], 8*40, 1, 0, plotconn, [], nconns_bridges, somatapositions,1,[1 cells_in_cluster 0 1]);


	% define the indices of 10% of the cells, these will be bridges
	bc = randperm(noneurons); % randomly permute cell indices
	z = zeros(noneurons,1) ; % initialize index vector
	z(bc(1:round(.1*noneurons))) = 1; % make 10% of the cells == bridges
	bc =z;
	bridge_idx = find(bc);

	% remove from curlie adjacency matrix all of those that will become bridges
	curlies.W = bsxfun(@times, curlies.W, ~z);
	curlies.W = bsxfun(@times, curlies.W, ~(z')); % multiply by the 'unitary' conductance
	curlies_bu = curlies; %_bu -> binary undirected
	curlies.W = curlies.W*gap_curlies;
	% curlies.stats.clusters(bridge_idx) = 0;
	
	cstats = connectivity_statistics(bridges);
	curlies.stats = cstats.stats ;

	% remove connections from curlies to bridges from bridge adjacency matrix
	bridges.W = bsxfun(@times, bridges.W, z);
	% create bridge cells connectivity 
	bridges.W = (bridges.W+bridges.W');
	bridges_bu = bridges; % _bu -> binary undirected
	bridges.W = bridges.W*gap_bridges;

	bstats = connectivity_statistics(bridges);
	bridges.stats = bstats.stats ;

	bridg_curlies.coords = curlies.coords;

	bridg_curlies.W = curlies.W + bridges.W;
	bridg_curlies.stats = connectivity_statistics(bridg_curlies);
	bridg_curlies.stats.clusters = curlies.stats.clusters;
	% bridg_curlies.stats.clusters(bridge_idx) = 0;

	clusteridx = bridg_curlies.stats.clusters;
	clusteridx(logical(bc)) = 70;
	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), clusteridx)

	brick = createW('3d_reconstruction', [], 8*40, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
	brick_bu = brick;
	brick.W = brick.W*gap_curlies;

end



% Wcluster150 = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], plotconn, [1 150 .9 .01],1);


% [=================================================================]
%  create cells
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
% netsize = [3 15 15];

% def_neurons = createDefaultNeurons(noneurons,'celltypes','param_sweep');
% def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized2');
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized3') 
% def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized4');

% randomized2 = 
% neurons.g_CaL = linspace(.5, 1, noneurons);



% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% gnoise = [.2 .3 0 5];
gnoise = [0 0 0 0];
sametoall = 0.05;


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]

% create overlapping ampa masks
if conjuctive_stimulation
	clusters = bridg_curlies.stats.clusters;
	bridges_in_cluster = single(bc .* clusters==41);
			neighbors_to_bridge = find(bridges_in_cluster'*(bridg_curlies.W>0));
			their_cluster = unique(clusters(neighbors_to_bridge));
			targeted_cluster_cells = ismember(clusters, their_cluster).*clusters;
			% plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), targeted_cluster_cells);
			plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), clusters==41 | clusters==34);
end

% numberofmasks = 10; 
onset_of_stim = [3005:5:3025];

% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

% pert.mask     {1} =  create_input_mask(netsize, 'dist_to_center','radius',2, 'synapseprobability', 1,'plotme',1);
pert.mask     {1} =  [curlies.stats.clusters==41];
pert.amplitude{1} = 1;
pert.triggers {1} = onset_of_stim;
pert.duration {1} = 10;
% pert.type	  {1} = 'gaba_soma';
pert.type	  {1} = 'ampa';





% [===========================================================================================================]
 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
 if ~exist('st_st','var')
	disp('calculating transients')

	 st_st = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W*0 ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 st_st.note = 'curlies and bridges'
end

% 	% st_st.Plist = Plist;
% end


if conjuctive_stimulation
	
			
pert.mask     {1} =  [curlies.stats.clusters==41];
			
 	sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{1}.note = '41_with_bridges'
	sim{1}.W = bridg_curlies;
	sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
	sim{1}.networkHistory.I_cx36 = [];
	sim{1}.networkHistory.backgroundnoise = [];


pert.mask     {1} =  [curlies.stats.clusters==41 | curlies.stats.clusters==34];
						
sim{2} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{2}.note = '41&34_with_bridges'
	sim{2}.W = bridg_curlies;
	sim{2}.networkHistory.V_soma = single(sim{2}.networkHistory.V_soma);
	sim{2}.networkHistory.I_cx36 = [];
	sim{2}.networkHistory.backgroundnoise = [];


pert.mask     {1} =  [curlies.stats.clusters==41];
			
 	sim{3} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{3}.note = '41_with_curlies'
	sim{3}.W = curlies;
	sim{3}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
	sim{3}.networkHistory.I_cx36 = [];
	sim{3}.networkHistory.backgroundnoise = [];


pert.mask     {1} =  [curlies.stats.clusters==41 | curlies.stats.clusters==34];
						
	sim{4} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{4}.note = '41&34_with_curlies'
	sim{4}.W = curlies;
	sim{4}.networkHistory.V_soma = single(sim{4}.networkHistory.V_soma);
	sim{4}.networkHistory.I_cx36 = [];
	sim{4}.networkHistory.backgroundnoise = [];

	sim{5} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{5}.note = 'nostim_curlies'
	sim{5}.W = curlies;
	sim{5}.networkHistory.V_soma = single(sim{5}.networkHistory.V_soma);
	sim{5}.networkHistory.I_cx36 = [];
	sim{5}.networkHistory.backgroundnoise = [];

	sim{6} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	sim{6}.note = 'nostim_curlies'
	sim{6}.W = curlies;
	sim{6}.networkHistory.V_soma = single(sim{6}.networkHistory.V_soma);
	sim{6}.networkHistory.I_cx36 = [];
	sim{6}.networkHistory.backgroundnoise = [];



		eval(['save curlies_bridges_stim_pair'  date ' -v7.3'])


end




if bridge_conductance_pspace
	s = 0;
	for bridge_conductance = [eps 0.01:0.01:0.1]
		s = s +1;

		bridg_conduct.W = curlies_bu.W*gap_curlies + bridges_bu.W*bridge_conductance;
		bridg_conduct.stats = connectivity_statistics(bridg_conduct);
		bridg_conduct.stats.clusters = curlies.stats.clusters;

		note = 'bridge conductance pspace';
	 	sims{s} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_conduct.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall, 'displaytext', [note '_' num2str(s)]);
	 	sims{s}.note = note;
	 	sims{s}.W = bridg_conduct;
	 	sims{s}.bridge_conductance = num2str(bridge_conductance);
	 end
	 eval(['save bridge_conductance_pspace'  date ' -v7.3'])
	 clear sims

end






if frombrick_to_clusters
%%% TODO: MAKE SURE THAT THE TOTAL GAP CONDUCTANCE IN THE NETWORK IS IDENTICAL

	s = 0;
	for alpha_W = [0:.25:1]
		s = s +1;

		mixture.W = alpha_W*curlies.W + (1-alpha_W)*brick.W;

		mixture.stats = connectivity_statistics(mixture);
		mixture.stats.clusters = curlies.stats.clusters;
		note = 'brick to clusters';

	 	sims{s} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', mixture.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall, 'displaytext' , [note '_' num2str(s)] );
	 	sims{s}.note = ['from brick to clusters, alpha = ' num2str(alpha_W)];
	 	sims{s}.W = mixture;
	 	
	 end
	 eval(['save brick_to_clusters_'  date ' -v7.3'])
	 clear sims
end





% [=================================================================]
%  GABA
% [=================================================================]

% BRIDGES AND CURLIES WITH PERTURBATION
if bridges_and_curlies_with_gaba 
	 sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{1}.note = 'curlies and bridges'
	 sim{1}.W = bridg_curlies;

	sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
	sim{1}.networkHistory.I_cx36 = single(sim{1}.networkHistory.I_cx36);
	sim{1}.networkHistory.backgroundnoise = [];


% ONLY CURLIES
	sim{2} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{2}.note = 'only curlies'
	 sim{2}.W = curlies;

	sim{2}.networkHistory.V_soma = single(sim{2}.networkHistory.V_soma);
	sim{2}.networkHistory.I_cx36 = single(sim{2}.networkHistory.I_cx36);
	sim{2}.networkHistory.backgroundnoise = [];


% DISCONNECTED NETWORK
	 sim{3} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W*0 ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{3}.note = 'gap is zero'
	 sim{3}.W = curlies.W*0;
	 

	sim{3}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
	sim{3}.networkHistory.I_cx36 = single(sim{3}.networkHistory.I_cx36);
	sim{3}.networkHistory.backgroundnoise = [];




% BRICK NETWORK
	 sim{4} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{4}.note = 'brick connectivity'
	 sim{4}.W = brick;
	 

	sim{4}.networkHistory.V_soma = single(sim{4}.networkHistory.V_soma);
	sim{4}.networkHistory.I_cx36 = single(sim{4}.networkHistory.I_cx36);
	sim{4}.networkHistory.backgroundnoise = [];


	eval(['save curlies_bridges_'  date ' -v7.3'])


end



if maskstimulation
	pert.mask{1} = create_input_mask(somatapositions, 'reconstruction','synapseprobability',.4, 'radius', 100, 'plotme',1)
	
	sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{1}.note = 'curlies random mask';
		sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
		sim{1}.W = curlies;

	sim{2} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{2}.networkHistory.V_soma = single(sim{2}.networkHistory.V_soma);
		sim{2}.note = 'bridge random mask'
		sim{2}.W = bridg_curlies;

	sim{3} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{3}.note = 'bridge no stim'
		sim{3}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
		sim{3}.W = bridg_curlies;

	sim{4} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{4}.note = 'brick no stim'
		sim{4}.networkHistory.V_soma = single(sim{4}.networkHistory.V_soma);
		sim{4}.W = brick;

		eval(['save curlies_bridges_randmaskstim'  date ' -v7.3'])

end



if nostimulation
	simtime = 15000;
	pert = [];
	
	sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{1}.note = 'curlies random mask';
		sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
		sim{1}.W = curlies;

	sim{2} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{2}.networkHistory.V_soma = single(sim{2}.networkHistory.V_soma);
		sim{2}.note = 'bridge random mask'
		sim{2}.W = bridg_curlies;

	sim{3} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{3}.note = 'bridge no stim'
		sim{3}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
		sim{3}.W  = bridg_curlies;

	sim{4} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
		sim{4}.note = 'brick no stim'
		sim{4}.networkHistory.V_soma = single(sim{4}.networkHistory.V_soma);
		sim{4}.W = brick;

		eval(['save curlies_bridges_nostim'  date ' -v7.3'])

end


