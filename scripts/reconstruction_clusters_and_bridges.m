% demo_clusters.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear
rng(0,'twister') % random seed

steadystate_time = 1000; %ms
simtime  = 5000; %ms
delta = .025;
gpu = 1;



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

nconns_curlies = 8;
nconns_bridges = 6;

gap_curlies = .1;
gap_bridges = .1;
plotconn = 1;
normalize = 1;

nconns_curlies = 5;
nconns_bridges = 4;
gap_curlies = .05;
gap_bridges = .05;
plotconn = 1;
normalize = 1;



load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);

if not(exist('curlies'))
	% curlies:
	% create a network with distance based connectivity for close by connections
	% this network is clusterized with about 20cells per cluster, according to a k-means algo.
	curlies = createW('3d_reconstruction', [], 4*40, 1, 0, 1, [], nconns_curlies, somatapositions,1,[1 20 1 0]);

	% create a network with distance based connectivity for further apart cells: bridges
	% these cells are not bound to specific clusters.
	bridges = createW('3d_reconstruction', [], 8*40, 1, 0, 1, [], nconns_bridges, somatapositions,1,[1 20 0 1]);

	% define the indices of 10% of the cells, these will be bridges
	bc = randperm(noneurons); % randomly permute cell indices
	z = zeros(noneurons,1) ; % initialize index vector
	z(bc(1:round(.1*noneurons))) = 1; % make 10% of the cells == bridges
	bc =z;

	% remove from curlie adjacency matrix all of those that will become bridges
	curlies.W = bsxfun(@times, curlies.W, ~z); 
	curlies.W = bsxfun(@times, curlies.W, ~(z'))*gap_curlies; % multiply by the 'unitary' conductance
	cstats = connectivity_statistics(bridges);
	curlies.stats = cstats.stats ;

	% remove connections from curlies to bridges from bridge adjacency matrix
	bridges.W = bsxfun(@times, bridges.W, z);
	% create bridge cells connectivity 
	bridges.W = (bridges.W+bridges.W')*gap_bridges;
	bstats = connectivity_statistics(bridges);
	bridges.stats = bstats.stats ;

	bridg_curlies.coords = curlies.coords;
	
	
	bridg_curlies.W = curlies.W + bridges.W;
	bridg_curlies.stats = connectivity_statistics(bridg_curlies);
	bridg_curlies.stats.clusters = curlies.stats.clusters;

	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridg_curlies.stats.clusters)

	brick = createW('3d_reconstruction', [], 4*40, 1, 0, 1, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
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
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized3');
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

% numberofmasks = 10; 
onset_of_stim = [1005:2:1025];

% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

% pert.mask     {1} =  create_input_mask(netsize, 'dist_to_center','radius',2, 'synapseprobability', 1,'plotme',1);
pert.mask     {1} =  [curlies.stats.clusters==5] | [curlies.stats.clusters==10] | [curlies.stats.clusters==20];
pert.amplitude{1} = 1;
pert.triggers {1} = onset_of_stim;
pert.duration {1} = 10;
pert.type	  {1} = 'gaba_soma';
% pert.type	  {1} = 'ampa';





% [===========================================================================================================]
 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
 if ~exist('st_st','var')
	disp('calculating transients')

	 st_st = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 st_st.note = 'curlies and bridges'
end

% 	% st_st.Plist = Plist;
% end


% [=================================================================]
%  GABA
% [=================================================================]

% BRIDGES AND CURLIES WITH PERTURBATION
if 1
	 sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{1}.note = 'curlies and bridges'
	 sim{1}.W = bridg_curlies;

	sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
	sim{1}.networkHistory.I_cx36 = single(sim{1}.networkHistory.V_soma);
	sim{1}.networkHistory.backgroundnoise = [];

end

% ONLY CURLIES
if 1
	sim{2} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{2}.note = 'only curlies'
	 sim{2}.W = curlies;

	sim{2}.networkHistory.V_soma = single(sim{2}.networkHistory.V_soma);
	sim{2}.networkHistory.I_cx36 = single(sim{2}.networkHistory.V_soma);
	sim{2}.networkHistory.backgroundnoise = [];

end


% DISCONNECTED NETWORK
if 1
	 sim{3} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', curlies.W*0 ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{3}.note = 'gap is zero'
	 sim{3}.W = curlies;
	 

	sim{3}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
	sim{3}.networkHistory.I_cx36 = single(sim{3}.networkHistory.V_soma);
	sim{3}.networkHistory.backgroundnoise = [];

end



% DISCONNECTED NETWORK
if 1
	 sim{4} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{4}.note = 'brick connectivity'
	 sim{4}.W = brick;
	 

	sim{4}.networkHistory.V_soma = single(sim{3}.networkHistory.V_soma);
	sim{4}.networkHistory.I_cx36 = single(sim{3}.networkHistory.V_soma);
	sim{4}.networkHistory.backgroundnoise = [];

end


% [=================================================================]
%  spontaneous
% [=================================================================]



eval(['save clusters_curlies_bridges_'  date ' -v7.3'])








