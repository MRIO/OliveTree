% demo_clusters.m


%  4 simulations
% : stimulate group 41 and one of its neighbors
% with and without bridges
% looking at kuramoto and mystical neighbor


% [================================================]
% 		 simulation parameters
% [================================================]
% clear

steadystate_time = 1000; %ms
simtime  = 3000; %ms
delta = .025;
gpu = 1;

if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
thisseed.Seed

% [================================================]
% simulations to perform
% [================================================]

frombrick_to_clusters = 0;
bridge_conductance_pspace = 0;
conjuctive_stimulation = 1;
	stim1 = 'gaba_dend'; stim2 = 'ampa_dend';
bridges_and_curlies_with_gaba = 0;
maskstimulation = 0;
two_maskstimulation = 0;
nostimulation = 0;

% [================================================]
% variables to report
% [================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

% variables to store
to_report = gapcur;


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
cells_in_cluster = 20;  % upperbound estimate from: % Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.
gap_curlies = .05;
gap_bridges = .05;
plotconn = 1;
normalize = 1;



load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);


%     __    __
%    / /_  / /___ _
%   / __ \/ / __ `/
%  / /_/ / / /_/ /
% /_.___/_/\__,_/


if not(exist('curlies'))
	% curlies:
	% create a network with distance based connectivity for close by connections
	% this network is clusterized with about 20cells per cluster, according to a k-means algo.

	% Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.

	% radius within which we allow connections
	% between curlies:
	median_soma_distance = 20;
	rad_cur = median_soma_distance * 6;
	% between bridges:
	rad_bri = median_soma_distance * 12;


	

	curlies = createW('3d_reconstruction', [], rad_cur, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[1 cells_in_cluster 1 0]);

	% create a network with distance based connectivity for further apart cells: bridges
	% these cells are not bound to specific clusters.
	bridges = createW('3d_reconstruction', [], rad_bri, 1, 0, plotconn, [], nconns_bridges, somatapositions,1,[1 cells_in_cluster 0 1]);


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

	brick = createW('3d_reconstruction', [], rad_bri, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
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
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized3', 'rng', thisseed) 
% def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized4');

% randomized2 = 
% neurons.g_CaL = linspace(.5, 1, noneurons);



% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% ounoise_params = [.2 .3 0 5];
ounoise_params = [0 0 0 double(thisseed.Seed)];
sametoall = 0.05;


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]


pert.mask     {1} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-50, -50, 0], 'synapseprobability', 1,'plotme',0)
pert.amplitude{1} = 1;

pert.duration {1} = 10;
pert.type	  {1} = 'gaba_dend';

pert.mask     {2} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-50, 50, 0], 'synapseprobability', .85,'plotme',0) % probability adjusted to make number of cells in cluster match
pert.amplitude{2} = 1;

pert.duration {2} = 1;
pert.type	  {2} = 'ampa_dend';

scatter3(somatapositions(:,1), somatapositions(:,2), somatapositions(:,3), 200, pert.mask{1} + pert.mask{2}*2 -1,'filled'), axis equal


cm = [150 147 130 ; 250 35 29 ; 60 35 250 ; 200 35 190]/255;
colormap(cm)







% apply some current to check the behavior of the cells
I_app = [];




% [===========================================================================================================]
 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
 if ~exist('st_st','var')
	disp('calculating transients')

	 st_st = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, ...
		   	'networksize', [1 1 noneurons] ,'time',steadystate_time ,'W', brick.W ,'ou_noise', ounoise_params , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 st_st.note = 'brick steady state'
end

% 	% st_st.Plist = Plist;
% end


if conjuctive_stimulation
	
	onset_of_inh = 3000;
	onset_of_exc = 3150;
	pert.triggers {1} = onset_of_exc;
	pert.triggers {2} = onset_of_inh;
			
 	sim{1} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', st_st.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', ounoise_params  ,...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

	sim{1}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
	sim{1}.networkHistory.I_cx36 = [];
	sim{1}.networkHistory.backgroundnoise = [];
	sim{1}.note = '41_with_bridges'


		eval(['save curlies_bridges_stim_pair_seeded_' num2str(seed) '_'  date ' -v7.3'])


end


