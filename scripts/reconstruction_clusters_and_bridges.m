% demo_clusters.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear
rng(0,'twister') % random seed

steadystate_time = 500; %ms
simtime  = 3000; %ms
delta = .02;
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
% netsize = [3 15 15];

% def_neurons = createDefaultNeurons(noneurons,'celltypes','param_sweep');
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized2');
% neurons.g_CaL = linspace(.5, 1, noneurons);


% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)

nconns_curlies = 5;
nconns_bridges = 5;
gap_curlies = .05;
gap_bridges = .05;
plotconn = 1;
normalize = 1;


load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);

if not(exist('curlies'))
	curlies = createW('3d_reconstruction', [], 5*40, 1, 0, 1, [], nconns_curlies, somatapositions,1,[1 20 1 0]);
	bridges = createW('3d_reconstruction', [], 10*40, 1, 0, 1, [], nconns_bridges, somatapositions,1,[1 20 0 1]);

	% 10% of cells are bridges
	bc = randperm(noneurons); 
	z = zeros(noneurons,1) ; 
	z(bc(1:round(.1*noneurons))) = 1;
	bc =z;

	curlies.W = bsxfun(@times, curlies.W, ~bc);
	curlies.W = bsxfun(@times, curlies.W, ~(bc'))*gap_curlies;

	bridges.W = bsxfun(@times, bridges.W, z);
	bridges.W = (bridges.W+bridges.W')*gap_bridges;

	bridg_curlies.coords = curlies.coords;
	bridg_curlies.stats = connectivity_statistics(bridg_curlies);
	bridg_curlies.stats.clusters = curlies.stats.clusters;
	bridg_curlies.W = curlies.W + bridges.W;

	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridg_curlies.stats.clusters)
end




% Wcluster150 = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], plotconn, [1 150 .9 .01],1);


% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% gnoise = [.2 .3 0 5];
gnoise = [0 0 0 0];
sametoall = 0.05


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]

% create overlapping ampa masks

% numberofmasks = 10; 
onset_of_stim = [505:5:525];

% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

% pert.mask     {1} =  create_input_mask(netsize, 'dist_to_center','radius',2, 'synapseprobability', 1,'plotme',1);
pert.mask     {1} =  curlies.stats.clusters==5;
pert.amplitude{1} = 1;
pert.triggers {1} = onset_of_stim;
pert.duration {1} = 5;
pert.type	  {1} = 'gaba_soma';
% pert.type	  {1} = 'ampa';





% [===========================================================================================================]
 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
% if ~exist('st_st','var')
% 	disp('calculating transients')
% 	st_st = IOnet_new('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
% 	% st_st.Plist = Plist;
% end


% [=================================================================]
%  GABA
% [=================================================================]

% 'tempState', st_st.lastState,
	 sim{1} = IOnet_new( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', bridg_curlies.W ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 sim{1}.note = 'curlies and bridges'


% [=================================================================]
%  spontaneous
% [=================================================================]



save bridges_curlies_wspace

generatemovies = 0;
if generatemovies
	for cc = 1:6
		animate_volume(simR{cc},[1:3000],1,1)
		eval(['!mv volume.mp4 ' num2str(cc) '_3p_nonoise_10con_norm_g05.mp4'])
	end
end













