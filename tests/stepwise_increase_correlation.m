% stepwise_increase_correlation.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

steadystate_time = 500; %ms
simtime  = 10000; %ms
delta = .025;
gpu = 1;

if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
thisseed.Seed



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

nconns_curlies   = 8;
nconns_bridges   = 8;
cells_in_cluster = 20;
gap_curlies      = .05;
gap_bridges      = .05;
plotconn 		 = 1;
normalize 		 = 1;

load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);

brick = createW('3d_reconstruction', [], 8*40, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
brick_bu = brick;
brick.W = brick.W*gap_curlies;



% [=================================================================]
%  create cells
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized', 'rng', thisseed) 



% [================================================]
%  Ornstein Uhlenbeck noise over masks
% [================================================]

% create non-overlapping ou-masks


onset_of_stim = 100;
stim_dur      = 9000;
offset_stim   = 100;
synapseprobability = 1;

th =	 1/5 ; % decay time parameter
mu = 	 -.6 ; % pA
sig = 	 .6 ; % pA
mix =    0;



pert.mask     {1} = create_input_mask(somatapositions, 'reconstruction','synapseprobability',synapseprobability, 'radius', 100, 'plotme',1)
pert.mask     {2} = double(not(pert.mask{1}))

pert.amplitude{1} = 1;
pert.triggers {1} = 0;
pert.duration {1} = stim_dur;
pert.type	  {1} = 'ou_noise';

pert.amplitude{2} = 1;
pert.triggers {2} = 0
pert.duration {2} = stim_dur;
pert.type	  {2} = 'ou_noise';

pert.param{1}  = [1/5 sig mu mix double(thisseed.Seed)];  % no correlation outside the mask
pert.param{2}  = [1/5 sig mu mix double(thisseed.Seed)]; % 50% correlation in the reveiving mask


%%================================================]
% 		 compute transients/steadystate
%=================================================]

sim{1} = IOnet( 'cell_parameters', def_neurons, ...
 		'perturbation', pert, ... 
	   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W , ... 
	   	'to_report', to_report ,'gpu', gpu , ...
	   	'cell_function', cell_function ,'delta',delta);


mix_increments = [0:.1:.5];
for ii = 2:6;

	pert.param{1}  = [1/5 sig mu mix_increments(ii) double(thisseed.Seed)];  % no correlation outside the mask

 	sim{ii} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', sim{ii-1}.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta);
	sim{ii}.note = ['incremental noise' num2str(ii)];
	sim{ii}.W = brick.W;
	sim{ii}.networkHistory.V_soma = single(sim{1}.networkHistory.V_soma);
	sim{ii}.networkHistory.I_cx36 = [];
	sim{ii}.networkHistory.backgroundnoise = [];
end



