% stepwise_increase_correlation.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

simtime  = 5000; %ms
delta = .01;
gpu = 1;

stim_dur      = 1500;
mix_increments = [.1 .2 .3 .4 .5 ];


if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister'); % random seed only for simulations (not for network cells)


% [================================================]
% variables to report
% [================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};
noise = {'V_soma' 'backgroundnoise'}

% variables to store
to_report = noise;



% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)

gap      = .1;
nconns   = 5;
cells_in_cluster = 20;

plotconn 		 = 1;
normalize 		 = 1;

load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);

brick = createW('3d_reconstruction', [], 4*40, 1, 0, plotconn, [], nconns, somatapositions,1,[0 0 0 0]);
brick_bu = brick;
brick.W = brick.W*gap;
W = brick;

% curlies = createW('3d_reconstruction', [], 4*40, 1, 0, plotconn, [], nconns, somatapositions,1,[1 cells_in_cluster 1 0]);

plotnetstruct(brick.W, somatapositions(:,1), somatapositions(:,2), somatapositions(:,3), ones(length(somatapositions),1))


% [=================================================================]
%  create cells
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
% def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized2', 'rng', thisseed) 
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized', 'rng', thisseed) 



% [================================================]
%  Ornstein Uhlenbeck noise over masks
% [================================================]

% create non-overlapping ou-masks


synapseprobability = 1;

th =	 1/20 ; % decay time parameter
mu = 	 -.8 ; % pA
sig = 	  .8 ; % pA
radius = 150;

pert.mask     {1} = create_input_mask(somatapositions, 'reconstruction','synapseprobability',synapseprobability, 'radius', 100, 'plotme',1)
pert.mask     {2} = not(pert.mask{1});

pert.amplitude{1} = 1;
pert.triggers {1} = 50;
pert.duration {1} = stim_dur;
pert.type	  {1} = 'ou_noise';

pert.amplitude{2} = 1;
pert.triggers {2} = 50;
pert.duration {2} = stim_dur;
pert.type	  {2} = 'ou_noise';

pert.param{1}  = [th sig mu mix_increments(1) ];  % no correlation outside the mask
pert.param{2}  = [th sig mu 0 ]; % 50% correlation in the reveiving mask


%%================================================]
% 		 compute transients/steadystate
%=================================================]

sim{1} = IOnet( 'cell_parameters', def_neurons, ...
 		'perturbation', pert, ... 
	   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', W , ... 
	   	'to_report', to_report ,'gpu', gpu , ...
	   	'cell_function', cell_function ,'delta',delta);


for ii = 2:length(mix_increments(2:end))+1

	pert.param{1}  = [th sig mu mix_increments(ii)];  % no correlation outside the mask

 	sim{ii} = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, 'tempState', sim{ii-1}.lastState, ...
		   	'networksize', [1 1 noneurons] ,'time',simtime ,'W', W , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta);
	sim{ii}.note = ['incremental noise' num2str(ii)];
	sim{ii}.W = brick.W;
	sim{ii}.networkHistory.V_soma = single(sim{ii}.networkHistory.V_soma);
	
end



