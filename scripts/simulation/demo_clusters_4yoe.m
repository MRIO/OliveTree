% demo_clusters.m

% [================================================]
% 		 simulation parameters
% [================================================]
% clear
rng(0,'twister') % random seed

steadystate_time = 500; %ms
simtime  = 2000; %ms
delta = .02;
gpu = 1;

computebrick		= 0;
computecluster50	= 1;
computeonecluster	= 0;


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
netsize = [3 10 10];
noneurons = prod(netsize);

% def_neurons = createDefaultNeurons(noneurons,'celltypes','param_sweep');
def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized');
% neurons.g_CaL = linspace(.5, 1, noneurons);


% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)

nconns = 10;
gap = .05;
plotconn = 1;
normalize = 1;

Wbrick 		= createW('3d_chebychev', netsize, 3, 1, 1, 1, [], nconns, [], plotconn, [0 0 0 0], normalize);
Wcluster50  = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], nconns, [], plotconn, [1 50 .5 0], normalize);
Wonecluster 		= createW('one_cluster', netsize, 3, 1, 1, 1, [], nconns, [], plotconn, [0 0 1 .1], normalize);

% Wcluster150 = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], plotconn, [1 150 .9 .01],1);


% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% gnoise = [.2 .3 0 5];
gnoise = [1/30 .5 -.5 0];
sametoall = 0.1;


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]

% create overlapping ampa masks

% numberofmasks = 10; 
stim_interval = [25 50 100 250 500];
onset_of_stim = [505:5:525];

% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

pert.mask     {1} =  create_input_mask(netsize, 'dist_to_center','radius',3, 'synapseprobability', 1,'plotme',1);
pert.amplitude{1} = 1;
pert.triggers {1} = onset_of_stim;
pert.duration {1} = 5;
pert.type	  {1} = 'gaba_soma';
pert.type	  {1} = 'ampa';
			% 



% [===========================================================================================================]
 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
	% st_st.Plist = Plist;
end


% [=================================================================]
%  GABA
% [=================================================================]


ccc = 1;
% W brick

if computebrick
	 simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
	 		'perturbation', pert, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wbrick.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	simR{ccc}.Wbrick = Wbrick;

end


if computecluster50
	ccc = ccc+1;

	simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
			'perturbation', pert, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wcluster50.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	simR{ccc}.Wcluster50 = Wcluster50;

	% W3
end

if computeonecluster
	ccc = ccc+1;

	simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
			'perturbation', pert, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wonecluster.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	simR{ccc}.Wonecluster = Wonecluster;
end


% [=================================================================]
%  spontaneous
% [=================================================================]

if exist('ccc','var')
	ccc = ccc+1;
else
	ccc = 2; V = []; st_st.lastState = st_st.lastState;
end


if computebrick
	 simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wbrick.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

	
end

if computecluster50
	ccc = ccc+1;
	simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wcluster50.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
			simR{ccc}.Wcluster50 = Wcluster50;
end

if computeonecluster

	ccc = ccc+1;
	% W3


	simR{ccc} = IOnet('tempState', st_st.lastState, 'cell_parameters', def_neurons, ...
		   	'networksize', netsize ,'time',simtime ,'W', Wonecluster.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
end



generatemovies = 0;
if generatemovies
	for cc = 1:6
		animate_volume(simR{cc},[1:3000],1,1)
		eval(['!mv volume.mp4 ' num2str(cc) '_3p_nonoise_10con_norm_g05.mp4'])
	end
end













