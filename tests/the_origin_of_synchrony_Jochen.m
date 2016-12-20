% Input correlations and center peak sharpnesss?

% pspace tests for distributions of ampa input over a wee time

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

cell_function = 'vanilla'; % 'devel'
steadystate_time = 1000; %ms
simtime  = 500; %ms
delta = .02;
gpu = 1;

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
selection = {'V_soma' 'I_cx36'};

% variables to store
to_report = selection ;


% [================================================]
% 		 create default neurons
% [================================================]
noneurons = 30;
netsize = [3 10 1];



rng(0,'twister');
gap_neurons = createDefaultNeurons(noneurons,'celltypes', 'randomized','gapcompensation', 1);
gap_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; %linspace(0,.3,noneurons);
gap_neurons.gbar_gaba_soma = ones(noneurons,1)*.5; %ones(noneurons,1)*.3;

rng(0,'twister')
def_neurons = createDefaultNeurons(noneurons,'celltypes', 'randomized','gapcompensation', 0);
def_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; %linspace(0,.3,noneurons);
def_neurons.gbar_gaba_soma = ones(noneurons,1)*.5; %ones(noneurons,1)*.3;


% [================================================]
% 		gap connections
% [================================================]
 W = zeros(noneurons);
 % W = createW(noneurons);

gaps = [0.0 0.025];
noconnections = 10;
symmetrize = 1;
randomize = 0;
normalize = 1;
plotconns = 1;
scaling = 1;

%topological
W1 = createW('3d_chebychev', netsize, 5,    scaling, randomize, plotconns, 0, noconnections, [],symmetrize, [0 0 0 0], normalize);
mean(sum(W1.W>0))
%all to all
W2 = createW('all to all'  , netsize, [],     scaling, randomize, plotconns, 0, noconnections, [],symmetrize, [0 0 0 0], normalize);
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)
mean(sum(W2.W>0))

% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% gnoise = [theta sigma mu seed];
gnoise = [1/5 0.1 0 0]; sametoall = 1;
% [================================================]
%  Ornstein Uhlenbeck Perturbation over masks
% [================================================]

% create overlapping ou-masks

numberofmasks = noneurons;

onset_of_stim = 100;
stim_dur      = simtime;
offset_stim   = 100;
% synapseprobability = .3;

th =	 1/5 ; % decay time parameter
mu = 	 0 ; % pA
sig = 	 .9; % pA
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


%============================= topological ==============================%
 clear simresults
 ns = 0;
 for g = 0.025
	 	ns = ns +1;

		neurons = gap_neurons;
		
		[simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W1.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


			if sum(simresults{ns}.failed)

				[simresults{ns}] = IOnet('tempState', simresults{part-1}.lastState ,'cell_parameters', neurons, ...
			   	'networksize', netsize, 'time',simtime ,'W', W1.W*g ,'ou_noise', gnoise , ...
			   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
			   	'cell_function', cell_function ,'delta',0.01,'sametoall', sametoall);
				
			end

		   

end

replayResults_3(simresults{1})
simresults{1}.topological = 1;

save topological simresults

%=============================flush==============================%

gpuDevice([])
gpuDevice(1)

%============================= all to all ==============================%

 for g = 0.025

		neurons = gap_neurons;
		
		[simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


			if sum(simresults{ns}.failed)

				[simresults{ns}] = IOnet('tempState', simresults{part-1}.lastState ,'cell_parameters', neurons, ...
			   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
			   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
			   	'cell_function', cell_function ,'delta',0.01,'sametoall', sametoall);
				
			end

		   

end

save alltoall simresults

%============================= no gaps ==============================%

clear simresults
 for g = 0.001
	 	ns = ns +1;

		neurons = gap_neurons;
		
		[simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


			if sum(simresults{ns}.failed)

				[simresults{ns}] = IOnet('tempState', simresults{part-1}.lastState ,'cell_parameters', neurons, ...
			   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
			   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
			   	'cell_function', cell_function ,'delta',0.01,'sametoall', sametoall);
				
			end

		   

end

replayResults_3(simresults{1})
save nogaps simresults


 for g = 0.025

		neurons = gap_neurons;
		
		[simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


			if sum(simresults{ns}.failed)

				[simresults{ns}] = IOnet('tempState', simresults{part-1}.lastState ,'cell_parameters', neurons, ...
			   	'networksize', netsize, 'time',simtime ,'W', W2.W*g ,'ou_noise', gnoise , ...
			   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
			   	'cell_function', cell_function ,'delta',0.01,'sametoall', sametoall);
				
			end

		   

end

replayResults_3(simresults{1})
save nonoise_alltoall simresults
