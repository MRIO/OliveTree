% Input correlations and center peak sharpnesss?

% pspace tests for distributions of ampa input over a wee time

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

cell_function = 'vanilla'; % 'devel'
steadystate_time = 50; %ms
simtime  = 1000; %ms
delta = .025;
gpu = 1;

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

% variables to store
to_report = vsoma;

% nice things


% [================================================]
% 		 create default neurons
% [================================================]
% netsize = [2 3 10];
netsize = [1 10 10];

noneurons = prod(netsize);
def_neurons = createDefaultNeurons(noneurons,'gapcompensation',1);



% [================================================]
% 		gap connections
% [================================================]
 W = zeros(noneurons);
 % W = createW(noneurons);
 % out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)

gaps = [0.05]; symmetrize = 1; nconn = 8; normalizeconn = 1; radius = 3;
W = createW('3d_chebychev', netsize, 3, 1, 1, 0, 0, nconn, [],symmetrize,[0 0 0 0],normalizeconn);

% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
gnoise = [1/5 0.2 0 0];
gnoise = [0 0 0 0];
% [ theta,  sigma  mu seed ]
% th = gaussnoise(1);  sig = gaussnoise(2); mu = gaussnoise(3);
% gnoise = [.2 2 0 0]
 	sametoall = 0;
% [=================================================================]
%  test mask manually
% [=================================================================]
mixturetype = 'manual';
% mixturetype = 'multiple';
switch mixturetype
	case 'manual'

		pert = [];
		nm = 1;
		pert.mask  	  {nm} = rand(noneurons,1)>.5;
		pert.amplitude{nm} = 0;
		pert.triggers {nm} = 100;
		pert.duration {nm} = 800;
		pert.type	  {nm} = 'ou_noise_pos';

		pert.param{nm}(1)  = 1/5 ; % th
		pert.param{nm}(2)  = 0  ;  % mu
		pert.param{nm}(3)  = .2 ;  % sig
		pert.param{nm}(4)  = 1;

		nm = 2;
		pert.mask  	  {nm} = rand(noneurons,1)>.5;
		pert.amplitude{nm} = 0;
		pert.triggers {nm} = 300;
		pert.duration {nm} = 600;
		pert.type	  {nm} = 'ou_noise_neg';

		pert.param{nm}(1)  = 1/20  ;	% th
		pert.param{nm}(2)  = 0  ;	% mu
		pert.param{nm}(3)  = .05 ;	% sig
		pert.param{nm}(4)  = 1;




% [================================================]
%  Ornstein Uhlenbeck Perturbation over masks
% [================================================]
	case 'multiple'
		% create overlapping ou-masks
		pert = [];


		numberofmasks = 10;

		onset_of_stim = 100;
		stim_dur      = 800;
		offset_stim   = 50;
		synapseprobability = .5;

		th =	 1/10 ; % decay time parameter
		sig = 	 .2 ; % pA
		mu = 	 0  ; % pA
		mix =    1 ; % amount of noise shared between neurons in the mask 

			for nm = 1:numberofmasks

				pert.mask  	  {nm} = rand(noneurons,1)>synapseprobability;
				% pert.mask  	  {nm}(nm) = 1;

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

	end


% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 

%%================================================]
% 		 compute transients/steadystate
%=================================================]
 ns = 0;
 for g = gaps
	 	ns = ns +1;


		% if g>0
		% 	neurons = gap_neurons;
		% else
			neurons = def_neurons;
		% end


		   [simresults{ns}] = IOnet_new('cell_parameters', neurons, ...
		   	'networksize', netsize, 'time',simtime ,'W', W.W*g ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);


end
% [===========================================================================================================]

save simresults_ounoisemask_gabaampa 

spks{1} = spikedetect(simresults{1});

