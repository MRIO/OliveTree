% 5 exc 30Hz pulses and 4 inhib pulses centered relative to excitatory pulses but also @30Hz. 
% The excitatory stim began 40ms prior to onset of inh.
% 
% single pulse, EPS following IPS by 50ms

% Questions:
% does it cancel if all cells have the same intrinsic frequency?
% does it cancel in all phases?
% does it cancel for any amplitude of stimulus?

% cell: 658 - 8H

% gap = 0.02;  noisesig = 0; noiseamp = 0 ; tau = 20; sametoall = 0.0; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;
gap = 0.04;  noisesig = 0; noiseamp = 0 ; tau = 20; sametoall = 0.0; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;
spot = 1;

dt = 0.02;
simtime = 2000;
gpu = 0;
% [=================================================================]
%  % create network
% [=================================================================]

if ~exist('conntype');  conntype  = 'iso'; end
if ~exist('spont'); 	spont  = 0; end
if ~exist('saveappliednoise');saveappliednoise = 1; end
if ~exist('sametoall');  sametoall = 0.2; end
if ~exist('rd'); 		 rd = 3        ; end
if ~exist('meannoconn'); meannoconn = 8;end
if ~exist('gap');	     gap = 0.02    ;end
if ~exist('tau'); 		 tau = 20      ;end
if ~exist('noiseamp'); 	 noiseamp = .5 ;end
if ~exist('noisesig'); 	 noisesig = 0  ;end

netsize = [1 3 3];
% netsize = [1 20 1];
	noneurons = prod(netsize);

plotthis  = 0;

normleak  = 1;
randomize = 1;
scaling   = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;


activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
activations =  {'V_soma','V_dend' ,'Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Potassium_n', 'Potassium_x_s', 'Hcurrent_q', 'Hcurrent_q','g_gaba_dend', 'g_ampa_soma'};
currents = {'V_soma','V_dend','I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s'};
soma = {'V_soma'}
% to_report = activations;
to_report = soma;



switch conntype
	case 'iso'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
	case 'cluster'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 25 .5 .05], normleak);
end

% W.W = [0 0 0 1 1 0 1 0 0;
% 	   0 0 0 0 0 0 0 0 0;
% 	   0 0 0 0 1 1 0 0 1;
% 	   1 0 0 0 1 0 1 0 0;
% 	   1 0 1 1 0 1 1 0 1;
% 	   0 0 1 0 1 0 0 0 1;
% 	   1 0 0 1 1 0 0 0 0;
% 	   0 0 0 0 0 0 0 0 0;
% 	   0 0 1 0 1 1 0 0 0]*gap/3;

% [=================================================================]
%  % create neurons
% [=================================================================]

rng(0,'twister')
neurons = createDefaultNeurons(noneurons, 'celltypes' , 'randomized', 'gapcompensation',gapcomp);

neurons.gbar_ampa_soma = .1*ones(noneurons,1);
neurons.gbar_gaba_soma = 0*ones(noneurons,1);

neurons.gbar_ampa_soma = .1*ones(noneurons,1);
neurons.gbar_gaba_dend = .8*ones(noneurons,1);
% neurons.g_CaL = ones(noneurons,1)*1.1;
% neurons.g_CaL(5) = 1.1;

neurons.g_CaL = linspace(.9,1.3,noneurons);
% neurons.g_h = linspace(.1,3,noneurons);

%============================= perturbation ==============================%



noise_level = [1/tau noisesig noiseamp 0];

if not(spont)
	npulsesI  = 4;
	npulsesE  = 5;
	bfreq    = 30; bT = round(1000/30); 
	isdelay  = 40;
	onset 	 = 630;
	pulsesI	 = round(linspace(0,(npulsesI-1)*bT,npulsesI));
	pulsesE	 = round(linspace(0,(npulsesE-1)*bT,npulsesE));
	


	% train = cumsum(ones(1,npulses)*bT)
	% pert.mask{1}  	  = [0 0 0 0 1 0 0 0 0]';
	% pert.mask{1}  	  = [0 0 0 0 1 1 1 1 1]';
	pert.mask{1}  	  = [1 1 1 1 1 1 1 1 1]';
	
	% pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa';
	pert.duration{1}  = 1;
	% pert.triggers{1} = [500 504 508 512 516 850 854 858];
	% pert.triggers{1} = [300:30:330 400:5:430];
	pert.triggers{1} = onset + isdelay + pulsesE;

	% pert.mask{2}	= [0 0 0 0 1 0 0 0 0]';
	% pert.mask{2}  	  = [0 0 0 0 1 1 1 1 1]';
	pert.mask{2}  	  = [1 1 1 1 1 1 1 1 1]';
	% pert.mask{2}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{2} = 2;
	pert.type{2}	  = 'gaba_dend';
	% pert.type{2}	  = 'gaba_soma';
	pert.duration{2}  = 1;
	pert.triggers{2} = onset + pulsesI;
else
	pert = [];
end




%          _                 __      __     
%    _____(_)___ ___  __  __/ /___ _/ /____ 
%   / ___/ / __ `__ \/ / / / / __ `/ __/ _ \
%  (__  ) / / / / / / /_/ / / /_/ / /_/  __/
% /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/ 
                                          


displaytext = ['singlesim 4 tycho'];


simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
	'cell_parameters',neurons,'W',W.W*gap ,...
	'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
	'saveappliednoise',saveappliednoise, 'displaytext',displaytext,'to_report', to_report);





