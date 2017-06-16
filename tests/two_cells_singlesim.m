gap = 0.0;  noisesig = 0; noiseamp = 0 ; tau = 20; sametoall = 0.2; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;

dt = 0.02;
simtime = 500;
gpu = 0;
% [=================================================================]
%  % create network
% [=================================================================]

if ~exist('conntype');  conntype  = 'iso'; end
if ~exist('spont'); 	spont  = 0; end
if ~exist('saveappliednoise');saveappliednoise = 1; end
if ~exist('sametoall');  sametoall = 0.2; end
if ~exist('rd'); 		 rd = 3; end
if ~exist('meannoconn'); meannoconn = 8;end
if ~exist('gap');	     gap = 0.02 ;end
if ~exist('tau'); 		 tau = 20 ;end
if ~exist('noiseamp'); 	 noiseamp = .5 ;end
if ~exist('noisesig'); 	 noisesig = 0 ;end

netsize = [1 2 1];
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
currentstoo =  {'V_soma','Calcium_l', 'Calcium_r', 'I_CaH' ,  'Ca2Plus'};

to_report = currentstoo;


switch conntype
	case 'iso'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
	case 'cluster'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 25 .5 .05], normleak);
end



% [=================================================================]
%  % create neurons
% [=================================================================]

rng(0,'twister')
neurons = createDefaultNeurons(noneurons, 'gapcompensation',gapcomp);
neurons.gbar_ampa_soma = .3*ones(noneurons,1);
neurons.gbar_gaba_soma = .25*ones(noneurons,1);
neurons.g_CaL = ones(noneurons,1)*.7;
neurons.g_CaL(1) = 2.;
neurons.g_CaL(2) = .5;
neurons.g_CaH(2) = .5;

%============================= perturbation ==============================%



noise_level = [1/tau noisesig noiseamp 0];

if not(spont)
	npulses  = 5;
	bfreq    = 3; bT = round(1000/bfreq);
	isdelay  = 0;
	isramp 	 = linspace(0,npulses*bfreq, npulses)
	onset 	 = 450;


	pert.mask{1}  	  = [1 1]';
	pert.type{1}	  = 'gaba_soma';
	pert.duration{1}  = 1;
	pert.triggers{1} = onset + isdelay + isramp	;
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



replayResults_3(simresults,'plotallfields',1)

