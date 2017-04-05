% noiseamp = .25; tau = 20; sametoall = 0.3; simtype = 'spont'; conntype = 'cluster'     
% noiseamp = .25; tau = 20; sametoall = 0.3; simtype = 'spont'; conntype = 'cluster'     
% noisesig = 1 ; noiseamp = -1; tau = 30; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;    
% gap = eps;  noisesig = .7; noiseamp = 0 ; tau = 20; sametoall = 0.2; simtype = 'spont'; conntype = 'iso' ;  gapcomp = -1.6; 
% gap = 0.04;  noisesig = .7; noiseamp = 0 ; tau = 20; sametoall = 0.2; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = 0.04;  noisesig = .5; noiseamp = -.5 ; tau = 20; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = 0.04;  noisesig = .6; noiseamp = -.6 ; tau = 20; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = eps;  noisesig = -.1; noiseamp = .1 ; tau = 10; sametoall = 0.0; spont = 0; conntype = 'iso' ;  gapcomp = 0;

dt = 0.025;
simtime = 5000;
gpu = 1;
% [=================================================================]
%  % apply defaults
% [=================================================================]

if ~exist('conntype');	 		conntype  = 'iso'; end
if ~exist('spont'); 	 		spont  = 0; end
if ~exist('saveappliednoise');  saveappliednoise = 1; end
if ~exist('sametoall');  		sametoall = 0.2; end
if ~exist('rd'); 		 		rd = 3; end
if ~exist('meannoconn'); 		meannoconn = 8;end
if ~exist('gap');	     		gap = 0.04 ;end
if ~exist('tau'); 		 		tau = 20 ;end
if ~exist('noiseamp'); 	 		noiseamp = .5 ;end
if ~exist('noisesig'); 	 		noisesig = 0 ;end
if ~exist('netsize');    		netsize = [2 1 1]; end
if ~exist('seed');  	 		seed = 0; end


% [=================================================================]
%  create network
% [=================================================================]

	noneurons = prod(netsize);
	noise_level = [1/tau noisesig noiseamp 0];

plotthis  = 0;

normleak  = 1;
randomize = 1;
scaling   = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;



switch conntype
	case 'iso'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
	case 'cluster'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 25 .5 .05], normleak);
end



% [=================================================================]
%  % create neurons
% [=================================================================]



if not(isempty(neurons))
	rng(seed,'twister')
	neurons = createDefaultNeurons(noneurons,'celltypes','randomized','gapcompensation',gapcomp);
	neurons.gbar_ampa_soma = .05*ones(noneurons,1) + .02*rand(noneurons,1);
end

%============================= perturbation ==============================%



if not(spont)
	rng(seed,'twister')
	pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
					'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa';
	pert.duration{1}  = 1;
	% pert.triggers{1} = [500 504 508 512 516 850 854 858];
	pert.triggers{1} = [1500];


	% pert.mask{2}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{2} = 2;
	% pert.type{2}	  = 'ampa';
	% pert.duration{2}  = 1;
	% pert.triggers{2} = [500 504 508 512 516 850 854 858];
else
	pert = [];
end




%          _                 __      __     
%    _____(_)___ ___  __  __/ /___ _/ /____ 
%   / ___/ / __ `__ \/ / / / / __ `/ __/ _ \
%  (__  ) / / / / / / /_/ / / /_/ / /_/  __/
% /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/ 
            
rng(0,'twister')                              


displaytext = ['singlesim'];
rng(seed,'twister')

simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
	'cell_parameters',neurons,'W',W.W*gap ,...
	'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
	'saveappliednoise',saveappliednoise, 'displaytext',displaytext);





