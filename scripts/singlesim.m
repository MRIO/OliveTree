% noiseamp = .25; tau = 20; sametoall = 0.3; simtype = 'spont'; conntype = 'cluster'     
% noiseamp = .25; tau = 20; sametoall = 0.3; simtype = 'spont'; conntype = 'cluster'     
% noisesig = 1 ; noiseamp = -1; tau = 30; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;    
% gap = eps;  noisesig = .7; noiseamp = 0 ; tau = 20; sametoall = 0.2; simtype = 'spont'; conntype = 'iso' ;  gapcomp = -1.6; 
% gap = 0.04;  noisesig = .7; noiseamp = 0 ; tau = 20; sametoall = 0.2; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = 0.04;  noisesig = .5; noiseamp = -.5 ; tau = 20; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = 0.04;  noisesig = .6; noiseamp = -.6 ; tau = 20; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = eps;  noisesig = -.1; noiseamp = .1 ; tau = 10; sametoall = 0.0; spont = 0; conntype = 'iso' ;  gapcomp = 0;


% gap = 0.04;  noisesig = .4; noiseamp = -.4 ; tau = 30; sametoall = 0.15; spont = 0; conntype = 'iso' ;  gapcomp = 0;
% gap = eps;  noisesig = .4; noiseamp = -.4 ; tau = 30; sametoall = 0.15; spont = 'spont'; conntype = 'iso' ;  gapcomp = 0;
% gap = eps;  noisesig = .3; noiseamp = -.3 ; tau = 30; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 10;

dt = 0.025;

gpu = 1;
to_report = {'V_soma' 'V_dend'};
% [=================================================================]
%  % apply defaults
% [=================================================================]
if ~exist('simtime');			simtime = 1000; end
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
if ~exist('netsize');    		netsize = [5 5 1]; end
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

	rng(seed,'twister')
	if gap > eps
		neurons = createDefaultNeurons(noneurons,'celltypes','randomized');
	else
		neurons = createDefaultNeurons(noneurons,'celltypes','randomized','nogapcompensation',gapcomp);
	end
	% neurons.gbar_ampa_soma = .05*ones(noneurons,1) + .02*rand(noneurons,1);
	% neurons.gbar_ampa_dend = .05*ones(noneurons,1) + .02*rand(noneurons,1);

%============================= perturbation ==============================%



if not(spont) & not(exist('pert'))
	rng(seed,'twister')
	pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
					'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa';
	pert.duration{1}  = 1;
	% pert.triggers{1} = [500 504 508 512 516 850 854 858];
	pert.triggers{1} = [301:4:312];


	pert.mask{2}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
					'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	pert.amplitude{2} = 2;
	pert.type{2}	  = 'ampa_dend';
	pert.duration{2}  = 1;
	pert.triggers{2} = [501:4:512];
elseif not(spont) & exist('pert')
	% use passed pert
	displaytest = ['using pert....'];
else
	pert = [];
end




%          _                 __      __     
%    _____(_)___ ___  __  __/ /___ _/ /____ 
%   / ___/ / __ `__ \/ / / / / __ `/ __/ _ \
%  (__  ) / / / / / / /_/ / / /_/ / /_/  __/
% /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/ 
            

displaytext = ['singlesim'];
rng(seed,'twister')

simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
	'cell_parameters',neurons,'W',W.W*gap ,...
	'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
	'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);





