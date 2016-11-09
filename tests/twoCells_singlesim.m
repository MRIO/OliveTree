
    % CAL     IH     g_int     g_l     p1     g_K_Ca    arb    freq_each     ampl     meanVm     spks    supth
    % ___    ____    _____    _____    ___    ______    ___    _________    ______    _______    ____    _____

    % 8      0.12     0.3      0.01    0.1    40        0.5    2.0075       103.07    -65.567    3       0.177
    % 4      0.12    0.05     0.015    0.1    45          1    2.5624        10.48    -61.681    0       0.911
    % 4       1.2    0.05      0.01    0.2    45          1    2.7262       8.3672     -58.12    0       0.683
    % 6       1.2     0.3      0.01    0.1    45          1    2.8419       8.5423    -65.029    0       0.484
    % 8      0.12    0.15      0.01    0.3    45        0.5    3.0229       104.47    -62.334    4       0.181
    % 8      0.12     0.3      0.01    0.1    50          1    3.5372       12.925    -67.466    0       0.891
    % 6      0.12     0.3     0.015    0.1    45        0.5    3.7255       79.369    -64.644    3       0.207
    % 8      0.12     0.3      0.01    0.3    45        0.5    3.8488        82.62    -65.647    3       0.145
    % 8       1.2    0.15     0.015    0.1    50          1    3.9049       10.304    -62.273    0       0.463
    % 6       1.2    0.05     0.015    0.3    40        0.5    4.6433        94.56    -30.468    5       0.868
    % 6      0.12    0.15      0.01    0.3    45        0.5    4.8059       102.36    -62.013    3       0.208
    % 6       1.2     0.3      0.01    0.2    50          1    5.0459       78.676    -63.626    5       0.237
    % 4      0.12    0.15      0.01    0.2    45        0.5    5.7248       76.379     -64.26    4       0.211
    % 4       1.2     0.3      0.01    0.3    40        0.5    5.8098        74.44    -63.933    5       0.206
    % 8       1.2     0.3      0.01    0.1    40          1    5.8912       8.7937    -63.581    0       0.297
    % 8       1.2     0.3     0.015    0.2    40          1    6.4458       99.061     -62.43    6       0.194
    % 8       1.2    0.15      0.01    0.1    45        0.5    6.7243        84.21    -59.505    6       0.237
    % 8       1.2    0.15      0.01    0.2    50        0.5    6.8069        94.75    -58.322    6       0.233
    % 4       1.2    0.15     0.015    0.3    50        0.5    7.0728       85.807    -55.067    8       0.291
    % 6       1.2     0.3     0.015    0.1    50        0.5    7.5533       76.058     -61.81    5       0.284
    % 8      0.12    0.05      0.01    0.2    50        0.5     7.748       72.656    -54.415    6       0.267
    % 8       1.2    0.15      0.01    0.3    40        0.5    7.8282       112.39    -55.957    7       0.244
    % 6       1.2    0.05      0.01    0.1    45        0.5    8.2362       24.384    -55.074    0       0.441
    % 8      0.12    0.05     0.015    0.2    50          1    8.7284       16.921    -55.244    0        0.49
    % 8       1.2    0.05     0.015    0.1    50        0.5     9.122       52.024    -49.349    9       0.364
    % 6       1.2    0.15     0.015    0.2    50          1    9.7515       10.624    -60.544    0       0.454


    % 8      0.12     0.3      0.01    0.1    40        0.5    2.0075       103.07    -65.567    3       0.177
    % 8      0.12     0.3      0.01    0.3    45        0.5    3.8488        82.62    -65.647    3       0.145

    % 6       1.2    0.15     0.015    0.2    50          1    9.7515       10.624    -60.544    0       0.454

displaytext = ['singlesim 2 cells to check for correlation lag as a function of Ca'];

gap = 0.2;  noisesig = .6; noiseamp = -.3 ; tau = 20; sametoall = 0.5; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;

dt = 0.02;
simtime = 5000;
gpu = 0;
% [=================================================================]
%  % create network
% [=================================================================]

if ~exist('conntype');  conntype  = 'iso'; end
if ~exist('spont'); 	spont  = 1; end
if ~exist('saveappliednoise');saveappliednoise = 1; end
if ~exist('sametoall');  sametoall = 0.2; end
if ~exist('rd'); 		 rd = 3; end
if ~exist('meannoconn'); meannoconn = 8;end
if ~exist('gap');	     gap = 0.02 ;end
if ~exist('tau'); 		 tau = 20 ;end
if ~exist('noiseamp'); 	 noiseamp = .5 ;end
if ~exist('noisesig'); 	 noisesig = 0 ;end


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
activations =  {'V_soma','Calcium_l'};

to_report = activations;


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

load dev_cells_pspace
[freqsorted oo] = sortrows(resultstable.allneurons,8);


select_cells = [8 24];

netsize = [1 length(select_cells) 1];
noneurons = prod(netsize);

neurons = createDefaultNeurons(length(select_cells));


c = 0;
for Ncell = select_cells
	c= c+1;
	neurons.g_CaL    (c,1) = table2array(freqsorted(Ncell,1));
	neurons.g_h      (c,1) = table2array(freqsorted(Ncell,2));
	neurons.g_int    (c,1) = table2array(freqsorted(Ncell,3));
	neurons.g_l      (c,1) = table2array(freqsorted(Ncell,4));
	neurons.p1       (c,1) = table2array(freqsorted(Ncell,5));
	neurons.g_K_Ca   (c,1) = table2array(freqsorted(Ncell,6));
	neurons.arbitrary(c,1) = table2array(freqsorted(Ncell,7));
end




neurons.gbar_ampa_soma = .1*ones(noneurons,1);
neurons.gbar_gaba_soma = .25*ones(noneurons,1);

%============================= perturbation ==============================%



noise_level = [1/tau noisesig noiseamp 0];

if not(spont)
	npulses  = 0;
	bfreq    = 30; bT = round(1000/30);
	isdelay  = 50;
	onset 	 = 650;


	% train = cumsum(ones(1,npulses)*bT)
	% pert.mask{1}  	  = [0 0 0 0 1 0 0 0 0]';
	pert.mask{1}  	  = [0 0 0 1 1 1 1 1 1]';
	
	% pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa';
	pert.duration{1}  = 1;
	% pert.triggers{1} = [500 504 508 512 516 850 854 858];
	% pert.triggers{1} = [300:30:330 400:5:430];
	pert.triggers{1} = onset + isdelay;

	% pert.mask{2}	= [0 0 0 0 1 0 0 0 0]';
	pert.mask{2}  	  = [0 0 0 1 1 1 1 1 1]';
	% pert.mask{2}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{2} = 2;
	pert.type{2}	  = 'gaba_dend';
	pert.duration{2}  = 1;
	pert.triggers{2} = onset ;
else
	pert = [];
end




%          _                 __      __     
%    _____(_)___ ___  __  __/ /___ _/ /____ 
%   / ___/ / __ `__ \/ / / / / __ `/ __/ _ \
%  (__  ) / / / / / / /_/ / / /_/ / /_/  __/
% /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/ 
                                          




simresults = IOnet_new('networksize', netsize,'time',simtime,'delta',dt,...
	'cell_parameters',neurons,'W',W.W*gap ,...
	'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
	'saveappliednoise',saveappliednoise, 'displaytext',displaytext,'to_report', to_report);





