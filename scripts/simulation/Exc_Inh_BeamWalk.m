% exciatation walking over inhibition
% According to Tycho and Sebastian
% single pulse, EPS following IPS by 50ms

% Questions:
% does it cancel if all cells have the same intrinsic frequency?
% does it cancel in all phases?
% does it cancel for any amplitude of stimulus?

% cell: 658 - 8H

% gap = 0.02;  noisesig = 0; noiseamp = 0 ; tau = 20; sametoall = 0.0; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;
gap = 0.0;  noisesig = 0; noiseamp = 0 ; tau = 20; sametoall = 0.0; simtype = 'burst'; conntype = 'iso' ;  gapcomp = 0;
spot = 1;
interstimT = 1000;

dt = 0.02;
gpu = 0;
% [=================================================================]
%  % create network
% [=================================================================]

if ~exist('conntype');  conntype  = 'iso'; end
if ~exist('spont'); 	spont  = 0; end
if ~exist('saveappliednoise');saveappliednoise = 1; end
if ~exist('sametoall');  sametoall = 0; end
if ~exist('rd'); 		 rd = 3        ; end
if ~exist('meannoconn'); meannoconn = 4;end
if ~exist('gap');	     gap = 0.02    ;end
if ~exist('tau'); 		 tau = 20      ;end
if ~exist('noiseamp'); 	 noiseamp = 0 ;end
if ~exist('noisesig'); 	 noisesig = 0.1  ;end

noneurons = 22;

netsize = [1 noneurons 1];
% netsize = [1 20 1];
	noneurons = prod(netsize);

cell_function = 'vanilla';

w = 298;

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
currents = {'V_soma','V_dend','I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_h', 'I_h_s'};
soma = {'V_soma'};
% to_report = activations;
to_report = soma;


if noneurons>1
	switch conntype
		case 'iso'
			W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
		case 'cluster'
			W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 25 .5 .05], normleak);
	end
else
	W.W = 0;
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
cellset = 'cellset_vanilla'
neurons = createDefaultNeurons(noneurons, 'celltypes' , cellset , 'addrand', 1);

neurons.gbar_ampa_dend = .1*ones(noneurons,1);
neurons.gbar_ampa_soma = .2*ones(noneurons,1);
neurons.gbar_gaba_dend = 1*ones(noneurons,1);

% neurons.gbar_gaba_soma = 0*ones(noneurons,1);
% neurons.g_CaL = ones(noneurons,1)*1.1;
% neurons.g_CaL(5) = 1.1;


% neurons.g_h = linspace(.1,3,noneurons);

%============================= perturbation ==============================%



% noise_level = [1/tau noisesig noiseamp 0];
noise_level = [0 0 0 0];


if not(spont)
	
	% inhibitory - 20ms (gaba)
	% excitatory - 

	pulsesI	 = [1000:interstimT:50000];
	% pulsesI	 = [500:700:2000];
	pulsesI = pulsesI + round(rand(size(pulsesI))*50);
	pulsesE	 = pulsesI + round(linspace(-100, 100  , length(pulsesI)));
	
	% train = cumsum(ones(1,npulses)*bT)
	% pert.mask{1}  	  = [0 0 0 0 1 0 0 0 0]';
	% pert.mask{1}  	  = [0 0 0 0 1 1 1 1 1]';
	pert.mask{1}  	  = ones(noneurons,1);
	
	% pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa_dend';
	pert.duration{1}  = 1;
	% pert.triggers{1} = [500 504 508 512 516 850 854 858];
	% pert.triggers{1} = [300:30:330 400:5:430];
	pert.triggers{1} = pulsesE;

	% pert.mask{2}	= [0 0 0 0 1 0 0 0 0]';
	% pert.mask{2}  	  = [0 0 0 0 1 1 1 1 1]';
	pert.mask{2}  	  = ones(noneurons,1);
	% pert.mask{2}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', 2, ...
	% 				'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',0);
	% pert.amplitude{2} = 2;
	% pert.type{2}	  = 'gaba_dend';
	pert.type{2}	  = 'gaba_soma';
	pert.duration{2}  = 1;
	pert.triggers{2} = pulsesI;
else
	pert = [];
end

simtime = max([pulsesI pulsesE])+w;
% simtime = 320;

%          _                 __      __     
%    _____(_)___ ___  __  __/ /___ _/ /____ 
%   / ___/ / __ `__ \/ / / / / __ `/ __/ _ \
%  (__  ) / / / / / / /_/ / / /_/ / /_/  __/
% /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/ 
                                          

if ~exist('simresults')
	displaytext = ['exc and inh over each other'];


	simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,'gpu', gpu, ...
		'cell_function', cell_function, 'cell_parameters',neurons,'W',W.W*gap , ...
		'ou_noise', noise_level , 'perturbation', pert, 'sametoall', sametoall, ...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext,'to_report', to_report);



end

% [=================================================================]
%   plot
% [=================================================================]
	V = simresults.networkHistory.V_soma;

	for c = 1:noneurons
		snip = @(t) V(c,t-w:t+w);		
		trigVonE{c} = cell2mat(arrayfun(snip, pert.triggers{1}(1:end-1), 'uniformoutput',0)');
		trigVonI{c} = cell2mat(arrayfun(snip, pert.triggers{2}(1:end-1), 'uniformoutput',0)');

		EbefI = find(pert.triggers{1}(1:end-1)-pert.triggers{2}(1:end-1)<0);
		IbefE = find(pert.triggers{1}(1:end-1)-pert.triggers{2}(1:end-1)>0);

		figure
		subplot(311)
		imagesc(trigVonE{c})
		colormap(cmap)
		subplot(312)
		plot(trigVonI{c}(IbefE,:)','b')
		title(num2str(c))
		subplot(313)
		plot(trigVonE{c}(EbefI,:)','r')
		title(num2str(c))
		% export_fig([num2str(c) '.png'])
		% close 


	end



resultstable = profile_sim(simresults);
R = resultstable;
sel_cel_idx = 1;
sel_fields = {'g_CaL', 'g_K_Ca', 'g_int', 'p1', 'p2', 'ampl', 'freq_each', 'maxV', 'meanVm'}
sel_table = R.allneurons(sel_cel_idx,sel_fields);
NDscatter(sel_table, 1)





