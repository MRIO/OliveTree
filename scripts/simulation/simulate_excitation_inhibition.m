% simulate_excitation_inhibition.m

% test significance of advance of neighboring neurons

% test significance of complex spike probability

prep_conditions    = 1;
compute_transients = 1; transienttime = 1000;
stimulate    = 1;
volumetric_activity = 0; % slow

savemovies = 0;
snapshot = 1;
gpu = 0;

rng(1, 'twister')
% use precomputed transients of 50x20
% x 5 cells
 % y 10 cells
  % z 20 cells

 to_report = {'V_soma','I_cx36', 'Ca2Plus', 'V_dend'};

% global 
dt =0.05;
fs = 1/dt;
sametoall = 0;

% parameterset = 'test2d';
% parameterset = 'test';
% parameterset = 'onlyexcitation';
parameterset = 'both';
% parameterset = 'onlyinhibition';
switch parameterset
	case 'test2d'

			simtime  = 600;	
			gap = 0.01;
			inputmasktype = 'random';
			% inputmasktype = 'all';
			% inputmasktype = 'column';
			% inputmasktype = 'mid_layer';
			% inputmasktype = 'centerneuron';

			depth = 3; breadth = 3; height = 3;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			% createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions)
			W = createW('all to all', netsize, [], gap, 1, 0, 0, 10, []);

			% W = zeros(noneurons); 
			% allcells = 1:prod(netsize);
			% centercell = ceil(netsize(2));
			% othercells = setdiff(allcells, centercell);
			% W(centercell, othercells) = 1*gap;
			% W = W+W';
			W_3d_trans.W = W*gap;

			
			def_neurons = createDefaultNeurons(noneurons);
				g_CaL = linspace(.3,1,noneurons);
				g_CaL = g_CaL(randperm(noneurons));
				% g_CaL(centercell) = 0.4;
				def_neurons.g_CaL = g_CaL;
				def_neurons.gbar_ampa_soma = linspace(0,.3,noneurons);
				def_neurons.gbar_gaba_soma = ones(noneurons,1)*.7;



			noise_level = [3 5 0 0]; % pA/ms per cell
			sametoall   = .5;
			
			firstpulse_inh = [200]; % start of stimulation
			excitation_lags = [-100:20:100];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [1]; % stimulation period in ms
			exc_n_of_pulses = 1;

			% inhibitory step
			gabarelease_soma = 1;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 1;
			
			% inhibitory step
			gabarelease_dend = 1;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 1;
			
			% ampa noise?
			ampa_noise = 0;
			 ampa_release_prob = .005*fs;

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			

		case 'test'

			simtime  = 300;	
			% inputmasktype = 'column';
			inputmasktype = 'mid_layer';
			depth = 3; breadth = 5; height = 10;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01;
			
			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.2];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',1); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;


			W_3d_trans = createW('3d_euclidean_rndwalk', netsize,2.5,gap, 1, 0);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, plotornot);

			noise_level = [3 5 0 0]; % pA/ms per cell
			
			firstpulse_inh = [50 100]; % start of stimulation
			excitation_lags = [-10 10];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [1]; % stimulation period in ms
			exc_n_of_pulses = 4;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			
			% inhibitory step
			gabarelease_dend = 0;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
			
			% 
			ampa_noise = 1;
			 ampa_release_prob = .005*fs;
			% syn_noise_ampa = 

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times
			

		case 'onlyexcitation'

			simtime  = 2000;	
			% inputmasktype = 'column';
			inputmasktype = 'mid_layer';
			depth = 5; breadth = 10; height = 20;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01; % mS/cm^2
			

			
			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.4];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',0); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;
				def_neurons.g_int = .1 + rand(noneurons,1)*.2;

			W_3d_trans = createW('3d_euclidean_rndwalk', netsize,1.5,gap, 1, 1);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, 1);

			noise_level = [3.5 3 1 1]; % pA/ms per cell - 3.5 3 1
			
			firstpulse_inh = [100]; % start of stimulation
			excitation_lags = [-80:10:80];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [1 2 4 8 12 16 28 20]; % stimulation period in ms
			exc_n_of_pulses = 3;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			
			% inhibitory step
			gabarelease_dend = 0;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
			
			% 
			ampa_noise = 0;
			 % ampa_release_prob = .005/(1/delta);

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times
			
		case 'onlyinhibition'

			simtime  = 2000;	
			% inputmasktype = 'column';
			% inputmasktype = 'dist_to_point';
			inputmasktype = 'mid_layer';
			depth = 5; breadth = 10; height = 20;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01;
			

			
			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.2];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',0); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;

			W_3d_trans = createW('3d_euclidean_rndwalk', netsize, 2, gap, 1, 0);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, 1);

			noise_level = [3.5 3 1]; % pA/ms per cell
			
			firstpulse_inh = [100 125 150]; % start of stimulation
			excitation_lags = 0;

			% excitatory step
			glurelease = 0;
			exc_pulse_dur = 1;
			exc_stim_T   = [1]; % stimulation period in ms
			exc_n_of_pulses = 3;

			% inhibitory step
			gabarelease_soma = 1;
			gaba_pulse_dur = 2;
			inh_stim_T   = [1 2 5 10 20 30 40 50]; % stimulation period in ms
			inh_n_of_pulses = 30;
			
			% inhibitory step
			gabarelease_dend = 1;
			gaba_pulse_dur = 1;
			% inh_stim_T   = 1; % stimulation period in ms
			% inh_n_of_pulses = 2;
			
			% 
			ampa_noise = 0;
			 ampa_release_prob = .005/20;

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times
	

		case 'both'

			simtime  = 1500;	
			% inputmasktype = 'column';
			inputmasktype = 'mid_layer';
			depth = 5; breadth = 10; height = 20;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01;
			

			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.2];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',0); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;


			W_3d_trans = createW('3d_euclidean_rndwalk', netsize,1.5,gap, 1, 0);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, 1);

			noise_level = [0 0 0 0]; % pA/ms per cell
			
			firstpulse_inh = [500:5:550]; % start of stimulation
			excitation_lags = [-80:10:80];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [2]; % stimulation period in ms
			exc_n_of_pulses = 2;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			
			% inhibitory step
			gabarelease_dend = 1;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
		
		% 
			ampa_noise = 0;
			 ampa_release_prob = .005/20

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times

	case 'both_random'

			simtime  = 1500;	
			% inputmasktype = 'column';
			inputmasktype = 'mid_layer';
			depth = 5; breadth = 10; height = 20;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01;
			

			
			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.2];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',0); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;
				

			W_3d_trans = createW('3d_euclidean_rndwalk', netsize,1.5,gap, 1, 0);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, 1);

			noise_level = [3.5 3 1 1]; % pA/ms per cell
			
			firstpulse_inh = [500:25:750]; % start of stimulation
			
			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [3]; % stimulation period in ms
			exc_n_of_pulses = 1;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			
			% inhibitory step
			gabarelease_dend = 1;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
			
			 
			ampa_noise = 1;
			 ampa_release_prob = .005/20


			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times


case 'both_tinynet'

			simtime  = 1500;	
			% inputmasktype = 'column';
			inputmasktype = 'mid_layer';
			depth = 3; breadth = 3; height = 3;
			netsize = [depth breadth height];
			noneurons = prod(netsize);
			
			gap = 0.01;
			

			def_neurons = createDefaultNeurons(noneurons)
				alp = 7;
				bet = 2;
				sup = [0.7 1.2];
				GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
				   noneurons, 'support', sup,'plot_distributions',0); 
				g_CaL = GCAL.sampleDraws{1}';				
				def_neurons.g_CaL = g_CaL;


			W_3d_trans = createW('3d_euclidean_rndwalk', netsize,1.5,gap, 1, 0);
			% W_3d = createW('3d_euclidean_rndwalk',netsize,0,gap, 1);

			noise_level = [3 5 0 0]; % pA/ms per cell
			
			firstpulse_inh = [500:10:550]; % start of stimulation
			excitation_lags = [-80:10:80];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [2]; % stimulation period in ms
			exc_n_of_pulses = 2;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			
			% inhibitory step
			gabarelease_dend = 1;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
		
		% 
			ampa_noise = 0;
			 ampa_release_prob = .005/20

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times


end



%    _   __                         __       
%   (_)_/_/   ____ ___  ____ ______/ /_______
%    _/_/    / __ `__ \/ __ `/ ___/ //_/ ___/
%  _/_/_    / / / / / / /_/ (__  ) ,< (__  ) 
% /_/ (_)  /_/ /_/ /_/\__,_/____/_/|_/____/  
                                           

if exist('inputmasktype')
	perturbation = create_input_mask(netsize,inputmasktype);
end


%                                                                    ___ __  _                 
%     ____  ________  ____  ____ _________     _________  ____  ____/ (_) /_(_)___  ____  _____
%    / __ \/ ___/ _ \/ __ \/ __ `/ ___/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%   / /_/ / /  /  __/ /_/ / /_/ / /  /  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
%  / .___/_/   \___/ .___/\__,_/_/   \___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
% /_/             /_/                                                                          

if prep_conditions

	conds = 1; clear condition


	for s = firstpulse_inh
		for lag = excitation_lags;
			for excT = exc_stim_T
				for inhT = inh_stim_T
				if gabarelease_dend
			        gaba_mask_dend = perturbation  ;
					gaba_onsets_dend = s + cumsum(repmat([0 inhT],1,inh_n_of_pulses)) ;
				else
					gaba_mask_dend = [];
					gaba_onsets_dend = [];
				end

				if gabarelease_soma
			        gaba_mask_soma = perturbation  ; %create_input_mask([5 10 20], 'dist_to_point');
					gaba_onsets_soma = s + cumsum(repmat([0 inhT],1,inh_n_of_pulses)) ;
				else
					gaba_mask_soma = [];
					gaba_onsets_soma = [];
				end

				if glurelease
			        glu_mask = perturbation;
					glu_onsets = s + cumsum(repmat([0 excT],1,exc_n_of_pulses)) + lag;
				else
					glu_onsets = [];
					glu_mask = [];

				end

				condition{conds}.description = 'inhibition and excitation';
				condition{conds}.perturbation_onsets{1} = [];
				condition{conds}.perturbation_mask{1} = [];
				condition{conds}.perturbation_type{1} = 'current';
				condition{conds}.perturbation_pulse_duration{1} = [0];

				condition{conds}.perturbation_onsets{2} = glu_onsets;
				condition{conds}.perturbation_mask{2} 	= glu_mask; % selection of stimulated neurons
				condition{conds}.perturbation_type{2} = 'ampa';
				condition{conds}.perturbation_pulse_duration{2} = exc_pulse_dur;
				
				condition{conds}.perturbation_onsets{3} = gaba_onsets_dend;
				condition{conds}.perturbation_mask{3} 	= gaba_mask_dend; % selection of stimulated neurons
				condition{conds}.perturbation_type{3} = 'gaba_dend';
				condition{conds}.perturbation_pulse_duration{3} =  gaba_pulse_dur;

				condition{conds}.perturbation_onsets{4} = gaba_onsets_soma;
				condition{conds}.perturbation_mask{4} 	= gaba_mask_soma; % selection of stimulated neurons
				condition{conds}.perturbation_type{4} = 'gaba_soma';
				condition{conds}.perturbation_pulse_duration{4} =  gaba_pulse_dur;					

				if ampa_noise
					condition{conds}.perturbation_onsets{5} = ampa_release_prob; 
					condition{conds}.perturbation_mask{5} 	= create_input_mask(netsize,'all'); % selection of stimulated neurons
					condition{conds}.perturbation_type{5} = 'ampa_noise';
					condition{conds}.perturbation_pulse_duration{5} =  1;					
				end

				condition{conds}.perturbation_amplitude = 1;

				condition{conds}.noise_level = noise_level; %pA per cell
				condition{conds}.offset = [];
				condition{conds}.gaussalpha = [];
				condition{conds}.simtime = simtime;

				% condition{conds}.g_CaL = g_CaL;
				% condition{conds}.W = W;

				condition{conds}.simtime = simtime;
	            
	            condition{conds}.I_app = 0;
                
	%}
	            conds = conds+1;
		        end
	        end
	    end
	end

	conds = [1:conds-1];
end

%                                     __          __                        _            __      
%   _________  ____ ___  ____  __  __/ /____     / /__________ _____  _____(_)__  ____  / /______
%  / ___/ __ \/ __ `__ \/ __ \/ / / / __/ _ \   / __/ ___/ __ `/ __ \/ ___/ / _ \/ __ \/ __/ ___/
% / /__/ /_/ / / / / / / /_/ / /_/ / /_/  __/  / /_/ /  / /_/ / / / (__  ) /  __/ / / / /_(__  ) 
% \___/\____/_/ /_/ /_/ .___/\__,_/\__/\___/   \__/_/   \__,_/_/ /_/____/_/\___/_/ /_/\__/____/  
%                    /_/                                                                         

if compute_transients
	
	transientnoise = 0;
	
    [transients] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',dt,'cell_parameters', def_neurons ,'W',W_3d_trans.W,'ou_noise', noise_level, 'sametoall',sametoall);


	
	transients.g_CaL = g_CaL;
	transients.condition.simtime = simtime;
	transients.condition.perturbation_mask = [];
	transients.W = W_3d_trans;
	
	transients.rows = breadth*depth;
	transients.columns = height;
    

    transients.condition.perturbation_amplitude = 0;
    transients.condition.perturbation_onsets = [];
    transients.condition.noise_level = transientnoise;
    transients.condition.offset = [0 0 0 0];
    transients.condition.gaussalpha = 0;
	transients.condition.perturbation_mask = [];
    transients.condition.I_app = 0;
    transients.condition.seed = rng(1, 'twister');

end


%          __  _                 __      __                               ___ __  _                 
%    _____/ /_(_)___ ___  __  __/ /___ _/ /____     _________  ____  ____/ (_) /_(_)___  ____  _____
%   / ___/ __/ / __ `__ \/ / / / / __ `/ __/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%  (__  ) /_/ / / / / / / /_/ / / /_/ / /_/  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
% /____/\__/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
                                                                                                  


flush = 1;
thisfig = figure;

if stimulate
	% for ccc = 22;
	for ccc = 1:length(condition);
		pert.mask  	   = condition{ccc}.perturbation_mask;
		pert.amplitude = condition{ccc}.perturbation_amplitude;
		pert.triggers  = condition{ccc}.perturbation_onsets;
		pert.duration  = condition{ccc}.perturbation_pulse_duration;
		pert.type	   = condition{ccc}.perturbation_type;
		

		sim3D = IOnet('networksize', netsize,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'tempState',transients.lastState,'W',W_3d_trans.W ,'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall);
		
		sim3D.perturbation = pert;
		sim3D.W = W_3d_trans;
        condition{ccc}.I_app = 'notwritinginput';
		sim3D.condition = condition{ccc};


			writeresults = 1;
			if writeresults
			sim = ccc;
		    
			
		    results{sim} = replayResults(sim3D, [1:2], 0,[],thisfig);
		    export_fig(num2str(sim)) 

		    orderparameter = measureGlobalSync(sim3D,[1:simtime],1);
		    export_fig(['sync' num2str(sim)]) 
		    
		    clf
		    

		    RESULTS(sim,1) = mean(mean(W_3d_trans.W(W_3d_trans.W~=0)));
		    RESULTS(sim,2) = noise_level(1);
		    % RESULTS(sim,3) = gamul;
		    % RESULTS(sim,4) = conntype;
		    RESULTS(sim,5 ) = results{sim}.popfrequency;
		    RESULTS(sim,6 ) = results{sim}.propspkneurons;
		    RESULTS(sim,7 ) = mean(results{sim}.spikespercell/transienttime*noneurons);
		    RESULTS(sim,8 ) = mean(results{sim}.medfreq(results{sim}.medfreq>0));
		    RESULTS(sim,9 ) = mean(sum(W_3d_trans.W))
			% RESULTS(sim,10) = median(W_3d_trans.stats.connections)
			% RESULTS(sim,11) = min(W_3d_trans.stats.connections)
			% RESULTS(sim,12) = max(W_3d_trans.stats.connections)
			RESULTS(sim,13) = orderparameter.stats.firstordersync(1);
			RESULTS(sim,14) = orderparameter.stats.secondordersync(1);
			RESULTS(sim,15) = orderparameter.stats.overallsync(1);

		    RESULTS(sim,16) = sim3D.failed;
		end



				    

		%               __                     __       _                    __  _       _ __       
		%  _   ______  / /_  ______ ___  ___  / /______(_)____   ____ ______/ /_(_)   __(_) /___  __
		% | | / / __ \/ / / / / __ `__ \/ _ \/ __/ ___/ / ___/  / __ `/ ___/ __/ / | / / / __/ / / /
		% | |/ / /_/ / / /_/ / / / / / /  __/ /_/ /  / / /__   / /_/ / /__/ /_/ /| |/ / / /_/ /_/ / 
		% |___/\____/_/\__,_/_/ /_/ /_/\___/\__/_/  /_/\___/   \__,_/\___/\__/_/ |___/_/\__/\__, /  
		%                                                                                  /____/   

		if volumetric_activity

						% which sim to display
						sim = sim3D;
					    animate_volume(sim,[1:simtime],savemovies)
					    if savemovies
							eval(['!mv volume.mp4 3d_pert_' num2str(condition{ccc}.perturbation_amplitude) 'pA.mp4'])
							close all
						end
		end


		if savemovies
			ons = condition{ccc}.perturbation_onsets;
			amps = condition{ccc}.perturbation_amplitude;

			% replayResults(sim3D,[1:simtime],0);
			% eval(['!mv sim.avi samp_' num2str(ons(1)) 'ms' num2str(amps) 'pA.mp4'])

			% close all
		end

		if snapshot
			results{ccc} = replayResults(sim3D, [1:2], 0)
		    % figure(1)
		    export_fig(num2str(ccc)) 
		    clf
		end


		if flush

			save(['sim3D_' num2str(ccc)],  'sim3D')
			clear sim3D
			% gpuDevice(1)
		end
	end


end
	
% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 

%    _   __    ______                        ______                                   __       __          
%   (_)_/_/   / ____/________  ____ ___     /_  __/___  ____ ___     ____ _____  ____/ /      / /___ _____ 
%    _/_/    / /_  / ___/ __ \/ __ `__ \     / / / __ \/ __ `__ \   / __ `/ __ \/ __  /  __  / / __ `/ __ \
%  _/_/_    / __/ / /  / /_/ / / / / / /    / / / /_/ / / / / / /  / /_/ / / / / /_/ /  / /_/ / /_/ / / / /
% /_/ (_)  /_/   /_/   \____/_/ /_/ /_/    /_/  \____/_/ /_/ /_/   \__,_/_/ /_/\__,_/   \____/\__,_/_/ /_/ 
                                                                                                         

% IO RESPONSES TO MEJ OR SUP PEDUNCLE STIM
% (rostral MAO and PO)

% BEFORE LESION

% 3 pulses Mesodiencephalon
% 4-8ms (~100% of neurons)
% 180ms (~50%)

% 3 pulses Superior Peduncle
% 9-15ms (66%)
% 180ms (~7%)

% AFTER LESION (no nucleo olivary pathway)

% 3 pulses Mesodiencephalon
% 4-8ms (~100% of neurons)
% 180ms (from ~50% (no lesion) to 75%)

% 3 pulses Superior Peduncle
% 9-15ms (66% - 33% = 33%)
% 180ms (from ~7% to 35%)

% short latency (MDJ) : 4ms
% short latency (SCP) : 9ms




