% simulate_3d_layer_perturbation.m

prep_conditions    = 1;
compute_transients = 1;
stimulate_layer    = 1;
volumetric_activity = 0;

savemovies = 0;



rng(1, 'twister')
% use precomputed transients of 50x20
% x 5 cells
 % y 10 cells
  % z 20 cells


parameterset = 'pset2';
switch parameterset
	case 'pset1'
		inputmasktype = 'random';
		depth = 5; breadth = 10; height = 20;
		noneurons = breadth*depth*height;
		
		% W = createW('no_torus', 32, 32, 0.005, 1);
		W_3d = createW('3d',[depth breadth height],0,0.005, 1);
		W = W_3d.W;

		% duration of the pulse
		stim_dur = 10; % du - frequency of stimulation (1000Hz pulses, 100Hz pulse duration)
		stim_T   = 2000; % stimulation period in ms
		simtime  = 170000;
		stim_off = 101000; % t of last perturbation 

		noise_level = 6; % pA per cell
		amplitudes = [-100 -50 50 100]; % TOTAL CURRENT TO THE STIMULATED NEURONS
		% stimulus times
		s1 = [1000]; % start of stimulation


	case 'pset1_test'
	
		inputmasktype = 'random';
		depth = 5; breadth = 10; height = 20;
		noneurons = breadth*depth*height;
		
		% W = createW('no_torus', 32, 32, 0.005, 1);
		W_3d = createW('3d',[depth breadth height],0,0.00, 1);
		W = W_3d.W;

		% duration of the pulse
		stim_dur = 5; % du - frequency of stimulation (1000Hz pulses, 100Hz pulse duration)
		stim_T   = 2000; % stimulation period in ms
		simtime  = 30000;
		stim_off = 101000; % t of last perturbation 

		noise_level = 4; % pA per cell
		amplitudes = [-50 50]; % TOTAL CURRENT TO THE STIMULATED NEURONS
		% stimulus times
		s1 = [1000]; % start of stimulation

	case 'pset2'
		breadth = 5; depth = 10; height = 20;
		noneurons = breadth*depth*height;
		inputmasktype = 'random';

		gaussalpha = [1.5 1.5 8.5 ]; % conditions for whole layer stimulation
		% gaussalpha = [4.5 2.5 100 ]; % conditions for subset stimulation
		offset = [ 0 0 0];

		% W = createW('no_torus', 32, 32, 0.005, 1);
		W_3d = createW('3d',[depth breadth height],0,0.005, 1);
		W = W_3d.W;

		% duration of the pulse
		stim_dur = 10; % du 
		stim_T   = 10; % stimulation period in ms
		simtime = 1000;
		stim_off = 100; % t of last perturbation 

		noise_level = 0; % pA per cell
		amplitudes = [-40]; % TOTAL CURRENT TO THE STIMULATED NEURONS
		% stimulus times
		s1 = [100]; %

		% alpha  = [2 ]; %edges for alpha parameter
		% beta   = [2 ];       %edges for beta parameter
		% res = 1;
		% sup = [0.13 0.17];
		% scale = 1;

		% GINT = BetaDistributions('alpha', alpha, 'beta', beta, 'no_draws',depth*height*breadth, 'support', sup,'plot_distributions',1); 
		% g_int = GINT.sampleDraws{1}';
		g_int = .13;


		% g_KCa  = 35;
		alpha  =2;
		beta   = 2;
		res = 1;
		sup = [30 40];
		scale = 1;
		GKCa = BetaDistributions('alpha', alpha, 'beta', beta, 'no_draws',depth*height*breadth, 'support', sup,'plot_distributions',1); 
		g_KCa = GKCa.sampleDraws{1}';

		
end


dt =0.05;
fs = 1/dt;

if ~exist('g_CaL','var')
	alpha  = [32.5]; %edges for alpha parameter
	beta   = [12];       %edges for beta parameter

	scale = 1.35; % 1.35
figure
	GCAL = BetaDistributions('alpha', alpha, 'beta', beta, 'no_draws',depth*height*breadth, 'scale', scale,'plot_distributions',false); 
	g_CaL = GCAL.sampleDraws{1}';
	distribution_parameters = GCAL.parameters{1};
end


% create masks for input
switch inputmasktype

	case 'mid_layer'

			perturbation = zeros(depth, breadth, height);
			perturbation(1:depth, 1:breadth, ceil(height/2)) = 1;

			perturbation = reshape(perturbation, breadth*height*depth,1); 

	case 'random'

			perturbation = zeros(depth, breadth, height);
			% HALF RECTIFIED GAUSSIAN
			perturbation(1:depth, 1:breadth, ceil(height/2)) =  abs(randn(depth,breadth)) >.5;

			perturbation = reshape(perturbation, breadth*height*depth,1); 


	case 'gauss'

			G = gausswin(breadth, gaussalpha(1))*gausswin(depth, gaussalpha(2))';
			g(1,1,:) = gausswin(height, gaussalpha(3));
			perturbation_gaus = bsxfun(@times,G, g );
			perturbation_norm = perturbation_gaus/sum(perturbation_gaus(:));
			perturbation_shift= circshift(perturbation_norm,offset);
			perturbation = reshape(perturbation_shift, breadth*height*depth,1); % scale to deliver a total of A pA;
end


%                                                                    ___ __  _                 
%     ____  ________  ____  ____ _________     _________  ____  ____/ (_) /_(_)___  ____  _____
%    / __ \/ ___/ _ \/ __ \/ __ `/ ___/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%   / /_/ / /  /  __/ /_/ / /_/ / /  /  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
%  / .___/_/   \___/ .___/\__,_/_/   \___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
% /_/             /_/                                                                          

if prep_conditions

	conds = 1; clear conditions

	for a = amplitudes
		for s = s1

			
	        perturbation_mask = abs(perturbation)> 0;
		
			pulses = [1:stim_dur]; % current pulses 
			onsets = [s:stim_T:(simtime-1)];

			% [pp oo] = meshgrid(pulses, onsets);

			% perturbation_pulses = pp' + oo';
			% perturbation_pulses = perturbation_pulses(:);

			pp = 1;
			for onse = onsets
				for puls = pulses
					perturbation_pulses(pp) = onse+puls;
					pp = pp+1;
				end
			end


			perturbation_triggers = perturbation_pulses*fs;
	        
				condition{conds}.perturbation_amplitude = a; % current per stimulated neuron in mask
				condition{conds}.perturbation_onsets = perturbation_pulses; % periodicity in ms
				condition{conds}.perturbation_mask = perturbation; % selection of stimulated neurons

				condition{conds}.noise_level = noise_level; %pA per cell
				condition{conds}.offset = [];
				condition{conds}.gaussalpha = [];
				condition{conds}.simtime = simtime;
				

				condition{conds}.simtime = simtime;
	            
	            condition{conds}.I_app = 0;
                condition{conds}.netsize = [depth breadth height]; 
	%}
	            conds = conds+1;
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
	transienttime = 1000;
	noise_level = 0;
	% [transients] = IOtoplevelCa_gpu('rows',breadth,'columns',depth*height,'appCurrent',noise,'time',simtime,'delta',dt,'g_CaL',g_CaL ,'W',W);
	% [transients] = IOtoplevelCa_gpu('rows',depth,'columns',breadth*height,'appCurrent',noise,'time',500,'delta',dt,'g_CaL',g_CaL ,'W',W);
    % [transients] = IOnet('rows',depth,'columns',breadth*height,'appCurrent',0,'time',transienttime,'delta',dt,'g_CaL',g_CaL ,'W',W,'ou_noise', [noise_level 0 1], 'somaDendriteLeak', g_int, 'dendriticCalciumActivatedPotassium', g_KCa);
	[transients] = IOnet('rows',depth,'columns',breadth*height,'appCurrent',0,'time',transienttime,'delta',dt,'g_CaL',g_CaL ,'W',W,'ou_noise', [noise_level 0 1]);


	
	transients.g_CaL = g_CaL;
	transients.condition.simtime = simtime;
	transients.condition.perturbation_mask = [];
	transients.condition.W = W;
	transients.condition.distribution_parameters = distribution_parameters;
	
	transients.rows = breadth*depth;
	transients.columns = height;
    

    transients.condition.perturbation_amplitude = 0;
    transients.condition.perturbation_onsets = [];
    transients.condition.noise_level = noise_level;
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

if stimulate_layer
	% for ccc = 22;
	for ccc = 1:length(condition);
		pert.mask  	   = condition{ccc}.perturbation_mask;
		pert.amplitude = condition{ccc}.perturbation_amplitude;
		pert.triggers  = condition{ccc}.perturbation_onsets;
		pert.type		= 'rectified_shifted_gauss';
		% pert.type		= 'same_to_all';

		sim3D = 			 IOnet('rows',depth,'columns',breadth*height,'time',simtime,'delta',dt,'g_CaL',transients.g_CaL,'tempState',transients.lastState,'W',W ,'ou_noise',[noise_level, 0, 1] , 'perturbation', pert);
		% sim3D = IOtoplevelCa_gpu('rows',depth,'columns',breadth*height,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',transients.g_CaL,'tempState',transients.lastState,'W',W);

        condition{ccc}.I_app = 'notwritinginput';
		sim3D.rows = breadth*depth;
		sim3D.columns = height;
		sim3D.condition = condition{ccc};

		%               __                     __       _                    __  _       _ __       
		%  _   ______  / /_  ______ ___  ___  / /______(_)____   ____ ______/ /_(_)   __(_) /___  __
		% | | / / __ \/ / / / / __ `__ \/ _ \/ __/ ___/ / ___/  / __ `/ ___/ __/ / | / / / __/ / / /
		% | |/ / /_/ / / /_/ / / / / / /  __/ /_/ /  / / /__   / /_/ / /__/ /_/ /| |/ / / /_/ /_/ / 
		% |___/\____/_/\__,_/_/ /_/ /_/\___/\__/_/  /_/\___/   \__,_/\___/\__/_/ |___/_/\__/\__, /  
		%                                                                                  /____/   

		if volumetric_activity

						% which sim to display
						sim = sim3D;
					    animate_volume(sim,[900:6000],savemovies)
					    if savemovies
							eval(['!mv volume.mp4 3d_rhythm_long' num2str(condition{ccc}.perturbation_amplitude) 'pA.mp4'])
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


		if flush
			save(['sim3D_' num2str(condition{ccc}.perturbation_amplitude) '_' num2str(condition{ccc}.perturbation_onsets(1)) 'ms_' num2str(condition{ccc}.noise_level) 'pA_long'],  'sim3D')
			clear sim3D
			gpuDevice(1)
		end
	end


end
	
% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 




