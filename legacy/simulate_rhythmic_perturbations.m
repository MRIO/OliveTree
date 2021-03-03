% 	   _   __                                              __     __  __    _                     __     _ ____               __     _                               _      _                                           __           
% 	  (_)_/_/   _________  ____ ___  ____ ___  ___  ____  / /_   / /_/ /_  (_)____   ____  __  __/ /_   (_) __/  ____  ____  / /_   (_)___     ____ ___  ____ ______(_)___ ( )_____   _________  ____ ___  ____  __  __/ /____  _____
% 	   _/_/    / ___/ __ \/ __ `__ \/ __ `__ \/ _ \/ __ \/ __/  / __/ __ \/ / ___/  / __ \/ / / / __/  / / /_   / __ \/ __ \/ __/  / / __ \   / __ `__ \/ __ `/ ___/ / __ \|// ___/  / ___/ __ \/ __ `__ \/ __ \/ / / / __/ _ \/ ___/
% 	 _/_/_    / /__/ /_/ / / / / / / / / / / /  __/ / / / /_   / /_/ / / / (__  )  / /_/ / /_/ / /_   / / __/  / / / / /_/ / /_   / / / / /  / / / / / / /_/ / /  / / /_/ / (__  )  / /__/ /_/ / / / / / / /_/ / /_/ / /_/  __/ /    
% 	/_/ (_)   \___/\____/_/ /_/ /_/_/ /_/ /_/\___/_/ /_/\__/   \__/_/ /_/_/____/   \____/\__,_/\__/  /_/_/    /_/ /_/\____/\__/  /_/_/ /_/  /_/ /_/ /_/\__,_/_/  /_/\____/ /____/   \___/\____/_/ /_/ /_/ .___/\__,_/\__/\___/_/     
% 	                                                                                                                                                                                                   /_/                           
if ~exist('default','var')
	addpath(genpath('/Users/M/Projects/Code/MatlabExtra/plots'))
	set(0, 'defaultfigurecolor', ones(3,1))
	set(0, 'defaultfigurecolormap', linspecer(64))
	
	default = 1;
end

if ~exist('s','var')
	s = 0;
end



% 	         _                   __                          
% 	   _____(_)___ ___  _____   / /_____     _______  ______ 
% 	  / ___/ / __ `__ \/ ___/  / __/ __ \   / ___/ / / / __ \
% 	 (__  ) / / / / / (__  )  / /_/ /_/ /  / /  / /_/ / / / /
% 	/____/_/_/ /_/ /_/____/   \__/\____/  /_/   \__,_/_/ /_/ 
% 	                                                         
% execute these scripts:
transients   = 1;
shot_noise   = 0;
	continue_sim = 1 ; % means that simIO_transients is in the workspace and will be used --  maybe we discover chaos, as lorentz did.
mix_stim = 0;
	continue_sim = 1 ;

% 	       __     ____            ____                                           __                
% 	  ____/ /__  / __/___ ___  __/ / /_   ____  ____ __________ _____ ___  ___  / /____  __________
% 	 / __  / _ \/ /_/ __ `/ / / / / __/  / __ \/ __ `/ ___/ __ `/ __ `__ \/ _ \/ __/ _ \/ ___/ ___/
% 	/ /_/ /  __/ __/ /_/ / /_/ / / /_   / /_/ / /_/ / /  / /_/ / / / / / /  __/ /_/  __/ /  (__  ) 
% 	\__,_/\___/_/  \__,_/\__,_/_/\__/  / .___/\__,_/_/   \__,_/_/ /_/ /_/\___/\__/\___/_/  /____/  
% 	                                  /_/                                                          


% default parameters

dt = .05; fs = 1/dt;
rows = 20; % 50*15
columns = 20;
gap_c = 0.00; % mS/cm^2
if ~exist('g_CaL','var')
	g_CaL = 1.7+ -1.8*betarnd(3,8,rows*columns,1)'; % as discussed with Jornt de Gruijl
end



% set for the whole workspace
set(0, 'defaultaxescolororder', linspecer(rows*columns))


% prelude

if transients
	simtime = 400; %ms 
	noise_level = 2; % pA per cell (* a normal, mu = 0, alpha = .5)
	noise = randn(rows*columns,fs*simtime)*noise_level;
	I_app = noise;

	% [simIO_gpu_transients] = IOtoplevelCa_M('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',g_CaL,'conductance',gap_c);
	[simIO_gpu_transients] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',g_CaL,'g_Gap',gap_c,'connections','no_torus');

	simIO_gpu_transients.g_CaL = g_CaL;
	simIO_gpu_transients.condition.perturb_amplitude = [0];
	simIO_gpu_transients.condition.offset = [0];
	simIO_gpu_transients.condition.perturb_onsets = [0];
	simIO_gpu_transients.condition.simtime = simtime;
	simIO_gpu_transients.condition.perturbation_map = [];
	simIO_gpu_transients.condition.noise_level = noise_level;
end


if shot_noise %impulse perturbation

	%read conditions from this file:
	create_conditions

	for sm = conds
		% gpuDevice(1)
		simtime = condition{sm}.simtime; %ms
		noise_level = condition{sm}.noise_level;
		perturb_triggers = condition{sm}.perturb_onsets;
		gap_c = condition{sm}.g_Gap;
		
		
			I_app = zeros(rows*columns,fs*simtime);
				perturbation_rand = rand(rows* columns,1);
				
					perturbation = condition{sm}.perturbation_map;
				
					perturbation = circshift(perturbation,[condition{sm}.offset 0]);

				noise = randn(rows*columns,fs*simtime)*noise_level;

			% perturb_triggers = [f:f+2 (100:1000/perturb_freq:simtime)]*fs;
			perturb_on = perturb_triggers*fs;
			I_app(:,perturb_on) = repmat(perturbation,[1,length(perturb_triggers)]);

			I_app = I_app + noise;


			if continue_sim
				[simIO_gpu{sm}] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',simIO_gpu_transients.g_CaL,'tempState',simIO_gpu_transients.lastState,'g_Gap',gap_c);
				simIO_gpu{sm}.condition = condition{sm};
			else
				[simIO_gpu] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',g_CaL);
			end
	end

	% [simIO_M] = IOtoplevelCa_M('rows',rows,'columns', columns,'appCurrent',I_app,'time',time,'delta',dt,'g_CaL',g_CaL);
end

% no stimulation, find peaks larger than -52 - find spread

% stimulation, find peaks larger than -52 in the next oscillation - find spread


if mix_stim


	condition{1}.perturb_amplitude = [5000 -5000];
	condition{1}.offset = [0 0];
	condition{1}.perturb_onsets = [90 110];

	condition{2}.perturb_amplitude = [-5000 5000];
	condition{2}.offset = [0 0];
	condition{2}.perturb_onsets = [90 110];

	condition{3}.perturb_amplitude = [5000 -5000];
	condition{3}.offset = [-200 200];
	condition{3}.perturb_onsets = [90 110];

	condition{4}.perturb_amplitude = [-5000 5000];
	condition{4}.offset = [-200 200];
	condition{4}.perturb_onsets = [90 110];


	for cond = 1:4

		simtime = 400; %ms
		perturb_amplitude_1 = condition{cond}.perturb_amplitude(1); %uA per group
		perturb_amplitude_2 = condition{cond}.perturb_amplitude(2); %uA per group
		noise_level = 5;
		
		pert1_on = condition{cond}.perturb_onsets(1);
		pert2_on = condition{cond}.perturb_onsets(2); 

		offset_1 = condition{cond}.offset(1);
		offset_2 = condition{cond}.offset(2);

		perturb_freq = 10; %Hz

		s = s+1;
			
		I_app = zeros(rows*columns,fs*simtime);
			perturbation_rand = rand(rows* columns,1);
			perturbation_gaus = reshape(gausswin(rows)*gausswin(columns)', rows*columns,1);
			perturbation_mask = reshape(perturbation_gaus>.95, rows*columns, 1);
		
				perturbation = perturbation_mask.*perturbation_gaus;
				perturbation = perturbation/sum(perturbation);
				
				perturbation_1 = circshift(perturbation,[offset_1 0])*perturb_amplitude_1;
				perturbation_2 = circshift(perturbation,[offset_2 0])*perturb_amplitude_2;
			
			noise = randn(rows*columns,fs*simtime)*noise_level;
		
		perturb_triggers_1 = [pert1_on:pert1_on+2 ]*fs;
		perturb_triggers_2 = [pert2_on:pert2_on+2 ]*fs;
		% perturb_triggers = [f:f+2 (100:1000/perturb_freq:simtime)]*fs;
		
		I_app(:,perturb_triggers_1) = repmat(perturbation_1,[1,length(perturb_triggers_1)]);
		I_app(:,perturb_triggers_2) = repmat(perturbation_2,[1,length(perturb_triggers_2)]);
		I_app = I_app + noise;


		if continue_sim
			[simIO] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',simIO_gpu_transients.g_CaL,'tempState',simIO_gpu_transients.lastState);
			% simIO_mix_depo{s}.perturb_freq = perturb_freq;
			simIO_mix_depo{s} = simIO;
			simIO_mix_depo{s}.condition = condition{cond};
			simIO_mix_depo{s}.perturbation_map = perturbation_1+perturbation_2;
			simIO_mix_depo{s}.noise_level = noise_level;


		else
			[simIO_mix_depo] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',g_CaL);
		end

	end
	
end


% no stimulation, find peaks larger than -52 - find spread

% stimulation, find peaks larger than -52 in the next oscillation - find spread

% 






































