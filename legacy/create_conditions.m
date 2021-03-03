
clear conds
clear condition

stim_dur = 10; % duration of the pulse
stim_T   = 500; % stimulation period in ms
simtime = 500; % overall simulation time
stim_off = 30; % the last perturbation onset time

noise_level = 2; % pA per cell
amplitudes = [-1000 1000]; 
gaps = [0.04 0.001];
s1 = [10:10:stim_off];

alpha = 1.5; % the inverse of the std determining the shape of the gaussian noise kernel


% s1 = 20;

% projection field size:
% proj_size of the gaussian stimulus measured in units of IO somata 
proj_size = [8 8]; % must be odd or even similar to net_size 

conds = 1;

for a = amplitudes
	for g = gaps
		for s = s1
		
			perturbation_gaus = reshape(padarray(gausswin(proj_size(1),alpha)*gausswin(proj_size(2),alpha)', [(rows-proj_size(1))/2,  (columns-proj_size(2))/2], 'both'), rows*columns,1);
				% perturbation_mask = reshape(perturbation_gaus>.95, rows*columns, 1);
				perturbation_gaus = perturbation_gaus/sum(sum(perturbation_gaus));
						
				% perturbation = perturbation_mask.*perturbation_gaus;
				% perturbation = perturbation/sum(perturbation);

				pulses = [1:stim_dur]; % current pulses 
				onsets = [s:stim_T:(simtime-1)];

				[pp oo] = meshgrid(pulses, onsets);

				perturb_pulses = pp' + oo';
				perturb_pulses = perturb_pulses(:);

				condition{conds}.perturb_amplitude = a; % average current per stimulated neuron in mask
				condition{conds}.perturb_onsets = perturb_pulses; % periodicity in ms

				condition{conds}.noise_level = noise_level; %pA per cell

				condition{conds}.offset = [0 0];
				condition{conds}.simtime = simtime;
								
				condition{conds}.g_Gap =  g;

				perturbation = perturbation_gaus * a; %*
				condition{conds}.perturbation_map = perturbation;
				
				conds = conds+1;			
		end

	end
	
end
conds = [1:conds-1];



% % groups
% stim_dur = 10; % duration of the pulse
% stim_T   = 400; % stimulation period in ms
% simtime = 250; 
% stim_off = 110;

% noise_level = 2;
% amplitudes = [-1000 1000]; % >-100 breaks the net
% gaps = [0.02 0.01 0.001];
% % gaps = [0.04 0.02 0.01 0.001];
% % gaps = [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
% % s1 = [10:20:100];  % first perturbation
% s1 = [10:20:stim_off];

% alpha = 1.5;


% % s1 = 20;

% proj_size = [8 8]; % must be odd or even similar to net_size 

% conds = 1;


% for a = amplitudes
% 	for g = gaps
% 		for s = s1
		
% 			perturbation_gaus = reshape(padarray(gausswin(proj_size(1),alpha)*gausswin(proj_size(2),alpha)', [(rows-proj_size(1))/2,  (columns-proj_size(2))/2], 'both'), rows*columns,1);
% 				% perturbation_mask = reshape(perturbation_gaus>.95, rows*columns, 1);
% 				perturbation_gaus = perturbation_gaus/sum(sum(perturbation_gaus));
						
% 				% perturbation = perturbation_mask.*perturbation_gaus;
% 				% perturbation = perturbation/sum(perturbation);

% 				pulses = [1:stim_dur]; % current pulses 
% 				onsets = [s:stim_T:(simtime-1)];

% 				[pp oo] = meshgrid(pulses, onsets);

% 				perturb_pulses = pp' + oo';
% 				perturb_pulses = perturb_pulses(:);

% 				condition{conds}.perturb_amplitude = a; % average current per stimulated neuron in mask
% 				condition{conds}.perturb_onsets = perturb_pulses; % periodicity in ms

% 				condition{conds}.noise_level = noise_level; %pA per cell

% 				condition{conds}.offset = [0 0];
% 				condition{conds}.simtime = simtime;
								
% 				condition{conds}.g_Gap =  g;

% 				perturbation = perturbation_gaus * a; %*
% 				condition{conds}.perturbation_map = perturbation;
				
% 				conds = conds+1;			
% 		end

% 	end
	
% end
% conds = [1:conds-1];



% % groups