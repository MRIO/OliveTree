% simulate_3d_layer_perturbation_images.m

prep_conditions    = 0;
compute_transients = 0;
stimulate_layer    = 1;
display_results    = 0;
volumetric_activity = 1;

savemovies = 1;

	IM = imread('Sierpinski_carpet.gif');

rng(1, 'twister')



% sim parameters
dt =0.025;
fs = 1/dt;

fieldstosave = {'V_soma'};%, 'V_dend'};
cell_function = 'vanilla';

steadystatetime = 40; simtime = 50;

% other parameters
gaps 	  = .01;
connset   = [3];   
radiuses  = [2.5 10];
alpha    = [7];    % alpha parameter
bet      = [3];    % beta parameter

betasup  = [.5 1.2]% support for the 

% net size
depth = 2; breadth = 25; height = 25;
noneurons = breadth*depth*height;
netsize = [depth breadth height];

% default neurons
def_neurons = createDefaultNeurons(noneurons);
	    GCAL = BetaDistributions('alpha', alpha, 'beta', bet, 'no_draws', ...
       noneurons, 'support', betasup,'plot_distributions',0);
       g_CaL = GCAL.sampleDraws{1}';
       distribution_parameters = GCAL.parameters{1};
       def_neurons.g_CaL = g_CaL;
       

% create image map

	figure

    if size(IM,3)>1
        IMM = squeeze(sum(IM,3));
    else 
        IMM = IM;
    end
	IMM(IMM==765)=0;
	IMM(IMM>0)=1;
	IMMR = imresize(IMM, [breadth ,height]);

	INPUT = permute(repmat(IMMR,[1, 1, depth]) , [3 1 2]);
	INPUT = flipdim(INPUT,3);
	
	h = vol3d('cdata',INPUT,'texture','3D');
	view(3)
	axis off
	axis tight;  daspect([1 1 1])


% create masks

	perturbation = double(INPUT);%/sum(INPUT(:));
	perturbation = reshape(perturbation, breadth*height*depth,1); % scale to deliver a total of A pA;
    perturbation_map = double(abs(perturbation)>quantile(abs(perturbation),.80));


% pulses

	firstpulse_inh = [3]; % start of stimulation
	excitation_lags = [0:20:100];


% create net
	radius = 2.5; gaps = 0.01;
	% W_3d_trans = createW('3d', netsize,2,radius, gaps, 1);
	W_3d = createW('3d',netsize,1,gaps, 1);
	W = W_3d.W;


% noise input

	noise_level = [3 2 0 0]; % pA/ms per cell
	sametoall   = .1;
	

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [1]; % stimulation period in ms
			exc_n_of_pulses = 1;

			% inhibitory step
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulatedulation period in ms
			inh_n_of_pulses = 1;
			
			% inhibitory step
			gabarelease_dend = 0;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 1;
			
			% ampa noise?
			ampa_noise = 0;
			 ampa_release_prob = .005*fs;

			amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS




%                     _     
%    ____ ___  ____ _(_)___ 
%   / __ `__ \/ __ `/ / __ \
%  / / / / / / /_/ / / / / /
% /_/ /_/ /_/\__,_/_/_/ /_/ 
                          

conds = 0; clear condition


for s = firstpulse_inh
	for lag = excitation_lags;
		for excT = exc_stim_T
			for inhT = inh_stim_T
				conds = conds+1;

				if gabarelease_dend
			        gaba_mask_dend = perturbation  ;
					gaba_onsets_dend = s + cumsum(repmat([0 inhT],1,inh_n_of_pulses)) ;
					condition{conds}.perturbation_onsets{3} = gaba_onsets_dend;
					condition{conds}.perturbation_mask{3} 	= gaba_mask_dend; % selection of stimulated neurons
					condition{conds}.perturbation_type{3} = 'gaba_dend';
					condition{conds}.perturbation_pulse_duration{3} =  gaba_pulse_dur;

				else
					gaba_mask_dend = [];
					gaba_onsets_dend = [];
				end

				if gabarelease_soma
			        gaba_mask_soma = perturbation  ; %create_input_mask([5 10 20], 'dist_to_point');
					gaba_onsets_soma = s + cumsum(repmat([0 inhT],1,inh_n_of_pulses)) ;
					condition{conds}.perturbation_onsets{4} = gaba_onsets_soma;
					condition{conds}.perturbation_mask{4} 	= gaba_mask_soma; % selection of stimulated neurons
					condition{conds}.perturbation_type{4} = 'gaba_soma';
					condition{conds}.perturbation_pulse_duration{4} =  gaba_pulse_dur;					

				else
					gaba_mask_soma = [];
					gaba_onsets_soma = [];
				end

				if glurelease
			        glu_mask = perturbation;
					glu_onsets = s + cumsum(repmat([0 excT],1,exc_n_of_pulses)) + lag;
					condition{conds}.perturbation_onsets{2} = glu_onsets;
					condition{conds}.perturbation_mask{2} 	= glu_mask; % selection of stimulated neurons
					condition{conds}.perturbation_type{2} = 'ampa';
					condition{conds}.perturbation_pulse_duration{2} = exc_pulse_dur;
					
				
				else
					glu_onsets = [];
					glu_mask = [];

				end

				

				if ampa_noise
					condition{conds}.perturbation_onsets{5} = ampa_release_prob; 
					condition{conds}.perturbation_mask{5} 	= create_input_mask(netsize,'all'); % selection of stimulated neurons
					condition{conds}.perturbation_type{5} = 'ampa_noise';
					condition{conds}.perturbation_pulse_duration{5} =  1;					
				end



	            
	        end
        end
    end
end

conds = [1:conds-1];


% transients

transientnoise = 0;
	
    [transients] = IOnet( 'networksize', netsize ,'appCurrent',0,'time',steadystatetime,'delta',dt,'cell_parameters', def_neurons ,'W',W,'ou_noise', noise_level, 'sametoall',sametoall);



flush = 1;

if stimulate_layer
	% for ccc = 22;
	for ccc = 1:length(condition);

		for ccc = 1:length(condition);
			pert.mask  	   = condition{ccc}.perturbation_mask;
			pert.amplitude = 1;
			pert.triggers  = condition{ccc}.perturbation_onsets
			pert.duration  = condition{ccc}.perturbation_pulse_duration;
			pert.type	   = condition{ccc}.perturbation_type;
			

			sim3D = IOnet('networksize', netsize,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'tempState',transients.lastState,'W',W ,'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall);
			
			sim3D.perturbation = pert;
			sim3D.condition = condition{ccc};




		if volumetric_activity

			% which sim to display
			sim = sim3D;
		    animate_volume(sim3D,[1:simtime],savemovies,0)
		    if savemovies
				eval(['!mv volume.mp4 3d_pert_' num2str(ccc) 'pA.mp4'])
				close all
			end
		end


		% if savemovies
		% 	ons = condition{ccc}.perturbation_onsets;
		% 	amps = condition{ccc}.perturbation_amplitude;

		% 	% replayResults(sim3D,[1:simtime],0);
		% 	% eval(['!mv sim.avi samp_' num2str(ons(1)) 'ms' num2str(amps) 'pA.mp4'])

		% 	% close all
		% end


			if flush
				save(['sim3D_image' num2str(ccc)],  'sim3D')
				clear sim3D
				gpuDevice(1)
			end
		end	
	end

end
	
  
%               __                     __       _                    __  _       _ __       
%  _   ______  / /_  ______ ___  ___  / /______(_)____   ____ ______/ /_(_)   __(_) /___  __
% | | / / __ \/ / / / / __ `__ \/ _ \/ __/ ___/ / ___/  / __ `/ ___/ __/ / | / / / __/ / / /
% | |/ / /_/ / / /_/ / / / / / /  __/ /_/ /  / / /__   / /_/ / /__/ /_/ /| |/ / / /_/ /_/ / 
% |___/\____/_/\__,_/_/ /_/ /_/\___/\__/_/  /_/\___/   \__,_/\___/\__/_/ |___/_/\__/\__, /  
%                                                                                  /____/   








if display_results

	% which sim to display
	sim = sim3D;
    animate_volume(sim,[],savemovies)

end


% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 




