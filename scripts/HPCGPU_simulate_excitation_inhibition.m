% simulate_excitation_inhibition.m

% test significance of advance of neighboring neurons

% test significance of complex spike probability

prep_conditions    = 1;
stimulate_layer    = 1; 
volumetric_activity = 0;

% script parameters
plotresults = 0; writeresults = 0;
savemovies = 0;
snapshot = 0;
gpu = 1;


% simulation parameters
% global 
dt =0.025;
fs = 1/dt;

% network structure
depth = 3; breadth = 10; height = 20;
netsize = [depth breadth height];
noneurons = prod(netsize);
connradi = 2.5;

	scaling = 1;
	% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions)
	connections = createW('3d_euclidean_rndwalk', netsize, connradi, scaling, 1, plotresults, 5, 10);

% cell parameters 
def_neurons = createDefaultNeurons(noneurons);
	
	alp = 7;
	bet = 3;
	sup = [0.5 1.1];
	GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
	   noneurons, 'support', sup,'plot_distributions',0); 
	g_CaL = GCAL.sampleDraws{1}';				
	def_neurons.g_CaL = GCAL.sampleDraws{1}';
	def_neurons.g_h   = .12 + rand(noneurons,1)*.24;;
	def_neurons.g_int = .10 + rand(noneurons,1)*.04;
	def_neurons.g_ls   = .013 + rand(noneurons,1)*.002;
	def_neurons.g_K_Ca   = 45 + rand(noneurons,1)*20;
	def_neurons.gbar_ampa_soma = .05 + rand(noneurons,1)*.2;


	def_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; % linspace(0,.3,noneurons); %
	def_neurons.gbar_gaba_soma = ones(noneurons,1)*1;

	% compensation for gap junctions
	gap_neurons = def_neurons;
	gap_neurons.g_CaL  = def_neurons.g_CaL + .1;
	gap_neurons.g_int  = def_neurons.g_int -0.04;
	gap_neurons.g_ld   = def_neurons.g_ld - 0.003;
	gap_neurons.g_K_Ca = def_neurons.g_K_Ca + 10;



% noise seeds
	rng(1, 'twister')
	noiseseeds = [1 2 3];
	sametoall = .1;


%default for simulations underneath. MAYBE OVERWRITTEN!
cell_function = 'vanilla';
to_report = {'V_soma', 'I_cx36'};
to_report = {'V_soma'};


% parameterset = 'test2d';
% parameterset = 'test';
parameterset = 'onlyexcitation'
% parameterset = 'spontaneous'
% parameterset = 'onlyinhibition';
% parameterset = 'both';

switch parameterset

	case 'spontaneous'

			simtime  = 240000;	
			
			gaps = [0 0.01 0.025]; % mS/cm^2
			
			noise_level = [3 5 0 0]; % [the sig mu seed] pA/ms per cell - 3.5 3 1
			
			firstpulse_inh = [1500]; % start of stimulation
			excitation_lags = [0];

			% excitatory step
			glurelease = 0;
			exc_pulse_dur = 1;
			exc_stim_T   = [0]; % stimulation period in ms
			exc_n_of_pulses = 0;

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

			amplitudes = 2; % TOTAL CURRENT TO THE STIMULATED NEURONS
			% stimulus times


		case 'onlyexcitation'

			simtime  = 3000;
			% inputmasktype = 'column';
			% inputmasktype = 'mid_layer';
			inputmasktype = 'random';

			gaps = [0.025];

			noise_level = [3 9 0 0]; % [the sig mu seed] pA/ms per cell - 3.5 3 1
            
			firstpulse_inh = [1500]; % first stimulation
			excitation_lags = [1];

			% excitatory step
			glurelease = 1;
			exc_pulse_dur = 1;
			exc_stim_T   = [1000]; % stimulation period in ms
			exc_n_of_stimuli = 100;
			exc_stim = [0 5 10];

			exc_stimulus_train = bsxfun(@plus, cumsum(exc_stim_T*ones(exc_n_of_stimuli,1)), exc_stim);
			exc_stimulus_train = sort(exc_stimulus_train(:));


			% inhibitory stim soma
			gabarelease_soma = 0;
			gaba_pulse_dur = 2;
			inh_stim_T   = 5; % stimulation period in ms
			inh_n_of_pulses = 3;
			inh_stim_train = cumsum(repmat([0 inhT],1,inh_n_of_pulses));
			
			% inhibitory stim dendrite
			gabarelease_dend = 0;
			gaba_pulse_dur = 1;
			inh_stim_T   = 1; % stimulation period in ms
			inh_n_of_pulses = 2;
			
			% 
			ampa_noise = 0;
			 % ampa_release_prob = .005/(1/delta);

			amplitudes = 2; % TOTAL CURRENT TO THE STIMULATED NEURONS
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

	for g = gaps;

		for nseed = noiseseeds

			for s = firstpulse_inh
				for lag = excitation_lags;
					for excT = exc_stim_T
						for inhT = inh_stim_T
							if gabarelease_dend
						        gaba_mask_dend = perturbation  ;
								gaba_onsets_dend = s + inh_stim_train ;
							else
								gaba_mask_dend = [];
								gaba_onsets_dend = [];
							end

							if gabarelease_soma
						        gaba_mask_soma = perturbation  ; %create_input_mask([5 10 20], 'dist_to_point');
								gaba_onsets_soma = s + inh_stim_train ;
							else
								gaba_mask_soma = [];
								gaba_onsets_soma = [];
							end

							if glurelease
						        glu_mask = perturbation;
								glu_onsets = s + exc_stimulus_train + lag;
							else
								glu_onsets = [];
								glu_mask = [];

							end
							
							condition{conds}.noiseseed = nseed;
							
							condition{conds}.description = 'inhibition and excitation';
							condition{conds}.perturbation_onsets{1} = [];
							condition{conds}.perturbation_mask{1} = [];
							condition{conds}.perturbation_type{1} = [];
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
							condition{conds}.gap = g;

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
			
		end
	end

end



%          __  _                 __      __                               ___ __  _                 
%    _____/ /_(_)___ ___  __  __/ /___ _/ /____     _________  ____  ____/ (_) /_(_)___  ____  _____
%   / ___/ __/ / __ `__ \/ / / / / __ `/ __/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%  (__  ) /_/ / / / / / / /_/ / / /_/ / /_/  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
% /____/\__/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
                                                                                                  


flush = 1;
thisfig = figure;

if stimulate_layer
	for ccc = 1:length(condition);
		pert.mask  	   = condition{ccc}.perturbation_mask;
		pert.amplitude = condition{ccc}.perturbation_amplitude;
		pert.triggers  = condition{ccc}.perturbation_onsets;
		pert.duration  = condition{ccc}.perturbation_pulse_duration;
		pert.type	   = condition{ccc}.perturbation_type;
		
		gap = condition{ccc}.gap;

		noise_level(4) = condition{ccc}.noiseseed;

		sim3D = IOnet('networksize', netsize,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'W',connections.W*gap ,'ou_noise', noise_level , 'perturbation', pert, 'sametoall', sametoall);

		sim3D.condition = condition{ccc};
		sim3D.connections = connections;

		
		if plotresults
			if ~sim3D.failed 
				results{ccc} = replayResults(sim3D, [1:2], 0,[],thisfig);
			    if snapshot
				    export_fig(num2str(ccc)) 
				end
			    orderparameter = measureGlobalSync(sim3D,[1:simtime])
			    sim3D.orderparameter = orderparameter;
			else
				results{ccc} = replayResults(sim3D, [1:2], 0,[],thisfig);
				if snapshot
				    export_fig(num2str(ccc)) 
				end

			    orderparameter.stats.firstordersync = 0;
			    orderparameter.stats.secondordersync = 0;
			    orderparameter.stats.overallsync = 0;
			    sim3D.orderparameter = orderparameter;
			end
		end

		
		if writeresults

		    RESULTS(ccc,1) = mean(mean(W_3d_trans.W(W_3d_trans.W~=0)));
		    RESULTS(ccc,2) = noise_level(1);
		    
		    RESULTS(ccc,5 ) = results{ccc}.popfrequency;
		    RESULTS(ccc,6 ) = results{ccc}.propspkneurons;
		    RESULTS(ccc,7 ) = mean(results{ccc}.spikespercell/transienttime*noneurons);
		    RESULTS(ccc,8 ) = mean(results{ccc}.medfreq(results{ccc}.medfreq>0));
		    RESULTS(ccc,9 ) = mean(sum(W_3d_trans.W))
			
			RESULTS(ccc,13) = orderparameter.stats.firstordersync(1);
			RESULTS(ccc,14) = orderparameter.stats.secondordersync(1);
			RESULTS(ccc,15) = orderparameter.stats.overallsync(1);

		    RESULTS(ccc,16) = sim3D.failed;
		end

		if flush
			save(['sim3D_long_Xcorr_' parameterset '_' num2str(ccc)],  'sim3D')
			clear sim3D
			% gpuDevice(1)
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




