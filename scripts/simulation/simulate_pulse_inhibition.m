% simulate_pulse_inhibition.m

prep_conditions    = 1;
compute_transients = 1;
stimulate_layer    = 1;
volumetric_activity = 0;

savemovies = 0;
snapshot = 1;


rng(1, 'twister')

% global 
dt =0.025;
fs = 1/dt;

gap = .05;

sametoall = .2; % unless overwritten

parameterset = 'pset1';
switch parameterset
	case 'pset1'
		inputmasktype = 'mid_layer';
		depth = 5; breadth = 10; height = 20;
		networksize = [depth breadth height];
		noneurons = prod(networksize);
		def_neurons = createDefaultNeurons(noneurons);
		% CaT distribution parameters
			alp = 7;
			bet = 3;
			sup = [0.7 1.2];
			GCAL = BetaDistributions('alpha', alp, 'beta', bet, 'no_draws', ...
		    depth*height*breadth, 'support', sup,'plot_distributions',0); 
	
			def_neurons.g_CaL = GCAL.sampleDraws{1}';
			def_neurons.g_int = [ones(1,500)*0.013 ; ones(1,500)*0.01];

		gap = 0.03;
		W_3d_trans = createW('3d_euclidean_rndwalk',[depth breadth height],2.5,gap, 1,0);


		noise_level = [3 8 0 0]; % pA/ms per cell
		
		firstpulse_inh = [300]; % start of stimulation
		
		% inhibitory step
		gabarelease = 1;
		inh_pulse_dur = [5]; % time to peak of [GABA] concetration
		inh_stim_T   = [5 10 20 50 100]; % stimulation period in ms
		inh_n_of_pulses = 10;
		
		glurelease = 0;
		ampanoise  = 0;		
		

		amplitudes = 1; % TOTAL CURRENT TO THE STIMULATED NEURONS
		% stimulus times
		
	
end	


%        __     ____            ____                __     _                  
%   ____/ /__  / __/___ ___  __/ / /_   _________ _/ /____(_)_  ______ ___    
%  / __  / _ \/ /_/ __ `/ / / / / __/  / ___/ __ `/ / ___/ / / / / __ `__ \   
% / /_/ /  __/ __/ /_/ / /_/ / / /_   / /__/ /_/ / / /__/ / /_/ / / / / / /  /
% \__,_/\___/_/  \__,_/\__,_/_/\__/   \___/\__,_/_/\___/_/\__,_/_/ /_/ /_/   \
                                                                              


% create masks for input
gaba_mask = create_input_mask(networksize, inputmasktype);

%                                                                    ___ __  _                 
%     ____  ________  ____  ____ _________     _________  ____  ____/ (_) /_(_)___  ____  _____
%    / __ \/ ___/ _ \/ __ \/ __ `/ ___/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%   / /_/ / /  /  __/ /_/ / /_/ / /  /  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
%  / .___/_/   \___/ .___/\__,_/_/   \___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
% /_/             /_/                                                                          

if prep_conditions

	conds = 1; clear conditions

		for iT = inh_stim_T


				if gabarelease
			        gaba_mask = gaba_mask;
					gaba_onsets = firstpulse_inh - iT + cumsum(repmat(iT , 1 , inh_n_of_pulses)) ;
				else
					gaba_mask = [];
					gaba_onsets = [];
				end

				
				if ampanoise
					ampanoise_mask = [1:noneurons];
					ampanoise_lambda = 10;
				end

					condition{conds}.perturbation_onsets = gaba_onsets;
					condition{conds}.perturbation_mask 	= gaba_mask; % selection of stimulated neurons
					
					condition{conds}.perturbation_amplitude = 1;

		            conds = conds+1;
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
	transienttime = 50;
	
    [transients] = IOnet( 'networksize',[depth breadth height],'time',transienttime,'delta',0.05,'cell_parameters', def_neurons ,'W',W_3d_trans.W,'ou_noise', [noise_level 0 1 0]);
end


%          __  _                 __      __                               ___ __  _                 
%    _____/ /_(_)___ ___  __  __/ /___ _/ /____     _________  ____  ____/ (_) /_(_)___  ____  _____
%   / ___/ __/ / __ `__ \/ / / / / __ `/ __/ _ \   / ___/ __ \/ __ \/ __  / / __/ / __ \/ __ \/ ___/
%  (__  ) /_/ / / / / / / /_/ / / /_/ / /_/  __/  / /__/ /_/ / / / / /_/ / / /_/ / /_/ / / / (__  ) 
% /____/\__/_/_/ /_/ /_/\__,_/_/\__,_/\__/\___/   \___/\____/_/ /_/\__,_/_/\__/_/\____/_/ /_/____/  
                                                                                                  


flush = 1;

if stimulate_layer
	
	for ccc = 1:length(condition);
		pert.mask{1}  	  = condition{ccc}.perturbation_mask;
		pert.amplitude{1} = condition{ccc}.perturbation_amplitude;
		pert.triggers{1}  = condition{ccc}.perturbation_onsets;
		pert.duration{1}  = 1;
		pert.type{1}		= 'gaba_soma';

		simtime = pert.triggers{1}(end) + 300;

		sim3D =  IOnet('networksize', networksize,'time',simtime,'delta',dt,'cell_parameters',def_neurons, 'tempState',transients.lastState,'W',W_3d_trans.W ,'ou_noise',noise_level , 'perturbation', pert, 'sametoall', sametoall);
		

        condition{ccc}.I_app = 'notwritinginput';
		
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
			results{ccc} = replayResults(sim3D)
		    figure(1)
		    export_fig(num2str(ccc)) 
		    clf
		end


		if flush

			save(['sim3D_gaba_' num2str(ccc)],  'sim3D')
			clear sim3D
			gpuDevice(1)
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




