% simulate_sequential_perturbation.m

prep_conditions    = 1;
compute_transients = 1;
stimulate_layer    = 0;
volume_animation   = 0;

savemovies = 1;



rng(1, 'twister')




alpha  = [32.5 32.5]; %edges for alpha parameter
beta   = [12 12];       %edges for beta parameter
res = 1;
scale = 1.35; % 1.35

GCAL = BetaDistributions('alpha', alpha, 'beta', beta, 'no_draws',breadth*depth*height, 'scale', scale,'plot_distributions',false); 
def_neurons.g_CaL = GCAL.sampleDraws{1}';



parameterset = 'pset2';
switch parameterset
	case 'pset1'
		depth = 5; breadth = 10; height = 20;
		networksize = [depth breadth height];
		noneurons = prod(networksize);

		def_neurons = createDefaultNeurons(noneurons)

		offset = [ 0 0 0];

		W_3d = createW('3d',[depth breadth height],0,0.005, 1);
		W = W_3d.W;

		% duration of the pulses
		stim_dur = 10; % du - frequency of stimulation (1000Hz pulses, 100Hz pulse duration)
		stim_T   = 5000; % stimulation period in ms
		simtime  = 2000;
		stim_off = 10000; % t of last perturbation 

		noise_level = 6; % pA per cell
		amplitudes = [-5]; % high values tend to break the explicit euler solver

		% layer 1 - inh on 100ms off 200ms
		% layer 3 - inh on 100ms off 180ms
		% layer 5 - inh on 100ms off 160ms
		% ...

		perturbation{1} = zeros(depth, breadth, height);
		perturbation{2} = zeros(depth, breadth, height);
		perturbation{3} = zeros(depth, breadth, height);
		perturbation{4} = zeros(depth, breadth, height);
		perturbation{5} = zeros(depth, breadth, height);
		perturbation{6} = zeros(depth, breadth, height);
		perturbation{7} = zeros(depth, breadth, height);
		perturbation{8} = zeros(depth, breadth, height);
		perturbation{9} = zeros(depth, breadth, height);
		perturbation{10} = zeros(depth, breadth, height);
	
		perturbation{1}(1:depth, 1:breadth,  1)  = 1;
		perturbation{2}(1:depth, 1:breadth,  3)  = 1;
		perturbation{3}(1:depth, 1:breadth,  5)  = 1;
		perturbation{4}(1:depth, 1:breadth,  7)  = 1;
		perturbation{5}(1:depth, 1:breadth,  9)  = 1;
		perturbation{6}(1:depth, 1:breadth,  11)  = 1;
		perturbation{7}(1:depth, 1:breadth,  13)  = 1;
		perturbation{8}(1:depth, 1:breadth,  15)  = 1;
		perturbation{9}(1:depth, 1:breadth,  17)  = 1;
		perturbation{10}(1:depth, 1:breadth, 19) = 1;

		perturbation{1} = reshape(perturbation{1} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{2} = reshape(perturbation{2} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{3} = reshape(perturbation{3} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{4} = reshape(perturbation{4} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{5} = reshape(perturbation{5} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{6} = reshape(perturbation{6} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{7} = reshape(perturbation{7} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{8} = reshape(perturbation{8} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{9} = reshape(perturbation{9} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{10} =reshape(perturbation{10}, breadth*height*depth,1)*amplitudes(1);

		stiminterv{1} =  [100:200];
		stiminterv{2} =  [100:190];
		stiminterv{3} =  [100:180];
		stiminterv{4} =  [100:170];
		stiminterv{5} =  [100:160];
		stiminterv{6} =  [100:150];
		stiminterv{7} =  [100:140];
		stiminterv{8} =  [100:130];
		stiminterv{9} =  [100:120];
		stiminterv{10} = [100:110];


case 'pset2'
		depth = 5; breadth = 10; height = 20;
		noneurons = breadth*depth*height;

		offset = [ 0 0 0];

		W_3d = createW('3d',[depth breadth height],0,0.005, 1);
		W = W_3d.W;

		% duration of the pulses
		stim_dur = 10; % du - frequency of stimulation (1000Hz pulses, 100Hz pulse duration)
		stim_T   = 5000; % stimulation period in ms
		simtime  = 2000;
		stim_off = 10000; % t of last perturbation 

		noise_level = 6; % pA per cell
		amplitudes = [-5]; % high values tend to break the explicit euler solver

		% layer 1 - inh on 100ms off 200ms
		% layer 3 - inh on 100ms off 180ms
		% layer 5 - inh on 100ms off 160ms
		% ...

		perturbation{1} = zeros(depth, breadth, height);
		perturbation{2} = zeros(depth, breadth, height);
		perturbation{3} = zeros(depth, breadth, height);
		perturbation{4} = zeros(depth, breadth, height);
		perturbation{5} = zeros(depth, breadth, height);
		perturbation{6} = zeros(depth, breadth, height);
		perturbation{7} = zeros(depth, breadth, height);
		perturbation{8} = zeros(depth, breadth, height);
		perturbation{9} = zeros(depth, breadth, height);
		perturbation{10} = zeros(depth, breadth, height);
	
		perturbation{1}(1:depth, 1:breadth,  1)  = 1;
		perturbation{2}(1:depth, 1:breadth,  3)  = 1;
		perturbation{3}(1:depth, 1:breadth,  5)  = 1;
		perturbation{4}(1:depth, 1:breadth,  7)  = 1;
		perturbation{5}(1:depth, 1:breadth,  9)  = 1;
		perturbation{6}(1:depth, 1:breadth,  11)  = 1;
		perturbation{7}(1:depth, 1:breadth,  13)  = 1;
		perturbation{8}(1:depth, 1:breadth,  15)  = 1;
		perturbation{9}(1:depth, 1:breadth,  17)  = 1;
		perturbation{10}(1:depth, 1:breadth, 19) = 1;

		perturbation{1} = reshape(perturbation{1} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{2} = reshape(perturbation{2} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{3} = reshape(perturbation{3} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{4} = reshape(perturbation{4} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{5} = reshape(perturbation{5} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{6} = reshape(perturbation{6} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{7} = reshape(perturbation{7} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{8} = reshape(perturbation{8} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{9} = reshape(perturbation{9} , breadth*height*depth,1)*amplitudes(1); 
		perturbation{10} =reshape(perturbation{10}, breadth*height*depth,1)*amplitudes(1);


		stimdur = 20;
		initstim = 100;
		istimi = 30;
		march  = 0;
		for si = 1:10

			stiminterv{si} =  [initstim+march:initstim+20+march];
			march = march + istimi;
		end


end


dt =0.05;
fs = 1/dt;

	




if prep_conditions

	conds = 1; clear conditions
	noise = randn(noneurons,fs*simtime)*noise_level;

        I_app = zeros(noneurons,fs*simtime) + randn(noneurons,fs*simtime)*noise_level;
        I_app_add = zeros(noneurons,fs*simtime);
		for s = 1:10
			
	        perturbation_triggers = stiminterv{s}*fs;
	        
	        I_app_add(:,perturbation_triggers) = repmat(perturbation{s},[1,length(perturbation_triggers)]);

	        I_app = I_app + I_app_add;

	  
		end
		condition{conds}.perturbation_amplitude = amplitudes; % average current per stimulated neuron in mask
		condition{conds}.perturbation_onsets = 'sequential_perturbations'; % periodicity in ms
		condition{conds}.noise_level = noise_level; %pA per cell
		condition{conds}.offset = [offset 0];
		condition{conds}.gaussalpha = [];
		condition{conds}.simtime = simtime;
		condition{conds}.perturbation_map = [];
		condition{conds}.simtime = simtime;            
        condition{conds}.I_app = I_app;
        condition{conds}.netsize = [depth breadth height]; 

			

	display_applied_current = 1;
	if display_applied_current
		figure
		I = reshape(sum(I_app(:,1:200),2),[depth breadth height]);
		% Iq = interp3(I, 3); Iq = Iq/simtime;
		% h = vol3d('cdata',Iq,'texture','3D');
		h = vol3d('cdata',I,'texture','3D');
		alphamap([0:0.01:.2]);
		view(3)
		axis off
		axis tight;  daspect([1 1 1])

	end

end



if compute_transients
	transienttime = 1000;
	desync_noise = randn(noneurons,fs*transienttime)*2;

	% [transients] = IOtoplevelCa_gpu('rows',breadth,'columns',depth*height,'appCurrent',noise,'time',simtime,'delta',dt,'g_CaL',g_CaL ,'W',W);
	% [transients] = IOtoplevelCa_gpu('rows',depth,'columns',breadth*height,'appCurrent',noise,'time',500,'delta',dt,'g_CaL',g_CaL ,'W',W);
    [transients] = IOnet('rows',depth,'columns',breadth*height,'appCurrent',noise,'time',transienttime,'delta',dt,'g_CaL',g_CaL ,'W',W);
	
	transients.g_CaL = g_CaL;
	transients.condition.simtime = simtime;
	transients.condition.perturbation_map = [];
	transients.condition.W = W;
	transients.condition.distribution_parameters = distribution_parameters;
	
	transients.rows = breadth*depth;
	transients.columns = height;
    

    transients.condition.perturbation_amplitude = 0;
    transients.condition.perturbation_onsets = [];
    transients.condition.noise_level = noise_level;
    transients.condition.offset = [0 0 0 0];
    transients.condition.gaussalpha = 0;
	transients.condition.perturbation_map = [];
    transients.condition.I_app = 0;
    transients.condition.seed = rng(1, 'twister');

end

flush = 1;

if stimulate_layer
	% for ccc = 22;
	for ccc = 1:length(condition);
		I_app = condition{ccc}.I_app;
		sim3D = IOnet('networksize', [depth breadth height], appCurrent',I_app,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'tempState',transients.lastState,'W',W);

        condition{ccc}.I_app = 'notwritinginput';
		sim3D.rows = breadth*depth;
		sim3D.columns = height;
		sim3D.condition = condition{ccc};

		if savemovies
			ons = condition{ccc}.perturbation_onsets;
			amps = condition{ccc}.perturbation_amplitude;

			replayResults(sim3D,[1:1000],1);
			eval(['!mv sim.avi seq_pert.mp4'])

			close all
		end

	if volume_animation

		% which sim to display
		sim = sim3D;
	    animate_volume(sim,[],savemovies)
	    if savemovies
			eval(['!mv volume.mp4 seq_pert_vol.mp4'])
			close all
		end
	end


		if flush
			save(['sim3D_sequential'],  'sim3D')
			clear sim3D
			gpuDevice(1)
		end
	end


end
	
  


% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 




