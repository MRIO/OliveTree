% generate_PRC.m

% 0. generate steady_state
% 1.   Generate unperturbed network
% 1.2. Measure phase
% 1.3. Partition one phase in 9 intervals from 0 to 1

% 1.4. Stimulate network 9 times
% 1.5. Measure phases of stimulated population
% 1.6. For every stimulated phase, find next peak

% Q: how does the prc correlate with different ion channel properties?

% assumptions:    
% 
% 	1. Cells have slightly different oscillatory properties as a
% 		function of their randomized calcium distribution 
%   2. Cells have slightly different after depolarization properties owing to varying I_h and I_int
% 
%   3. 100 cells in a deformed 2d mesh, connected with a 2 cell radius (#decisionpoint)
% 
%   4. add some additive noise, uncorrelated between cells
% 
% 
% initial conditions:
% 		All simulations have the same initial state:
%		   Syncrhonized transients in the absence of noise due to 1 and 2 and no gap junctions.
% 
% manipulation:
% 
% 		Ampa or GABA pulse to subpopulation at increasing phases 
% 
% experimental conditions:
% 
%	1. With gaps and without gaps
%	2. With input and without input
% 
% measurements:
%   1. Phase difference with respect to non-stimulated network
% 	2. Time to for whole network to resync
% 

% [=================================================================]
%  script parameters
% [=================================================================]

% cell_function = 'devel';
cell_function = 'vanilla';
% synapse_type = 'gaba_soma';
synapse_type = 'ampa';
phase_partitions  = 15;

prep_conditions    = 1;
compute_transients = 1; transienttime = 500;
stimulate    = 1;
volumetric_activity = 0; % slow

savemovies = 0;
snapshot = 1;
gpu = 0;

% rseed = 5; % this seed leads to wave propagations and solitons
rseed = 1;
rng(rseed, 'twister')
to_report = {'V_soma'};

% global 
	dt =0.02;
	fs = 1/dt;
	simtime  = 450;

% network parameters
	depth = 2; breadth = 1; height = 1;
	netsize = [depth breadth height];
	noneurons = prod(netsize);
	
	gap = eps; % mS/cm^2
	
	
	def_neurons = createDefaultNeurons(noneurons);
	% def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized');
		% def_neurons.gbar_gaba_soma = ones(noneurons,1);
		% def_neurons.gbar_gaba_dend = ones(noneurons,1);
		def_neurons.gbar_ampa_soma = ones(noneurons,1)*.1;
		def_neurons.g_CaL = ones(noneurons,1)*1.1;

	% W_3d = createW('3d', netsize,2.5,gap, 1, 0,5);
	W_3d.W = 0;

	noise_parameters = [0 0 0 rseed]; % pA/ms per cell - 3.5 3 1

% perturbation
	pert.mask  	  {1} = create_input_mask(netsize, 'dist_to_center','radius',3, 'synapseprobability', 1,'plotme',1);
	pert.amplitude{1} = 1;
	pert.duration {1} = 1;
	pert.type	  {1} = synapse_type;
	pert.triggers{1} = [100];

% noise
	sametoall = 0.0;

% calc steady state
if ~exist('steady_state')
	[steady_state] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',0.02,'cell_parameters', def_neurons , 'cell_function', cell_function, 'W',W_3d.W,'ou_noise', noise_parameters, 'sametoall',sametoall);
end

if ~exist('unperturbed_state')
	
	[unperturbed_state] = IOnet( 'networksize', netsize ,'time',simtime, 'perturbation', pert, 'cell_function', cell_function, ...
								'delta',dt,'cell_parameters', def_neurons ,'W',W_3d.W,'ou_noise', noise_parameters, ...
	 							'sametoall',sametoall, 'tempState', steady_state.lastState);
end

% measure unperturbed phases
phaseResults = measureGlobalSync(unperturbed_state, 'duration', [1:simtime],'plotme', 1);
phases = phaseResults.hilbert.hilbert;

% find original peaks
[PKS LOCS] = findpeaks(mean(phases,1), 'minpeakdistance',50);
LOCS(1) = [];
meanT = LOCS(2) - LOCS(1);

% return

pertphases = round(linspace(LOCS(1),LOCS(2),phase_partitions));


% add perturbation at different phases and compute
VV = unperturbed_state.networkHistory.V_soma;
k = 0;
for pertphase = pertphases
	k = k+1;

	pert.triggers {1} = pertphase;

	simPert{k} = IOnet( 'networksize', netsize , 'time',simtime,'delta',dt, ...
		'cell_parameters', def_neurons , 'cell_function', cell_function, 'W',W_3d.W,'ou_noise', noise_parameters,...
		'sametoall',sametoall, 'perturbation', pert, 'tempState', steady_state.lastState);

	replayResults_3(simPert{k})
	drawnow
	
	simPert{k}.condition.perturbation_mask = pert.mask{1};
	simPert{k}.condition.perturbation_onsets = pert.triggers{1};

	VV = [VV ; simPert{k}.networkHistory.V_soma];

end



% compute phases and check
k = 0; t = [LOCS(2)-50:LOCS(2)+50]; PP = mean(phases);
for pertphase = pertphases
	k = k+1;

	perturberdPhaseResults = measureGlobalSync(simPert{k},'duration', [1:simtime], 'plotme', 0);
		perturbedPhases{k} = perturberdPhaseResults.hilbert.hilbert;

		PP = [PP; perturbedPhases{k}];

	% find new peaks
	[newPKS newLOCS{k}] = findpeaks(mean(perturbedPhases{k}(find(pert.mask{1}) , t)), 'minpeakdistance',40);

	peakDelta(k) = LOCS(2) - (t(1) + newLOCS{k}(1));
	PRC(k) = peakDelta(k)/meanT;


end

PRC = peakDelta/meanT;


% k = 0; t = [LOCS(2)-50:LOCS(2)+50];
% for pertphase = pertphases
% 	k = k+1;

% 	perturberdPhaseResults = measureGlobalSync(simPert{k}, 'duration', [1:simtime],'plotme', 1);
% 		perturbedPhases{k} = perturberdPhaseResults.hilbert.hilbert;
		
% 		sync_stim(k) = perturberdPhaseResults.stats.firstordersync(1);
% 		sync_all(k)  = perturberdPhaseResults.stats.overallsync(1);


% 	% find new peaks
% 	meanV_stim(k,:)   = mean(simPert{k}.networkHistory.V_soma(find(pert.mask{1}) , :));
% 	meanV_nostim(k,:) = mean(simPert{k}.networkHistory.V_soma(find(~pert.mask{1}) , :));


% end


