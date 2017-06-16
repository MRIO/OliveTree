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
phase_partitions  = 13;
gbar_ampa = .25;
reset_pulse = 90;


prep_conditions    = 1;
compute_transients = 1; transienttime = 300;
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
	simtime  = 350;

% network parameters
	netsize = [10 10 2];
	noneurons = prod(netsize);

	
	gap = .04; % mS/cm^2
	% gap = eps; % mS/cm^2
	
	
	% [=================================================================]
	%  % create net and neurons
	% [=================================================================]

	rng(0,'twister')

	plotthis  = 0;
	rd = 2;
	meannoconn = 8;

	normleak  = 1;
	randomize = 1;
	scaling   = 1;
	maxiter	  = 1;
	somatapositions = [];
	randomize = 1;
	symmetrize = 1;

	W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
		W.W = W.W*gap;

	def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized');
		def_neurons.gbar_ampa_soma = gbar_ampa*ones(noneurons,1);
		% def_neurons.g_CaL = 1.4*ones(noneurons,1);
		

	

	% W_3d = createW('3d', netsize,2.5,gap, 1, 0,5);
	% W_3d.W = 0;


% [=================================================================]
%  % noise levels
% [=================================================================]

 tau = 20  		;
 noisemu = -0.6	;
 noisesig = 0.6	;

% noise_level = [1/tau noisesig noisemu 0];
noise_level = [0 0 0 0];
noise_parameters = [0 0 0 0];

% [=================================================================]
%  perturbations
% [=================================================================]
inpurrad = 3;
% perturbation
	pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', inpurrad,'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',plotthis);
	% pert.mask  	  {1} = create_input_mask(netsize, 'dist_to_center','radius',3, 'synapseprobability', 1,'plotme',1);
	pert.amplitude{1} = 1;
	pert.duration {1} = 1;
	pert.type	  {1} = synapse_type;
	pert.triggers {1} = [reset_pulse];

% noise
	sametoall = 0.0;

% [=================================================================]
% calc steady state
% [=================================================================]

if ~exist('steady_state')
	[steady_state] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',0.02,'cell_parameters', def_neurons , 'cell_function', cell_function, 'W',W.W,'ou_noise', noise_parameters, 'sametoall',sametoall);
	% replayResults_3(steady_state)
	% drawnow
end

if ~exist('unperturbed_state')
	
	[unperturbed_state] = IOnet( 'networksize', netsize ,'time',simtime, 'perturbation', pert, 'cell_function', cell_function, ...
								'delta',dt,'cell_parameters', def_neurons ,'W',W.W,'ou_noise', noise_parameters, ...
	 							'sametoall',sametoall, 'tempState', steady_state.lastState);
	% replayResults_3(unperturbed_state)
	% drawnow

end


% measure unperturbed phases
phaseResults = measureGlobalSync(unperturbed_state, 'duration', [1:simtime],'plotme', 1);
phases = phaseResults.hilbert.hilbert;

% find original peaks
[PKS LOCS] = findpeaks(mean(phases,1), 'minpeakdistance',50);
LOCS(LOCS<reset_pulse)  = [];
meanT = LOCS(2) - LOCS(1);

% return


estimatedT = 240;

pertphases = round(linspace(reset_pulse,reset_pulse+estimatedT,phase_partitions));


% add perturbation at different phases and compute
def_neurons.gbar_ampa_soma = gbar_ampa*ones(noneurons,1);
VV = unperturbed_state.networkHistory.V_soma;
k = 0;
for pertphase = pertphases
	k = k+1;

	pert.triggers {1} = [ reset_pulse pertphase];

	simPert{k} = IOnet( 'networksize', netsize , 'time',simtime,'delta',dt, ...
		'cell_parameters', def_neurons , 'cell_function', cell_function, 'W',W.W,'ou_noise', noise_parameters,...
		'sametoall',sametoall, 'perturbation', pert, 'tempState', steady_state.lastState);

	% replayResults_3(simPert{k})
	% drawnow
	
	simPert{k}.condition.perturbation_mask = pert.mask{1};
	simPert{k}.condition.perturbation_onsets = pert.triggers{1};

	VV = [VV ; simPert{k}.networkHistory.V_soma];

end



% compute phases and check
k = 0; t = [LOCS(2)-60:LOCS(2)+60]; PP = mean(phases);
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

%============================= plots ==============================%




% membrane potential of stimulated cells
no_stimcells = length(find(pert.mask{1}));
VVV = VV(find(repmat(pert.mask{1},1, length(pertphases))),:);

% phase of stimulated cells
PPP = PP(find(repmat(pert.mask{1},1, length(pertphases))),:);

% spike count of membrane potential of stimulated cells
SSS = sum(VVV(:,150:350)>-30,2)>0;

for c = 1:length(pertphases)
	spkcells(c) = sum(SSS(no_stimcells*(c-1)+1:no_stimcells*c));
end


if ~exist('cbrewer')
		lc = jet(length(pertphases));
	else
		lc = cbrewer('qual', 'Set1', length(pertphases))
end


%============================= membrane potential of stimulated groups==============================%
figure
	waterfall(VVV(no_stimcells+1:end,:)), zlabel('mV'), ylabel('perturbation group'), xlabel('ms')
	set(gca,'ytick', [], 'xtick', [90 190], 'xticklabel', [0 100] )
	colormap(bone)
	hold on
	
	for c = 1:length(pertphases)
		line([simtime simtime] , [no_stimcells*(c-1)+1 no_stimcells*c], [-70 -70], 'color', lc(c,:), 'linewidth',5);
	end	


figure
	subplot(2,1,1)
	imagesc(VVV(no_stimcells+1:end,:)), ylabel('perturbation group'), xlabel('ms')
	set(gca,'ytick', [],'clim',[-80 -30])
	colormap(bone)
	hold on	
	for c = 1:length(pertphases)
		line([simtime simtime] , [no_stimcells*(c-1)+1 no_stimcells*c],  'color', lc(c,:), 'linewidth',5);
	end
	%=============================unperturbed membrane potential==============================%
	subplot(2,1,2)
	imagesc(VVV(1:no_stimcells,:)), ylabel('perturbation group'), xlabel('ms')
	set(gca,'ytick', [],'clim',[-80 -30])
	

%=============================spikes per cell histogram==============================%
figure
	 bar((pertphases-pertphases(1))/meanT*2*pi,spkcells/no_stimcells)
	 ylabel('probability of CS')
	 xlabel('stimulation phase')
	hold on
	set(gca,'xtick', [0:.5*pi:4*pi])

%=============================PRC ==============================%
%  only valid for the first period
figure
	subplot(3,1,1)
	plot((pertphases-pertphases(1))/meanT*2*pi, PRC)
	axis tight
	ylabel('perturbation phase (radians)')
	xlabel('\Delta phase')
	title('PRC')
	set(gca,'xtick', [0:.5*pi:4*pi])
	%============================= means ==============================%
	subplot(3,1,2)
	plot_mean_and_std(VVV(1:54,:),[1:350])
	xlim([pertphases(1) pertphases(end)])
	xlabel('ms')
	ylabel('mV')

	subplot(3,1,3)
	for c = 1:length(pertphases)
		plot_mean_and_std( VVV(1+(c-1)*no_stimcells:c*no_stimcells,:),[1:350],'color', lc(c,:))
		hold on
	end
	xlabel('time (ms)')
	ylabel('mV')




figure
	for c = 1:length(pertphases)
		line(mean(PPP(1+(c-1)*no_stimcells:c*no_stimcells,:)),[1:350],'color', lc(c,:))
		hold on
	end
	xlabel('time (ms)')
	ylabel('mV')

	xlim([pertphases(1) pertphases(end)])
	xlabel('ms')
	ylabel('phase')




figure
	for c = 1:length(pertphases)
		line([1:350], mean(VVV(1+(c-1)*no_stimcells:c*no_stimcells,:)),'color', lc(c,:))
		hold on
	end
	xlabel('time (ms)')
	ylabel('mV')

	xlim([pertphases(1) pertphases(end)])
	set(gca, 'xtick', [90 190], 'xticklabel', [0 100] )









