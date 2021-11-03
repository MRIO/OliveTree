function out =  IO_PRC(neuron, pert, phase_partitions)

% given a single InfOli neuron, generate it's PRC for the specified synapse type.


% 0. generate steady_state
% 1.   Generate unperturbed network
% 1.2. Measure phase
% 1.3. Partition one phase in N intervals from 0 to 1

% 1.4. Stimulate network N times
% 1.5. Measure phases of stimulated population
% 1.6. For every stimulated phase, find next peak




cell_function = 'vanilla';


simtime  = 800;
compute_transients = 1; transienttime = 500;
simtime = 700
stimulate    = 1;
volumetric_activity = 0; % slow

savemovies = 0;
snapshot = 1;
gpu = 0;

rseed = 1;
rng(rseed, 'twister')
to_report = {'V_soma','V_dend'};

% global 
	dt =0.02;
	fs = 1/dt;
	
% network parameters
	gap = eps; % mS/cm^2
	
% network

    netsize = [1 1 1];

% perturbation parameters

cell_function = 'vanilla';


% noise
	sametoall = 0.0;
	noise_parameters = [0 0 0 rseed]; % pA/ms per cell - 3.5 3 1


% calc steady state

    description = 'steady_state';
% 	  Vclamp = -60; [steady_state] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',0.02,'cell_parameters', neuron , 'cell_function', cell_function, 'W',0,'ou_noise', noise_parameters, 'sametoall',sametoall, 'appVoltage', Vclamp, 'displaytext', description);
    [steady_state] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',dt,'cell_parameters', neuron , 'cell_function', cell_function, 'W',0,'ou_noise', noise_parameters, 'sametoall',sametoall, 'displaytext', description);


	description = 'unperturbed_state';
	[unperturbed_state] = IOnet( 'networksize', netsize ,'time',simtime, 'perturbation', [], 'cell_function', cell_function, ...
								'delta',dt,'cell_parameters', neuron,'W',0,'ou_noise', noise_parameters, ...
	 							'sametoall',sametoall, 'tempState', steady_state.lastState, 'displaytext', description);

    vs_amplitude = max(unperturbed_state.networkHistory.V_soma) - min(unperturbed_state.networkHistory.V_soma);



% measure unperturbed phases
phaseResults = measureGlobalSync(unperturbed_state, 'duration', [1:simtime],'plotme', 0);
phases = phaseResults.hilbert.hilbert;

% find original peaks
[PKS LOCS] = findpeaks(phases, 'minpeakdistance',50, 'minpeakheight', 5);
LOCS(1) = []; % use second period to perturb (safer)


if isempty(LOCS)
    disp('IO_PRC found no peaks')
    keyboard
else
    periods=diff(LOCS);
    meanT = mean(periods);
end


pertphases = round(linspace(LOCS(1),LOCS(2),phase_partitions));



%% add perturbation at different phases and compute

VV = unperturbed_state.networkHistory.V_soma;

parfor k = 1:length(pertphases)

    pertphase = pertphases(k);
   
    VS{k} = perturb_net( pertphase, netsize , simtime,dt, ...
		neuron , cell_function, 0 ,noise_parameters,...
        		sametoall, pert,  steady_state);

	VV = [VV ; VS{k}.networkHistory.V_soma];


end

%% compute phases and check
k = 0;  PP = phases;

% Only search for peaks after this interval (ms).  Factor to account for causality
% violation in estimation of phase.
pec = 10;

for pertphase = pertphases
	k = k+1;

        perturbedPhaseResults = measureGlobalSync(VS{k},'duration', [1:simtime], 'plotme', 0);
		perturbedPhases{k} = perturbedPhaseResults.hilbert.hilbert;

		PP = [PP; perturbedPhases{k}];

	% find new peaks

        
        t = [pertphase+pec:LOCS(3)+meanT];

        try
            [newPKS{k} newLOCS{k}] = findpeaks(perturbedPhases{k}(t), 'minpeakdistance',40,  'npeaks', 1, 'minpeakheight', 5.5);
        catch
            newPKS{k} = NaN;
            newLOCS{k} = [];
        end

        if isempty(newLOCS{k})
            peakDelta(k) = NaN;
            PRC(k) = NaN;
        else
            peakDelta(k) = LOCS(2) - (t(1) + newLOCS{k}(1));
            PRC(k) = peakDelta(k)/meanT;
        end


            
        
end
%%
PRC = peakDelta/meanT;

out.PRC = [(pertphases - LOCS(1))/meanT*2*pi ; PRC];
out.peaktimes = LOCS;
out.peakamps  = newPKS;
out.newPeaks  = newLOCS;
out.pertphases = pertphases;
out.VS = VV;
out.steadystate = steady_state;
out.neuron = neuron;

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

end

%% workaround parfor limitation with structures

function VS = perturb_net(pertphase, netsize ,simtime, dt, ...
		 neuron , cell_function, W , noise_parameters,...
		sametoall, pert , steady_state)

    pert.triggers{1} = pertphase;
    
    if length(pert.triggers) == 1
        pert.triggers{2} = [];
    else
        pert.triggers{2} = pert.triggers{2} + pertphase;
    end
    
	VS = IOnet( 'networksize', netsize , 'time',simtime,'delta',dt, ...
		'cell_parameters', neuron , 'cell_function', cell_function, 'W',0 ,'ou_noise', noise_parameters,...
		'sametoall',sametoall, 'perturbation', pert, 'tempState', steady_state.lastState);
    
    
end
    