% measureGroupSync.m

function results = measureGroupSync(varargin)
% measureGlobalSync.m
% presumes sample rate 1000Hz
%
% given a simulation, computes, using the DAMOCO toolbox:
% 	the hilber transform (damoco)
% 	the instantaneous frequency of all neurons
% 	the order parameter of all neurons and that of a potential group
% 
% 
% Parameter: ('duration', []) 
% Parameter: ('plotme', 1) 
% Parameter: ('trigger', []) 
% Parameter: ('group',[])


	p = inputParser;
	p.addRequired('sim')  
	
	p.addParameter('duration', []) 
	p.addParameter('plotme', 1) 
	p.addParameter('trigger', []) 
	p.addParameter('group',[])
	
	p.parse(varargin{:});

	sim = p.Results.sim;
	duration = p.Results.duration;
	plotme = p.Results.plotme;
	trigger = p.Results.trigger;
	group = p.Results.group;
	



VS = sim.networkHistory.V_soma;
t = size(VS,2);

noneurons = prod(sim.networksize);

if ~isempty(duration)
	tt = [duration(1):duration(end)];
else
	tt = 1:t;
end

VS = VS(:,tt);
t = size(VS,2);



% retrieve fields
if isempty(trigger)
	if isfield('mask', sim.perturbation) 
			pert_mask			= sum(cell2mat(sim.perturbation.mask),2)>0;
			perturbation_onsets = reshape(cell2mat(sim.perturbation.triggers),1,[]);
	else
		pert_mask = [];
		perturbation_onsets = [];
	end
else
	pert_mask			= sim.perturbation.mask{trigger};
	perturbation_onsets = sim.perturbation.triggers{trigger};	
end	


W = sim.networkParameters.connectivityMatrix;



if isempty(group)
	group = randi(noneurons,20);
end


% [================================================]
% 		 compute hilbert transform (DAMOCO)
% [================================================]
	warning off
 	H = hilbert_of_membranepotential(VS(group,:)  ); 
 	warning on


 
U = H.hilbert;


% [================================================]
%  compute instantaneous freuquencies (per neuron)
% [================================================]

% assuming sampling rate of 1000Hz!

cyclecount = unwrap(U'); 
cyclecount = cyclecount(end,:);
Freq = cyclecount /(2*pi) / length(tt);

InstFreq = diff(unwrap(U)') /(2*pi) * 1e3; % #check


[HistInstFreq x] = hist(InstFreq',[-2000:2000]); % per cell
[HistInstFreq_overall x] = hist(InstFreq(:),x); %InstFreq


% [=================================================================]
%  partition
% [=================================================================]

GA = U;
MA = circ_mean(GA)+pi;
VA = circ_var(GA);


% KURAMOTO ORDER PARAMETER
% mean of exp( e^(i*(theta_k,p(t) - theta_syn(t)) ; p is group, k is neuron

order_parameter_GA = mean( exp(i*(bsxfun(@minus, GA, MA))));


% keyboard

results.sync =     [mean(mean(abs(order_parameter_GA)))   ; var(abs(order_parameter_GA))   ];
results.order_parameter = abs(order_parameter_GA)';
results.hilbert = H;
results.instantaneousFrequency = InstFreq;
results.frequency = Freq;
results.spectrum.overall = {HistInstFreq_overall x};
results.spectrum.percell = {HistInstFreq x};

if plotme

	figure
	subplot(2,1,1)
	plot(abs(order_parameter_GA),'y', 'linewidth',2)
	hold on

	% plot(mean(abs(order_parameter_g_to)),'b')
	hold off
	xlabel('ms')
	ylabel('kuramoto parameter')
	legend({'all' 'group' 'first order neighbors' 'second order'})
	axis tight


	ylim([0 1])

	
	subplot(2,1,2)
	imagesc(hist(U,100))
	xlabel('ms')
	ylabel('phase')

	figure
	imagesc([1:t], x, HistInstFreq)
	ylim([-30 30])
	ylabel('Hz')
	xlabel('ms')

	figure
	plot([1:t], VS(group,:))


end
% keyboard



