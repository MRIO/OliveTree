function results = measureGlobalSync(varargin)
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
try
	warning off
 	H = hilbert_of_membranepotential(VS); 
 	warning on

catch E


  warning('could not compute hilbert')
  keyboard

  results.stats.firstordersync = 0;
	results.stats.secondordersync = 0;
	results.stats.overallsync =     0;
	results.stats.order_parameter_g_fo = 0;
	results.stats.order_parameter_g_so = 0;
	results.hilbert = 0;


  return 

end

 
U = H.hilbert;



% [================================================]
%  compute instantaneous freuquencies (per neuron)
% [================================================]

% assuming sampling rate of 1000Hz!

cyclecount = unwrap(U'); 
cyclecount = cyclecount(end,:);
Freq = cyclecount /(2*pi) / length(tt);

InstFreq = diff(unwrap(U)') /(2*pi) * 1e3; % 

[HistInstFreq x] = hist(InstFreq',[-2000:2000]); % per cell
[HistInstFreq_overall x] = hist(InstFreq(:),x); %InstFreq


% [=================================================================]
%  partition
% [=================================================================]


allneu = 1:noneurons;
if not(isempty(group))
	allneu = setdiff(allneu,find(group));
end



GA = U(allneu,:);
MA = circ_mean(GA)+pi;
VA = circ_var(GA);

GR = U(group,:);
MR = circ_mean(GR)+pi;
VR = circ_var(GR);



% KURAMOTO ORDER PARAMETER
% mean of exp( e^(i*(theta_k,p(t) - theta_syn(t)) ; p is group, k is neuron

order_parameter_GA = mean( exp(i*(bsxfun(@minus, GA, MA))));
order_parameter_GR = mean( exp(i*(bsxfun(@minus, GR, MR))));

order_parameter_g_fo = zeros(noneurons, length(tt));
order_parameter_g_so = zeros(noneurons, length(tt));
for neuron = 1:noneurons
	Ind = zeros(noneurons,1);
	Ind(neuron) = 1;

	fo_neighbors = W * Ind;
	so_neighbors = W * fo_neighbors;
	% to_neighbors = W * so_neighbors;


	fo_neighbors_ind = [neuron ; find(W * Ind)];
	so_neighbors_ind = unique([neuron ; fo_neighbors_ind ; find(so_neighbors)]);
	% to_neighbors_ind = unique([neuron ; fo_neighbors_ind ; so_neighbors_ind ; find(to_neighbors)] );


	mean_fo = circ_mean(U(fo_neighbors_ind,:));
	mean_so = circ_mean(U(so_neighbors_ind,:));
	% mean_to = circ_mean(U(to_neighbors_ind,:));
	
try

 
	order_parameter_g_fo(neuron,:) = mean( exp(i*(bsxfun(@minus, U(fo_neighbors_ind,:), mean_fo))));
	order_parameter_g_so(neuron,:) = mean( exp(i*(bsxfun(@minus, U(so_neighbors_ind,:), mean_so))));
	 

catch E

  keyboard

end

	% order_parameter_g_fo(neuron,:) = mean( exp(i*(bsxfun(@minus, U(fo_neighbors_ind,:), mean_fo)+ pi)));
	% order_parameter_g_fo(neuron,:) = mean( exp(i*(bsxfun(@minus, U(fo_neighbors_ind,:), mean_fo)+ pi)));
	% order_parameter_g_so(neuron,:) = mean( exp(i*(bsxfun(@minus, U(so_neighbors_ind,:), mean_so)+ pi)));
	% order_parameter_g_to(neuron,:) = mean( exp(i*(bsxfun(@minus, U(to_neighbors_ind,:), mean_to))));

end

% keyboard

results.stats.firstordersync =  [mean(mean(abs(order_parameter_g_fo))) ; var(mean(abs(order_parameter_g_fo))) ];
results.stats.secondordersync = [mean(mean(abs(order_parameter_g_so))) ; var(mean(abs(order_parameter_g_so))) ];
results.stats.overallsync =     [mean(mean(abs(order_parameter_GA)))   ; var(mean(abs(order_parameter_GA)))   ];
results.stats.order_parameter_all = abs(order_parameter_GA)';
results.stats.order_parameter_group = abs(order_parameter_GR)';
results.stats.order_parameter_group_fo = abs(order_parameter_g_fo)';
results.stats.order_parameter_group_so = abs(order_parameter_g_so)';
results.hilbert = H;
results.instantaneousFrequency = InstFreq;
results.frequency = Freq;
results.spectrum.overall = {HistInstFreq_overall x};
results.spectrum.percell = {HistInstFreq x};

if plotme

	figure
		aaxx(1) = subplot(2,1,1)
		plot(abs(order_parameter_GA),'color', [1 1 1]*.5, 'linewidth',2)
		hold on
		plot(abs(order_parameter_GR),'r', 'linewidth',2)
		hold on
		plot(mean(abs(order_parameter_g_fo)),'r','linestyle','--')
		hold on
		% plot(mean(abs(order_parameter_g_so)),'g')
		% hold on

		% plot(mean(abs(order_parameter_g_to)),'b')
		hold off
		xlabel('ms')
		ylabel('kuramoto parameter')
		legend({'all' 'group' 'first order neighbors' })
		axis tight


		ylim([0 1])

		
		aaxx(2) = subplot(2,1,2)
		% [HHH XXX] = hist(U,100);
		imagesc(hist(U,100), [0 noneurons/10])
		xlabel('ms')
		ylabel('phase')
		title('phase concentration')
		colorbar

		linkaxes(aaxx,'x')

	figure
		imagesc([1:t], x, HistInstFreq , [0 noneurons]/10)
		ylim([-30 30])
		ylabel('Hz')
		xlabel('ms')
		colorbar
		title('instantaneous frequency distribution')

end
% keyboard







