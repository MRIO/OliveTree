% gaba_pspace.m


% ampa_soma_dend_comparison.m


gap = eps;  noisesig = 0; noiseamp = 0 ; tau = 30; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
simtime = 1000;
netsize = [1 1 1];
W.W = 0;
noise_level = [1/tau .005 0 0];
saveappliednoise = 0;
displaytext = ['ampa gaba pspace loyola']
to_report = {'V_soma', 'V_dend'}


M1 = zeros(prod(netsize),1); M1(1) = 1;
M2 = zeros(prod(netsize),1); M2(3) = 1;


 % singlesim

neurons = createDefaultNeurons(1);

dt = 0.01

gaba_conds = linspace(.1,4.1,9);
ampa_conds = linspace(.1,.6,9);


i = 0;
for gcal = [.5 1.1];

	neurons.g_CaL = gcal;
	neurons.g_CaK = 1.5;


	for gbg = gaba_conds
		i = i+1;
		neurons.gbar_gaba_dend = gbg;


		pert.mask{1}  	= M1;
		pert.type{1}	  = 'gaba_dend';
		pert.duration{1}  = 4;
		pert.triggers{1} = [500];

		simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		'cell_parameters',neurons,'W',W.W*gap ,...
		'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		dend_vdend(i,:) = simresults.networkHistory.V_dend;
		dend_vsoma(i,:) = simresults.networkHistory.V_soma;


		

		pert.mask{2}  	  = M1;
		% pert.amplitude{2} = 2;
		pert.type{2}	  = 'gaba_soma';
		pert.duration{2}  = 4;
		pert.triggers{2} = [500];

		simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		'cell_parameters',neurons,'W',W.W*gap ,...
		'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		soma_vdend_gaba(i,:) = simresults.networkHistory.V_dend;
		soma_vsoma_gaba(i,:) = simresults.networkHistory.V_soma;

		[gaba_loc{i} gaba_val{i}] = findpeaks(soma_vdend_gaba(i,:), 'minpeakdistance', 80);

	end



end



i = 0;
for gcal = [.5 1.1];

	neurons.g_CaL = gcal;
	neurons.g_CaK = 1.5;


	for gbg = ampa_conds
		i = i+1;
		neurons.gbar_ampa_soma = gbg;


		pert.mask{1}  	= M1;
		pert.type{1}	  = 'ampa';
		pert.duration{1}  = 1;
		pert.triggers{1} = [550];

		simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		'cell_parameters',neurons,'W',W.W*gap ,...
		'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		vdend_ampa(i,:) = simresults.networkHistory.V_dend;
		vsoma_ampa(i,:) = simresults.networkHistory.V_soma;


		[ampa_loc{i} ampa_val{i}] = findpeaks(vsoma_ampa(i,:), 'minpeakdistance', 80);

	end


end



f = figure;

a(1) = subplot(2,2,1)
plot(soma_vsoma_gaba(1:9,:)', 'linewidth',1)
title('GABA @ dendrite: dampened oscillator')

a(2) = subplot(2,2,2)
plot(soma_vsoma_gaba(10:18,:)', 'linewidth',1)
title('GABA @ dendrite: intrinsic oscillator')
legend(num2str(gaba_conds'))



a(3) = subplot(2,2,3)
plot(vsoma_ampa(1:i/2,:)', 'linewidth',1)
title('AMPA @ soma: dampened oscillator')

a(4) = subplot(2,2,4)
plot(vsoma_ampa(i/2+1:i,:)', 'linewidth',1)
title('AMPA @ soma: intrinsic oscillator')

legend(num2str(ampa_conds'))
f.Color = [1 1 1];;
a(1).ColorOrder = cbrewer('seq', 'Reds', length(gaba_conds));
a(2).ColorOrder = cbrewer('seq', 'Reds',length(gaba_conds));
a(3).ColorOrder = cbrewer('seq', 'Blues', length(ampa_conds));
a(4).ColorOrder = cbrewer('seq', 'Blues', length(ampa_conds));
a(1).XLabel.String = 'ms';
a(2).XLabel.String = 'ms';
a(3).XLabel.String = 'ms';
a(4).XLabel.String = 'ms';
a(1).YLabel.String = 'mV';
a(2).YLabel.String = 'mV';
a(3).YLabel.String = 'mV';
a(4).YLabel.String = 'mV';

a(1).LineWidth = 2
a(2).LineWidth = 2
a(3).LineWidth = 2
a(4).LineWidth = 2


a(1).XLim = [300 600]
a(2).XLim = [300 600]
a(3).XLim = [350 650]





