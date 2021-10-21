% gaba_pspace.m


% ampa_soma_dend_comparison.m


gap = eps;  noisesig = 0; noiseamp = 0 ; tau = 30; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
simtime = 500;
netsize = [1 1 1]

M1 = zeros(prod(netsize),1); M1(1) = 1;
M2 = zeros(prod(netsize),1); M2(3) = 1;


 % singlesim

dt = 0.25

gaba_conds = [.1:.5:4.1];
ampa_conds = [.1:.1:.8];


i = 0;
for gcal = [.5 1.1];

	neurons.g_CaL = gcal;
	neurons.g_CaH = 1.5;


	for gbg = gaba_conds
		i = i+1;
		neurons.gbar_gaba_dend = gbg;


		pert.mask{1}  	= M1;
		pert.type{1}	  = 'gaba_dend';
		pert.duration{1}  = 4;
		pert.triggers{1} = [1000];

		simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		'cell_parameters',neurons,'W',W.W*gap ,...
		'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		dend_vdend(i,:) = simresults.networkHistory.V_dend;
		dend_vsoma(i,:) = simresults.networkHistory.V_soma;


		

		% pert.mask{2}  	  = M1;
		% % pert.amplitude{2} = 2;
		% pert.type{2}	  = 'gaba_soma';
		% pert.duration{2}  = 5;
		% pert.triggers{2} = [1000];

		% simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		% 'cell_parameters',neurons,'W',W.W*gap ,...
		% 'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		% 'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		% % soma_vdend(i,:) = simresults.networkHistory.V_dend;
		% soma_vsoma_gaba(i,:) = simresults.networkHistory.V_soma;

		[gaba_loc(i) gaba_val(i)] = findpeaks(soma_vdend_gaba(i,1000:end), 'minpeakdistance', 80);

	end


end


neurons.g_CaL = 1.5;


i = 0;
for param = [0 2];

	neurons.g_CaCC = param;

	for gba = [.1:.1:.8]
		i = i+1;
		neurons.gbar_ampa_dend = gba;


		pert.mask{1}  	= M1;
		pert.type{1}	  = 'ampa_dend';
		pert.duration{1}  = 2;
		pert.triggers{1} = [200];

		simresults = IOnet('networksize', netsize,'time',simtime,'delta',dt,...
		'cell_parameters',neurons,'W',W.W*gap ,...
		'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,...
		'saveappliednoise',saveappliednoise, 'displaytext',displaytext , 'to_report', to_report);

		% dend_vdend_dendampa(i,:) = simresults.networkHistory.V_dend;
		dend_vsoma_dendampa(i,:) = simresults.networkHistory.V_soma;


		[ampa_loc(i) ampa_val(i)] = findpeaks(dend_vdend_dendampa(i,550:end), 'minpeakdistance', 80);


	end

end




figure
subplot(2,2,3)
title('dampened oscillator')
waterfall(dend_vsoma(1:9,:))
subplot(2,2,4)
title('intrinsic oscillator')
waterfall(dend_vsoma(10:18,:))
legend(num2str(gaba_conds'))



subplot(2,2,1)
title('dampened oscillator')
waterfall(dend_vsoma_dendampa(1:i/2,:))
subplot(2,2,2)
title('intrinsic oscillator')
waterfall(dend_vsoma_dendampa(i/2+1:i,:))
legend(num2str(ampa_conds'))



% ==



f = figure;

a(1) = subplot(2,2,3)
plot(dend_vsoma(1:9,:)', 'linewidth',1)
title('dendGABA @ dampened oscillator')

a(2) = subplot(2,2,4)
plot(dend_vsoma(10:18,:)', 'linewidth',1)
title('GABA @ intrinsic oscillator')
legend(num2str(gaba_conds'))



a(3) = subplot(2,2,1)
plot(dend_vsoma_dendampa(1:i/2,:)', 'linewidth',1)
title('AMPA @ dampened oscillator')

a(4) = subplot(2,2,2)
plot(dend_vsoma_dendampa(i/2+1:i,:)', 'linewidth',1)
title('AMPA @ intrinsic oscillator')

legend(num2str(ampa_conds'))
f.Color = [1 1 1];;
a(1).ColorOrder = cbrewer('seq', 'Blues', 9);
a(2).ColorOrder = cbrewer('seq', 'Purples',9);
a(3).ColorOrder = cbrewer('seq', 'Greens', 9);
a(4).ColorOrder = cbrewer('seq', 'Oranges', 9);
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


a(1).XLim = [500 1500]
a(2).LineWidth = 2
a(3).LineWidth = 2
a(4).LineWidth = 2
