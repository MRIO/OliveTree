% Warnaar et al. compare oscillation properties between gap and no gaps

clear;

onesec_vs_30s = 0;
plotcellscatters_gap_gapless  = 0; % compare cell scatters in nets
joinedcellscatter = 0;
triggeredphase =0;
hist_sto_freq_amp_w_wo_gaps = 0;

frequency_drift_singlesim_prc = 0;

calculatesynchrony = 0;
sto_and_propfiring_histograms = 0;

interperiodintervals = 1;

sortedsampletraces = 1;


tslice = 1:50000;

if 0
	addpath('/Users/M/Projects/Experiments/Olive/model/simresults/periodic_ampa/')
	addpath('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/periodic_ampa/')

	F1 = 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017';
	
	F1 = 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_spont_50000_4_17-Jan-2017';

	load(F1)
	disp('loaded.')

	runs_with   = [5:8];
	runs_without = [1:4];

	Joinedsim{1}  = joinsim(simresults,runs_without); 
	Joinedsim{2}  = joinsim(simresults,runs_with); 


	statevar{1} = Joinedsim{1}.networkHistory.V_soma(:,tslice);
	statevar{2} = Joinedsim{2}.networkHistory.V_soma(:,tslice);

	% replayResults_clusters(sim{1});
	% replayResults_clusters(sim{1});
	disp('profiling...1')
	R{1} = profile_sim(Joinedsim{1},'tslice',tslice);
	disp('profiling...2')
	R{2} = profile_sim(Joinedsim{2},'tslice',tslice);

	set(0,'defaultaxescolororder', linspecer(10))
	set(0,'defaultfigurecolormap', linspecer(10))


end


if onesec_vs_30s
		Ronesec_withgap = profile_sim(Joinedsim{2},'tslice', [4000:5000]);
		stacked = vertcat(R{2}.allneurons,Ronesec_withgap.allneurons );
		G = [ones(200,1)*2 ;ones(200,1)*1];
		sel_fields = {'ampl', 'freq_each', 'spks', 'g_CaL', 'g_h'};
		sel_fields = {'ampl', 'freq_each', 'g_CaL', 'g_h'};
		NDscatter(stacked(:,sel_fields), G);
end

	

if plotcellscatters_gap_gapless 
	load(noiseless_200)
	% sel_fields = {'g_CaL', 'g_h',  'ampl', 'freq_each', 'meanVm' 'spks'};
	sel_fields = {'ampl', 'freq_each' , 'g_CaL'};
	% sel_fields = {'g_CaL', 'ampl', 'freq_each'}
	sel_table_1 = R{1}.allneurons(:,sel_fields);
	sel_table_2 = R{2}.allneurons(:,sel_fields);
	
	stacked = vertcat(sel_table_1, sel_table_2);
	G = [ones(200,1) ; ones(200,1)*2];

	NDscatter(stacked, G, 'colors', [0 163 218; 201 28 35]/255  )
	legend({'WT' ; 'Gdj2'})

end

sampletraces = 1;
if sampletraces
	replayResults_3(sims{1})
	replayResults_3(sims{2})
end



if sortedsampletraces

	load noiseless_200.mat

	colorder = flipud(cbrewer('seq', 'Greys', 30));
	colrmap = flipud(cbrewer('seq', 'Greys', 30));

	V1 = sims{1}.networkHistory.V_soma;
	V2 = sims{2}.networkHistory.V_soma;

	ampV1 = max(V1')-min(V1')
	ampV2 = max(V2')-min(V2')

	[V ordV1] = sort(ampV1)
	[V ordV2] = sort(ampV2)
	figure 
	imagesc(V1(ordV1,1000:1500))
	colormap(colrmap);
	caxis([-70 -40])
	axis off

	figure 
	imagesc(V2(ordV2,1000:1500))
	colormap(colrmap);
	caxis([-70 -40])
	axis off
	

	
	figure 
	set(0,'defaultaxescolororder', colorder)
	plot(V1(ordV1,1000:1500)')
	hold on , plot(mean(V1(ordV1,1000:1500)),'b','linewidth', 2)
	axis off, axis tight
	ylim([-70, -40])
	add_x_scalebar(gca, 100)


	figure 
	plot(V2(ordV2,1000:1500)')
	hold on , plot(mean(V2(ordV2,1000:1500)),'r','linewidth', 2)
	axis tight, axis off
	ylim([-70, -40])
	
	add_x_scalebar(gca, 100)

	saveallfigs('prefix', 'nonoise_traces_w_wo', 'style', '4x4')


end


noisesnippets = 1;

if noisesnippets
	for i = 1:5
		ounoise(i,:) = OUnoise('seed', i, 'plotme', 0, 'simtime', 500, 'dt', 1 ,'thetas', 1/30);
		figure
		plot(subplus(ounoise(i,:)),'r'); hold on;plot(-subplus(-ounoise(i,:)),'b');
		axis tight
		ylim([-10 10])
	end
	saveallfigs('prefix',  'noisesnips', 'style', '4x4');
end


if triggeredphase


	STPD{1} = stim_trig_phase_dist(Joinedsim{1});
	STPD{2} = stim_trig_phase_dist(Joinedsim{2});

	figure
	 plot_mean_and_std(STPD{1}.R1.kp_mask,'color' ,[1 0 0]); hold on
	 plot_mean_and_std(STPD{2}.R1.kp_mask,'color' ,[0 0 1])
	 alpha(.5)
	 legend({'MT' 'MT' 'WT' 'WT'})
	 title('kuramoto (stimulated mask)')
	 xlim([800 1500])

	 figure
	 plot_mean_and_std(STPD{1}.R1.Vm_othr,'color' ,[.5 .5 .5]); hold on
	 plot_mean_and_std(STPD{1}.R1.Vm_mask,'color' ,[1 0 0])
	 plot_mean_and_std(STPD{1}.R1.Vm_neig,'color' ,[0 1 0])
	 title('Stim Trig Average Vm (MT)')
	 alpha(.5)
	 xlim([800 1500])
	 
	 legend({'Other' 'Other' 'Mask' 'Mask' 'Neighbors' 'Neighbors'  })

	figure
	 plot_mean_and_std(STPD{2}.R1.Vm_othr,'color' ,[.5 .5 .5]); hold on
	 plot_mean_and_std(STPD{2}.R1.Vm_mask,'color' ,[1 0 0])
	 plot_mean_and_std(STPD{2}.R1.Vm_neig,'color' ,[0 1 0])
	 title('Stim Trig Average Vm (WT)')
	 alpha(.5)
	 xlim([800 1500])
	 
	legend({'Other' 'Other' 'Neighbors' 'Neighbors' 'Mask' 'Mask' })


end




if hist_sto_freq_amp_w_wo_gaps

	freqbins = [0:1:10];
	freqhistwithout  = hist(table2array(R{1}.allneurons(:,'freq_each')), freqbins)
	freqhistwith  = hist(table2array(R{2}.allneurons(:,'freq_each')), freqbins)

	ampbins = [0:2:25];
	amplitudehistwithout = hist(table2array(R{1}.allneurons(:,'ampl')),ampbins)
	amplitudehistwith 	 = hist(table2array(R{2}.allneurons(:,'ampl'))  ,ampbins)


	figure
		subplot(1,2,1)
		bar(freqbins, [ freqhistwithout ; freqhistwith]',1)
		xlabel('Freq (Hz) ')
		ylabel('Cells')
		title('STO freq')
		
		subplot(1,2,2)
		bar(ampbins, [amplitudehistwithout ; amplitudehistwith]',1)
		xlabel('Amplitude (mV) ')
		ylabel('Cells')
		title('STO amplitude')
end



if calculatesynchrony
	
	
		sim1_sync = measureGlobalSync(Joinedsim{1}, 'duration', 2000:25000, 'plotme',1, 'group', Joinedsim{1}.perturbation.mask{1});
		sim2_sync = measureGlobalSync(Joinedsim{2}, 'duration', 2000:25000, 'plotme',1, 'group', Joinedsim{2}.perturbation.mask{1});

		xlim([8600 9800])
		ninetyfive_1 = quantile(sim1_sync.stats.order_parameter_all(1:4900),.95)
		ninetyfive_2 = quantile(sim2_sync.stats.order_parameter_all(1:4900),.95)

		fig3 = figure;


		[HOPA_1 XOPA_1] = hist([sim1_sync.stats.order_parameter_all' ;  sim2_sync.stats.order_parameter_all']')
		[HOPA_2 XOPA_2] = hist([sim1_sync.stats.order_parameter_group' ;  sim2_sync.stats.order_parameter_group']')

		subplot(2,1,1)
		bar(XOPA_1,HOPA_1/(23000))
		title('order parameter (all)')
		legend({'MT'  'WT'})

		subplot(2,1,2)
		bar(XOPA_2,HOPA_2/(23000))
		title('order parameter (group)')
		legend({'MT'  'WT'})
		

end



if sto_and_propfiring_histograms


	F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
	% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0_1Hz_50000_4_25-Sep-2016.mat';


	addpath('/Users/M/Synced/Titan/Bench2/periodic_ampa/')
	addpath('/Users/M/Synced/Titan/Bench2/')
	addpath('/Users/M/Synced/Titan/Bench/')

	load(F1)
	numruns = 4;
	results_part{1} = profile_sim(simresults{1});
	results_part{2} = profile_sim(simresults{2});
	results_part{3} = profile_sim(simresults{3});
	results_part{4} = profile_sim(simresults{4});

	% STO FREQUENCIES
	stofreq_bins = [0:15];
	[HFsto_1 bins] = hist(table2array(results_part{1}.allneurons(:,'freq_each')),stofreq_bins); 
	[HFsto_2 bins] = hist(table2array(results_part{2}.allneurons(:,'freq_each')),stofreq_bins); 
	[HFsto_3 bins] = hist(table2array(results_part{3}.allneurons(:,'freq_each')),stofreq_bins); 
	[HFsto_4 bins] = hist(table2array(results_part{4}.allneurons(:,'freq_each')),stofreq_bins); 

	subplot(1,2,1)
	bar(stofreq_bins, [HFsto_1 ; HFsto_2 ; HFsto_3 ; HFsto_4]'/4,'stacked')
	title('STO in population')

	% SPIKE PROBABILITIES
	spkfreq_bins = [0:20];
	[HFspks_1 bins] = hist(table2array(results_part{1}.allneurons(:,'spks'))/50,spkfreq_bins); % 50second each part
	[HFspks_2 bins] = hist(table2array(results_part{2}.allneurons(:,'spks'))/50,spkfreq_bins); % 50second each part
	[HFspks_3 bins] = hist(table2array(results_part{3}.allneurons(:,'spks'))/50,spkfreq_bins); % 50second each part
	[HFspks_4 bins] = hist(table2array(results_part{4}.allneurons(:,'spks'))/50,spkfreq_bins); % 50second each part
	subplot(1,2,2)
	bar(spkfreq_bins, [HFspks_1 ; HFspks_2 ; HFspks_3 ; HFspks_4]'/4,'stacked')
	title('spike frequency')

	% SPIKE PROBABILITIES
	ampl_bins = [0:50];
	[HFamp_1 bins] = hist(table2array(results_part{1}.allneurons(:,'ampl')),ampl_bins); % 50second each part
	[HFamp_2 bins] = hist(table2array(results_part{2}.allneurons(:,'ampl')),ampl_bins); % 50second each part
	[HFamp_3 bins] = hist(table2array(results_part{3}.allneurons(:,'ampl')),ampl_bins); % 50second each part
	[HFamp_4 bins] = hist(table2array(results_part{4}.allneurons(:,'ampl')),ampl_bins); % 50second each part
	figure
	plot(ampl_bins, [HFamp_1 ; HFamp_2 ; HFamp_3 ; HFamp_4]'/4,'stacked')
	bar(ampl_bins, [HFamp_1 ; HFamp_2 ; HFamp_3 ; HFamp_4]'/4,'stacked')

	title('amplitude')


end 


if interperiodintervals
	netsize = [1 1 1];
	neurons = createDefaultNeurons(1);
	
	spont = 0; gap = eps;  noisesig = 0; noiseamp  = 0 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
	singlesim
	R{1} = simresults;
	
	spont = 1; gap = eps;  noisesig = .1; noiseamp = -.1 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
	singlesim
	R{2} = simresults;
	
	spont = 1; gap = eps;  noisesig = .2; noiseamp = -.2 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
	singlesim
	R{3} = simresults;
	
	spont = 1; gap = eps;  noisesig = .3; noiseamp = -.3 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
	singlesim
	R{4} = simresults;

	
	neurs = [1:(prod(netsize))];
	tslice = [1000:5000];
	H{1} = hilbert_of_membranepotential(R{1}.networkHistory.V_soma(neurs,tslice));
	H{2} = hilbert_of_membranepotential(R{2}.networkHistory.V_soma(neurs,tslice));
	H{3} = hilbert_of_membranepotential(R{3}.networkHistory.V_soma(neurs,tslice));
	H{4} = hilbert_of_membranepotential(R{4}.networkHistory.V_soma(neurs,tslice));

	N(1,:) = R{1}.networkHistory.backgroundnoise(1,:);
	N(2,:) = R{2}.networkHistory.backgroundnoise(1,:);
	N(3,:) = R{3}.networkHistory.backgroundnoise(1,:);
	N(4,:) = R{4}.networkHistory.backgroundnoise(1,:);

	for s = [1:4]
		ISIs = [];
		for n = [1]
			[pks tpks] = findpeaks(H{s}.hilbert(n,:),'minpeakdistance', 40, 'minpeakheight', 6);
			ISI = diff(tpks);
			isi_hist{s}(n,:) = hist(ISI,[0:10:200]);
		end
	end

	figure
	ax(1) = subplot(3,2,[1:2])
	ax(1).ColorOrder =  jet(3);
	plot([H{1}.hilbert(1,:); H{2}.hilbert(1,:); H{3}.hilbert(1,:)]')
	axis tight

	subplot(3,2,[3:4])
	VS = [R{1}.networkHistory.V_soma(1,:);R{2}.networkHistory.V_soma(1,:);R{3}.networkHistory.V_soma(1,:)];
	plot(VS')

	subplot(3,2,5)
	plot([0:10:200], isi_hist{1}','b')
	hold on
	plot([0:10:200], isi_hist{2}','g')
	plot([0:10:200], isi_hist{3}','r')
	
	figure
	X(1,:) = xcorr(VS(1,:), 'coeff');
	X(2,:) = xcorr(VS(2,:), 'coeff');
	X(3,:) = xcorr(VS(3,:), 'coeff');
	XN(1,:) = xcorr(N(1,:), 'coeff');
	XN(2,:) = xcorr(N(2,:), 'coeff');
	XN(3,:) = xcorr(N(3,:), 'coeff');

	plot(X'); hold on; 	
	plot(X'./max(X)'); hold on; 	


end
