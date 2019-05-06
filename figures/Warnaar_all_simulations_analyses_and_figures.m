% this script reproduces multiple simulations / figures and analysis for Warnaar et al.
% 

% SET THE PATH TO FILES HERE:
pathtodata = '/Users/M/Synced/Projects/Experiments/Olive/model/simresults/periodic_ampa/'
% pathtodata = '.';
addpath(pathtodata)

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSE BELOW what analysis to REPRODUCE
%
% --- RECOMMENDATION to AVOID WORKSPACE MESS:
% on a clear workspace set onlye 
% one of the following variables to ONE AT A TIME!
% 
% also: uncomment lines in 'load_data' code block 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

produce_data = 0; % one of these two lines
load_data  = 1; % must be set to 1

profile_simulations = 0;
onesec_vs_30s = 0; % compare network behavior for different simulation durations
plotcellscatters_gap_gapless  = 0; % compare cell scatters in nets
joinedcellscatter = 0; % scatter of cell properties for joined simulations (200s)
triggeredphase = 0; % calculate triggered phase for 1Hz stimulation with and without gaps
hist_sto_freq_amp_w_wo_gaps = 0;
frequency_drift_singlesim_prc = 0; % demonstration of frequency drifts 
calculatesynchrony = 0;
sto_and_propfiring_histograms = 0;
interperiodintervals = 0;
sortedsampletraces = 0;
sampletraces = 0;
noisesnippets = 0; 
comparison_boundarycells = 0;
comparison_resonance_cav31 = 0;



if produce_data
	disp('To produce original data, uncomment the lines below (this runs HPCGPU_periodic_ampa with appropriate parameters);')

% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa	 	
% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		 
% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% load data and join  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if load_data


	% to produce figures with stimulation, uncomment:
		% F1 = 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017.mat';
	% to produce figures for spontaneous, uncomment:
		F1 = 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_spont_50000_4_17-Jan-2017.mat';


	load(F1)
	disp('loaded.')

	tslice = 1:50000;

	runs_with   = [5:8];
	runs_without = [1:4];

	Joinedsim{1}  = joinsim(simresults,runs_without); 
	Joinedsim{2}  = joinsim(simresults,runs_with); 

	statevar{1} = Joinedsim{1}.networkHistory.V_soma(:,tslice);
	statevar{2} = Joinedsim{2}.networkHistory.V_soma(:,tslice);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	%% profile simulations 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% profile simulations for this interval


	% replayResults_clusters(sim{1});
	% replayResults_clusters(sim{1});
	disp('profiling...without gaps')
	R{1} = profile_sim(Joinedsim{1},'tslice',tslice);
	disp('profiling...with gaps')
	R{2} = profile_sim(Joinedsim{2},'tslice',tslice);

	set(0,'defaultaxescolororder', linspecer(10))
	set(0,'defaultfigurecolormap', linspecer(10))


end



if onesec_vs_30s
		tslice = 1:50000;

		Ronesec_withgap = profile_sim(Joinedsim{2},'tslice', [4000:5000]);
		stacked = vertcat(R{2}.allneurons,Ronesec_withgap.allneurons );
		G = [ones(200,1)*2 ;ones(200,1)*1];
		sel_fields = {'ampl', 'freq_each', 'spks', 'g_CaL', 'g_h'};
		sel_fields = {'ampl', 'freq_each', 'g_CaL', 'g_h'};
		NDscatter(stacked(:,sel_fields), G);
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


	% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
	% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0_1Hz_50000_4_25-Sep-2016.mat';


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
	% plot(ampl_bins, [HFamp_1 ; HFamp_2 ; HFamp_3 ; HFamp_4]'/4,'stacked')
	bar(ampl_bins, [HFamp_1 ; HFamp_2 ; HFamp_3 ; HFamp_4]'/4,'stacked')

	title('amplitude')


end 


% if interperiodintervals
% 	netsize = [1 1 1];
% 	neurons = createDefaultNeurons(1);
	
% 	spont = 0; gap = eps;  noisesig = 0; noiseamp  = 0 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
% 	singlesim
% 	R{1} = simresults;
	
% 	spont = 1; gap = eps;  noisesig = .1; noiseamp = -.1 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
% 	singlesim
% 	R{2} = simresults;
	
% 	spont = 1; gap = eps;  noisesig = .2; noiseamp = -.2 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
% 	singlesim
% 	R{3} = simresults;
	
% 	spont = 1; gap = eps;  noisesig = .3; noiseamp = -.3 ; tau = 20; sametoall = 0.0; spont = 1; conntype = 'iso' ;  gapcomp = 0;
% 	singlesim
% 	R{4} = simresults;

	
% 	neurs = [1:(prod(netsize))];
% 	tslice = [1000:5000];
% 	H{1} = hilbert_of_membranepotential(R{1}.networkHistory.V_soma(neurs,tslice));
% 	H{2} = hilbert_of_membranepotential(R{2}.networkHistory.V_soma(neurs,tslice));
% 	H{3} = hilbert_of_membranepotential(R{3}.networkHistory.V_soma(neurs,tslice));
% 	H{4} = hilbert_of_membranepotential(R{4}.networkHistory.V_soma(neurs,tslice));

% 	N(1,:) = R{1}.networkHistory.backgroundnoise(1,:);
% 	N(2,:) = R{2}.networkHistory.backgroundnoise(1,:);
% 	N(3,:) = R{3}.networkHistory.backgroundnoise(1,:);
% 	N(4,:) = R{4}.networkHistory.backgroundnoise(1,:);

% 	for s = [1:4]
% 		ISIs = [];
% 		for n = [1]
% 			[pks tpks] = findpeaks(H{s}.hilbert(n,:),'minpeakdistance', 40, 'minpeakheight', 6);
% 			ISI = diff(tpks);
% 			isi_hist{s}(n,:) = hist(ISI,[0:10:200]);
% 		end
% 	end

% 	figure
% 	ax(1) = subplot(3,2,[1:2])
% 	ax(1).ColorOrder =  jet(3);
% 	plot([H{1}.hilbert(1,:); H{2}.hilbert(1,:); H{3}.hilbert(1,:)]')
% 	axis tight

% 	subplot(3,2,[3:4])
% 	VS = [R{1}.networkHistory.V_soma(1,:);R{2}.networkHistory.V_soma(1,:);R{3}.networkHistory.V_soma(1,:)];
% 	plot(VS')

% 	subplot(3,2,5)
% 	plot([0:10:200], isi_hist{1}','b')
% 	hold on
% 	plot([0:10:200], isi_hist{2}','g')
% 	plot([0:10:200], isi_hist{3}','r')
	
% 	figure
% 	X(1,:) = xcorr(VS(1,:), 'coeff');
% 	X(2,:) = xcorr(VS(2,:), 'coeff');
% 	X(3,:) = xcorr(VS(3,:), 'coeff');
% 	XN(1,:) = xcorr(N(1,:), 'coeff');
% 	XN(2,:) = xcorr(N(2,:), 'coeff');
% 	XN(3,:) = xcorr(N(3,:), 'coeff');

% 	plot(X'); hold on; 	
% 	plot(X'./max(X)'); hold on; 	


% end





if plotcellscatters_gap_gapless 
	load noiseless_200
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

if sampletraces
	t_slice = [9750:11350];
	t_slice = [11750:13350];
	t_slice = [12750:14350];
	close all
	load periodic_ampa_1Hz_20s_ou_input_24Jun-2017.mat
	figure
	replayResults_3(simresults{1}, t_slice)
	load periodic_ampa_1Hz_20s_no_noise_24Jun-2017.mat
	figure
	replayResults_3(simresults{1}, t_slice)
	saveallfigs('prefix', 'Pascal_traces', 'style', '1col')
	close all
end



if sortedsampletraces

	load noiseless_200.mat

	colorder = flipud(cbrewer('seq', 'Greys', 30));
	colrmap = flipud(cbrewer('seq', 'Greys', 60));

	V1 = sims{1}.networkHistory.V_soma;
	V2 = sims{2}.networkHistory.V_soma;

	ampV1 = max(V1')-min(V1')
	ampV2 = max(V2')-min(V2')

	[V ordV1] = sort(ampV1)
	[V ordV2] = sort(ampV2)
	
	interv = [2000 2500];
	figure 
	imagesc(V1(ordV1,:))
	colormap(colrmap);
	caxis([-70 -40])
	axis off
	xlim(interv)

	figure 
	imagesc(V2(ordV2,:))
	colormap(colrmap);
	caxis([-70 -40])
	axis off
	xlim(interv)
	

	
	figure 
	set(0,'defaultaxescolororder', colorder)
	plot(V1(ordV1,:)')
	hold on , plot(mean(V1(ordV1,:)),'b','linewidth', 2)
	axis off, axis tight
	ylim([-70, -40])
	addScalebar(gca, [100 10])
	xlim(interv)


	figure 
	plot(V2(ordV2,:)')
	hold on , plot(mean(V2(ordV2,:)),'r','linewidth', 2)
	axis tight, axis off
	ylim([-70, -40])
	xlim(interv)
	
	addScalebar(gca, [100 10])

	saveallfigs('prefix', 'nonoise_traces_w_wo', 'style', '4x4')


end



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



if comparison_boundarycells



	netsize = [2 10 10];
		noneurons = prod(netsize);

	plotthis  = 0;
	rd = 2;
	meannoconn = 8;
	normleak  = 1;
	randomize = 1;
	scaling   = 0.04;
	maxiter	  = 1;
	somatapositions = [];
	randomize = 1;
	symmetrize = 1;

	W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);


	plotnetstruct(W.W,W.coords(:,1),W.coords(:,2),W.coords(:,3),W.stats.clustercoeff.bu)

	borders = W.coords(:,2)==1 | W.coords(:,2)==10 | W.coords(:,3)==1 | W.coords(:,3)==10 

	figure
	BU = W.stats.clustercoeff.bu;
	[h1 x] = hist(BU(borders));
	h2 = hist(BU(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('cluster coefficient')
	legend({'border' 'center'})
	[sig p] = kstest2(BU(borders), BU(~borders))
	title({'cluster coefficient (binary undirected)' ;  ['sig: ' num2str(sig) ' ; p-val: ' num2str(p)]})

	figure
	gapleaks = sum(W.W);
	[h1 x] = hist(gapleaks(borders));
	h2 = hist(gapleaks(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('gap leaks')
	legend({'border' 'center'})
	[sig p] = kstest2(gapleaks(borders), gapleaks(~borders))
	title({'gap leaks'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})
	
	figure
	conns = sum(W.W>0);
	[h1 x] = hist(conns(borders));
	h2 = hist(conns(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('number of connections')
	[sig p] = kstest2(conns(borders), conns(~borders))
	legend({'border' 'center' })
	title({'connections'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})
	
	figure
	medfr = simresults{1}.spikes.medfreq;
	[h1 x] = hist(medfr(borders));
	h2 = hist(medfr(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('median frequency')
	[sig p] = kstest2(medfr(borders), medfr(~borders))
	legend({'border' 'center' })
	title({'median frequency (Hz)'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})

	figure
	firingrate = simresults{1}.spikes.spikespercell/(50*4); % (4 runs of 50s each )
	[h1 x] = hist(firingrate(borders));
	h2 = hist(firingrate(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('Firing rate (Hz)')
	[sig p] = kstest2(firingrate(borders), firingrate(~borders))
	legend({'border' 'center' })
	title({'firing rate (Hz)'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})


end


if comparison_resonance_cav31

	clear
		load periodic_ampa_CaH_bump_1_iso_1Hz_50000_1_1_14-Dec-2018.mat
		X_1 = xcorr_summa(simresults{1});
		save tmp X_1

	clear
		load tmp X_1
 		load periodic_ampa_Control_1_iso_1Hz_50000_1_1_13-Dec-2018.mat
		X_2 = xcorr_summa(simresults{1}, 'selectedneurons', X_1.selectedneurons)
		load tmp
		save tmp X_1 X_2

	close all

	plot(X_1.XcorrNoAc,'b')
 	hold on
	plot(X_2.XcorrNoAc,'r')
	legend({'CaH bump' 'control'})

	figure
	subplot(1,2,1)
	imagesc([-400,400], [1:190, X_1.XC{1})
	title('CaH bump')
	colorbar

	subplot(1,2,2)
	imagesc(X_2.XC{1})
	title('control')
	colorbar


end


