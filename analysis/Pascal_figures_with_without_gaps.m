% Pascal compare oscillation properties between gap and no gaps

clear;

oscillatingcells_comp = 1;
calculatesynchrony = 1;
sto_and_propfiring_histograms = 1;

addpath('/Users/M/Projects/Experiments/Olive/model/simresults/periodic_ampa')

if oscillatingcells_comp
	tslice = 1001:4000;

	if not(exist('sims'))
		% addpath('/Users/M/Projects/Experiements/Olive/model/simresults')
		% load('periodic_ampa_2_iso_0.04_spont_50000_2_12-Jun-2016.mat')
		load('periodic_ampa_replay_06_12_16_with_spont_gaptest2_iso_spont_5000_1_')
	end



	statevar{1} = simresults{1}.networkHistory.V_soma(:,tslice);
	statevar{2} = simresults{2}.networkHistory.V_soma(:,tslice);

	% replayResults_clusters(sim{1});
	% replayResults_clusters(sim{1});

	R{1} = profile_sim(simresults{1},'tslice',tslice);
	R{2} = profile_sim(simresults{2},'tslice',tslice);

	set(0,'defaultaxescolororder', linspecer(10))
	set(0,'defaultfigurecolormap', linspecer(10))



	freqbins = [0:1:30];
	freqhistwithout  = hist(table2array(R{1}.allneurons(:,'freq_each')), freqbins)
	freqhistwith  = hist(table2array(R{2}.allneurons(:,'freq_each')), freqbins)

	ampbins = [0:2:25];
	amplitudehistwithout = hist(table2array(R{1}.allneurons(:,'ampl')),ampbins)
	amplitudehistwith 	 = hist(table2array(R{2}.allneurons(:,'ampl'))  ,ampbins)



figure
	subplot(1,2,1)
	bar(freqbins, [ freqhistwithout ; freqhistwith]',1)
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO freq')
	
	subplot(1,2,2)
	bar(ampbins, [amplitudehistwithout ; amplitudehistwith]',1)
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO amplitude')


	if calculatesynchrony
		fig3 = figure;;
		
			sim1_sync = measureGlobalSync(simresults{1},'plotme',1);
			sim2_sync = measureGlobalSync(simresults{2},'plotme',1);
	end

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


