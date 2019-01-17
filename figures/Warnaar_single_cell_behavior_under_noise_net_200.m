% Warnaar_single_cell_behavior_under_noise_net_200.m

	load 'periodic_ampa_MT_1_iso_1Hz_5000_1_1_16-Jan-2019.mat'
	S{1} = simresults{1};
	res{1} = profile_sim(S{1},'plotme', 1);

	load 'periodic_ampa_WT_1_iso_1Hz_5000_1_1_16-Jan-2019.mat'
	S{2} = simresults{1};
	res{2} = profile_sim(S{2},'plotme', 1);

	r = vertcat(res{1}.allneurons , res{2}.allneurons)

	res{1} = profile_sim(S{1},'plotme', 1,'tslice',[1000:5000]);
	res{2} = profile_sim(S{2},'plotme', 1,'tslice',[1000:5000]);
	r = vertcat(res{1}.allneurons , res{2}.allneurons)
	close all, NDscatter(r(:,{'ampl', 'freq_each', 'g_CaL'}),[ones(1,200) ones(1,200)+1]')
	legend({'MT', 'WT'})



