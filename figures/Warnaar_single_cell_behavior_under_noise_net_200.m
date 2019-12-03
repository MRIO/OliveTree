% Warnaar_single_cell_behavior_under_noise_net_200.m


% to produce sims:

	      % val = 0.6; nameprefix = 'WT' ; seed = 0; Cm_bump = 0; simtime = 5000; tau = 20; gapcomp = 0; noisesig = val; noisemu = -val; sametoall = 0.2; simtype = '1Hz'; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa	
% gaps = 0; val = 0.7; nameprefix = 'MT' ; seed = 0; Cm_bump = 0; simtime = 5000; tau = 20; gapcomp = 1; noisesig = val; noisemu = -val; sametoall = 0.2; simtype = '1Hz'; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa	


	load 'periodic_ampa_test_WT_4_iso_spont_50000_4_1_30-Apr-2018'
	S{1} = simresults{1};
	res{1} = profile_sim(S{1},'plotme', 1);

	load 'periodic_ampa_test_MT_4_iso_spont_50000_4_1_30-Apr-2018'
	S{2} = simresults{1};
	res{2} = profile_sim(S{2},'plotme', 1);

	r = vertcat(res{1}.allneurons , res{2}.allneurons)

	res{1} = profile_sim(S{1},'plotme', 1,'tslice',[1000:5000]);
	res{2} = profile_sim(S{2},'plotme', 1,'tslice',[1000:5000]);
	r = vertcat(res{1}.allneurons , res{2}.allneurons)
	close all, NDscatter(r(:,{'ampl', 'freq_each', 'spks',  'g_ls' , 'g_ld'}),[ones(1,200)+1 ones(1,200)]')
	legend({'WT', 'MT'})



