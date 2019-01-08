% Warnaar_all_prc_figures.m

whattodo = 'prc'

switch whattodo

	case 'gallop_with_without_noise'
		% [=================================================================]
		%  stimulation with and without contextual input
		% [=================================================================]

		X_README = '1Hz_Periodic_4_Pascal_without_noise';
		nameprefix = '1Hz_10s_no_noise';
		seed = 0; tau = 30; noisesig =  0; noisemu = 0     ; sametoall = 0; simtype = '1Hz' ; gaps = [0.04] ; simtime = 20000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;
		clear('all');

		X_README = '1Hz_Periodic_4_Pascal_with_noise';
		nameprefix = '1Hz_10s_with_noise';
		seed = 0; tau = 30; noisesig =  -.4; noisemu = -.4 ; sametoall = .15; simtype = '1Hz' ; gaps = [0.04] ; simtime = 20000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;

	case 'prc'

		% [=================================================================]
		%  Phase response analysis
		% [=================================================================]

		generate_PRC_Pascal

	case 'structure and example activity'

		figure_model_Pascal

	case 'stimulus triggered phases'

		% 1. load data

		STPD{1} = stim_triggered_phase_dist(sim{1});
		STPD{2} = stim_triggered_phase_dist(sim{2});


	case 'profile cells'
		profile_sim(sim{1})
		profile_sim(sim{1})

	
	case 'network pspace'

		NetPspace;
		analysis_NetPspace;


	case 'cell scatters'


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

end
