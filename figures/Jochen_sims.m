% Jochen_sims.m


X_README = 'XCorr_Jochen_ncorr_pspace'
thisnameprefix = 'xcorr_stim_ncorr_pspace'


	
	for val = [0 .15 .3 .45]
		nameprefix = [thisnameprefix '_ncorr_' num2str(val)  '_MT'];
		seed = 0; gaps = [eps]; simtime = 50000; tau = 30; noisesig = .4; noisemu = -.4; sametoall = val; simtype = '1Hz'; conntype = 'iso' ; numruns = 4;  nogapcomp = 15; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
		nameprefix = [thisnameprefix '_ncorr_' num2str(val)  '_WT'];
		seed = 0; gaps = [0.04]; simtime = 50000; tau = 30; noisesig = .4; noisemu = -.4; sametoall = val; simtype = '1Hz'; conntype = 'iso' ; numruns = 4; nogapcomp = 0; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
	end


X_README = 'FR variability tests'
thisnameprefix = 'firing rate variability'
	
	for val = [1 2 3 4]
		nameprefix = [thisnameprefix '_ncorr_' num2str(val)  '_MT'];
		seed = val; gaps = [eps]; simtime = 50000; tau = 30; noisesig = .4; noisemu = -.4; sametoall = .15; simtype = '1Hz'; conntype = 'iso' ; numruns = 4;  nogapcomp = 15; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
		nameprefix = [thisnameprefix '_ncorr_' num2str(val)  '_WT'];
		seed = 0; gaps = [0.04]; simtime = 50000; tau = 30; noisesig = .4; noisemu = -.4; sametoall = .15; simtype = '1Hz'; conntype = 'iso' ; numruns = 4; nogapcomp = 0; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
	end

