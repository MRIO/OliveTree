% periodic_sims_4_Pascal.m


X_README = '1Hz_Periodic_4_Pascal';
thisnameprefix = '1Hz_10s_no_noise';
seed = 0; tau = 30; noisesig =  0; noisemu = 0     ; sametoall = 0; simtype = '1Hz' ; gaps = [0.04] ; simtime = 10000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;

X_README = '1Hz_Periodic_4_Pascal';
thisnameprefix = '1Hz_10s_with_noise';
seed = 0; tau = 30; noisesig =  -.4; noisemu = -.4 ; sametoall = .15; simtype = 'spont' ; gaps = [0.04] ; simtime = 10000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;

