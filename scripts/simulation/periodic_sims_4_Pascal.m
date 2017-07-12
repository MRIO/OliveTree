% periodic_sims_4_Pascal.m


X_README = '1Hz_Periodic_4_Pascal_without_noise';
nameprefix = '1Hz_10s_no_noise';
seed = 0; tau = 30; noisesig =  0; noisemu = 0     ; sametoall = 0; simtype = '1Hz' ; gaps = [0.04] ; simtime = 20000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;
clear('all');

X_README = '1Hz_Periodic_4_Pascal_with_noise';
nameprefix = '1Hz_10s_with_noise';
seed = 0; tau = 30; noisesig =  -.4; noisemu = -.4 ; sametoall = .15; simtype = '1Hz' ; gaps = [0.04] ; simtime = 20000; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa;

