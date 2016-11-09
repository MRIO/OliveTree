
% comparison between phases in triggered with and without gaps

% load /Users/M/Cabinet/SyncBox/Bench2/sim3D_long_Xcorr_onlyexcitation_1.mat
% S1 = sim3D;
% load /Users/M/Cabinet/SyncBox/Bench2/sim3D_long_Xcorr_onlyexcitation_12.mat
% S2 = sim3D;

for n = 1:4:12
	n
	% eval(['load /Users/M/Cabinet/SyncBox/Bench2/sim3D_long_Xcorr_onlyexcitation_' num2str(n) '.mat'])
	eval(['load /Users/M/Cabinet/SyncBox/Bench/stim5pa_1hz/sim3D_long_Xcorr_onlyexcitation_' num2str(n) '.mat'])
	sum(sum(sim3D.networkParameters.connectivityMatrix))

	R{n} = stim_trig_phase_dist(sim3D);
	D{n} = conductance_contribution_to_firingrates(sim);
	
end



if 0
	muISI1 = [D{1}.meanISI D{2}.meanISI D{3}.meanISI D{4}.meanISI];
	muISI2 = [D{5}.meanISI D{6}.meanISI D{7}.meanISI D{8}.meanISI];

	figure(1), hist(muISI1, 20)
	figure(2), hist(muISI2, 20)
end



