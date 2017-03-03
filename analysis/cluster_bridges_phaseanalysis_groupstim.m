
  % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
  % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
 % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');
 % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_10-Feb-2017.mat');
% sim{1}.networkHistory.V_soma = gather(sim{1}.networkHistory.V_soma);
% sim{2}.networkHistory.V_soma = gather(sim{2}.networkHistory.V_soma);
% sim{3}.networkHistory.V_soma = gather(sim{3}.networkHistory.V_soma);



animate = 1;
duration = [1500:1900];
fname = 'clusters_curlies_bridges_22-Jan-2017';
fname = 'clusters_curlies_bridges_11-Feb-2017.mat';
fname = 'clusters_curlies_bridges_24-Feb-2017';


load([fname '.mat']);

for s = [1:4];

	phase_distribution_over_time(sim{s}, duration,1,  'group', find(sim{1}.W.stats.clusters==20),'animate',animate);
	if animate; eval(['!mv volume.mp4 ' fname '.mp4']);end
	if animate; eval(['!mv volume.avi ' fname '.avi']);end


	saveallfigs('prefix', [fname '_phase_group_sim1g20'])
	close all


end
