
  % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
  % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
 % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');
 % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_10-Feb-2017.mat');
% sim{1}.networkHistory.V_soma = gather(sim{1}.networkHistory.V_soma);
% sim{2}.networkHistory.V_soma = gather(sim{2}.networkHistory.V_soma);
% sim{3}.networkHistory.V_soma = gather(sim{3}.networkHistory.V_soma);




animate = 0;
duration = [500:1500];




load('clusters_curlies_bridges_22-Jan-2017.mat');

% phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==5),'animate', animate);
% % !mv volume.mp4 sim1group5.mp4
% saveallfigs('prefix', 'phase_group_sim1g5_22Jan17')
% close all

% phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==10),'animate',animate);
% % !mv volume.mp4 sim1group10.mp4
% saveallfigs('prefix', 'phase_group_sim1g10_22Jan17')
% close all

phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==20),'animate',animate);
% !mv volume.mp4 sim1group20.mp4
saveallfigs('prefix', 'phase_group_sim1g20_22Jan17')
close all

% phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==5),'animate', animate);
% % !mv volume.mp4 sim1group5.mp4
% saveallfigs('prefix', 'phase_group_sim2g5_22Jan17')
% close all

% phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==10),'animate',animate);
% % !mv volume.mp4 sim1group10.mp4
% saveallfigs('prefix', 'phase_group_sim2g10_22Jan17')
% close all

phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==20),'animate',animate);
% !mv volume.mp4 sim1group20.mp4
saveallfigs('prefix', 'phase_group_sim2g20_22Jan17')
close all


load('clusters_curlies_bridges_11-Feb-2017.mat');


% phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==5),'animate', animate);
% % !mv volume.mp4 sim1group5.mp4
% saveallfigs('prefix', 'phase_group_sim1g5_11Feb17')
% close all

% phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==10),'animate',animate);
% % !mv volume.mp4 sim1group10.mp4
% saveallfigs('prefix', 'phase_group_sim1g10_11Feb17')
% close all

phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==20),'animate',animate);
% !mv volume.mp4 sim1group20.mp4
saveallfigs('prefix', 'phase_group_sim1g20_11Feb17')
close all

% phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==5),'animate', animate);
% % !mv volume.mp4 sim1group5.mp4
% saveallfigs('prefix', 'phase_group_sim2g5_11Feb17')
% close all

% phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==10),'animate',animate);
% % !mv volume.mp4 sim1group10.mp4
% saveallfigs('prefix', 'phase_group_sim2g10_11Feb17')
% close all

phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==20),'animate',animate);
% !mv volume.mp4 sim1group20.mp4
saveallfigs('prefix', 'phase_group_sim2g20_11Feb17')
close all

