
  % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
  % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
 % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');
 % load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_10-Feb-2017.mat');
% sim{1}.networkHistory.V_soma = gather(sim{1}.networkHistory.V_soma);
% sim{2}.networkHistory.V_soma = gather(sim{2}.networkHistory.V_soma);
% sim{3}.networkHistory.V_soma = gather(sim{3}.networkHistory.V_soma);

load('clusters_curlies_bridges_10-Feb-2017.mat');


R{1} = profile_sim(sim{1});
R{2} = profile_sim(sim{2});
R{3} = profile_sim(sim{3});

sel_fields = {'g_CaL', 'g_int', 'g_h', 'ampl', 'freq_each', 'meanVm','minV'};
sel_fields = { 'freq_each' 'g_CaL', 'g_h', 'ampl' }


NDscatter(R{1}.allneurons(:,sel_fields),1)
NDscatter(R{2}.allneurons(:,sel_fields),1)
NDscatter(R{3}.allneurons(:,sel_fields),1)

duration = [500:1500];
phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==5));
!mv volume.mp4 sim1group5.mp4
phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==5));
!mv volume.mp4 sim2group5.mp4
% phase_distribution_over_time(sim{3}, duration,1,  'group', find(sim{1}.W.stats.clusters==5));
% !mv volume.mp4 sim3group5.mp4


phase_distribution_over_time(sim{1}, duration,1,  'group', find(sim{1}.W.stats.clusters==10));
!mv volume.mp4 sim1group10.mp4
phase_distribution_over_time(sim{2}, duration,1,  'group', find(sim{1}.W.stats.clusters==10));
!mv volume.mp4 sim2group10.mp4
% phase_distribution_over_time(sim{3}, duration,1,  'group', find(sim{1}.W.stats.clusters==10));
% !mv volume.mp4 sim3group10.mp4

 phase_distribution_over_time(sim{1}, duration, 1, 'group', find(sim{1}.W.stats.clusters==20));
 !mv volume.mp4 sim1group20.mp4
 phase_distribution_over_time(sim{2}, duration, 1, 'group', find(sim{2}.W.stats.clusters==20) );
!mv volume.mp4 sim2group20.mp4

 % phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==20), [480:800],1, 'fname', 'sim3group20');
