
 load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');

 load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');





 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==5),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==5),'duration', [480:800])
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==5),'duration', [480:800])

 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==10),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==10),'duration', [480:800])
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==10),'duration', [480:800])

 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==20),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==20),'duration', [480:800])
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==20),'duration', [480:800])