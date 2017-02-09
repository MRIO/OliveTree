
 % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
 % load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');
 load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_09-Feb-2017.mat');
sim{1}.networkHistory.V_soma = gather(sim{1}.networkHistory.V_soma);
% sim{2}.networkHistory.V_soma = gather(sim{2}.networkHistory.V_soma);
% sim{3}.networkHistory.V_soma = gather(sim{3}.networkHistory.V_soma);



R{1} = profile_sim(sim{1});
% R{2} = profile_sim(sim{2})
% R{3} = profile_sim(sim{3})

sel_fields = {'g_CaL', 'g_int', 'g_h', 'ampl', 'freq_each', 'meanVm','minV'}

NDscatter(R{1}.allneurons(:,sel_fields),1)
% NDscatter(R{2}.allneurons(:,sel_fields),1)
% NDscatter(R{3}.allneurons(:,sel_fields),1)


 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==5),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==5),'duration', [480:800]);
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==5),'duration', [480:800]);

 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==10),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==10),'duration', [480:800]);
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==10),'duration', [480:800]);

 phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==20),'duration', [480:800]);
 phase_distribution_over_time(sim{2}, 'group', find(sim{2}.W.stats.clusters==20),'duration', [480:800]);
 phase_distribution_over_time(sim{3}, 'group', find(sim{3}.W.stats.clusters==20),'duration', [480:800]);

