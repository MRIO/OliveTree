
% load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
% load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
% load('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_7-Feb-2017.mat');
% load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_10-Feb-2017.mat');
% sim{1}.networkHistory.V_soma = gather(sim{1}.networkHistory.V_soma);
% sim{2}.networkHistory.V_soma = gather(sim{2}.networkHistory.V_soma);
% sim{3}.networkHistory.V_soma = gather(sim{3}.networkHistory.V_soma);



animate = 0;
duration = [500:8000];
% fname = 'clusters_curlies_bridges_22-Jan-2017';
% fname = 'clusters_curlies_bridges_11-Feb-2017.mat';
% fname = 'clusters_curlies_bridges_24-Feb-2017';


fname = 'bridge_conductance_pspace01-Mar-2017';

load([fname '.mat']);

if 0
	for s = [1:11];
		sims{s}.W.coords = JM394_horizontal_coordinates;
		P{1} = phase_distribution_over_time(sims{s}, duration,1500,  'group', find(sims{s}.W.stats.clusters==20),'animate',animate);
		if animate; eval(['!mv volume.mp4 ' fname '.mp4']); end
		if animate; eval(['!mv volume.avi ' fname '.avi']); end


		% saveallfigs('prefix', [fname '_phase_group_sim1g20'])
		% close all
	end
end

if 1
	for s = 1:11
		for g = 1:50;
			groupsync{s,g} = measureGroupSync(sims{s}, 'duration',duration, 'plotme', 0, 'group', find(sims{s}.W.stats.clusters==g));

			% find(sims{s}.W.stats.clusters==g)

			G(s,g) = groupsync{s,g}.sync(1);
			V(s,g) = groupsync{s,g}.sync(2);
			MF(s,g) = mean(groupsync{s,g}.frequency);
			VF(s,g) = var(groupsync{s,g}.frequency);
			MK(s,g) = mean(mean(groupsync{s,g}.order_parameter));
			VK(s,g) = var(groupsync{s,g}.frequency);

		end
	end
end


% -- group 33, 35 desync



% only use those from groups with oscillating cells

% for nsims = 1:5

%     z = sims{nsims}.networkHistory.V_soma;

%     tt = 2001:10000;

%     % Remove the beginning for better processing
%     z1 = double(z(:, tt)');
%     % Data without a dc component
%     z2 = bsxfun(@minus, z1, mean(z1,1));

     
%     % Here we first remove the non-oscllating cells.
%     z3 = hilbert(z2);
%     logamp = log(mean(abs(z3)));
%     if plot_all
%         hist(logamp, 500)
%         box off
%     end

%     % Two peaks at the end are from non-oscillating cells. We remove them from our data...

     
%     igood = find(logamp>=-5.4);
%     ibad = find(logamp<-5.4);

