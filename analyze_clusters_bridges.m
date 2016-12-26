% analyze_clusters_bridges.m

tslice = 1001:3000;

load('/Users/M/Synced/Titan/clusters_curlies_bridges_20-Dec-2016.mat')

sim{1}.W = bridg_curlies;
sim{2}.W = curlies;

statevar{1} = sim{1}.networkHistory.V_soma(:,tslice);

% replayResults_clusters(sim{1});
% replayResults_clusters(sim{1});

R{1} = profile_sim(sim{1},'tslice',tslice);
% R{2} = profile_sim(sim{2});


% clusters
clusters = sim{1}.W.stats.clusters;
if ~exist('cbrewer')
		lc = jet(no_clusters);
	else
		lc = cbrewer('qual', 'Set1', no_clusters);
end


% bridge cells and neighbors:
bridgecells = bc;

	for c = find(bc)'
		figure
		thisbridge = zeros(size(bc));
		thisbridge(c) = 1;
		neighbors{c} = find(thisbridge'*(sim{1}.W.W>0));

		plot(tslice, statevar{1}(neighbors{c},:),'color', [1 1 1]*.9)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(tslice, statevar{1}(c,:),'linewidth', 2)
		pause

		quantile(table2array(R{1}.allneurons(neighbors{c},'freq_each')), [.25, .5, .75])

		hist(table2array(R{1}.allneurons(:,'freq_each')))
		hist(table2array(R{1}.allneurons(:,'ampl')))

	end


hist(table2array(R{1}.allneurons(logical(~bridgecells),'ampl')))
hist(table2array(R{1}.allneurons(logical(bridgecells),'ampl')))




% relationship between cluster amplitude and oscillator amplitude

% relationship between single cell frequency and oscilltor frequency

% probability of non oscillating cluster becoming oscillating cluster

% cluster synchrony distribution with and without bridges

	% mean and std of 2000s of order parameter


% activity of bridge cells and neighbors


	for c = 1:no_clusters
		c
		% clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
		% clustered{c}.no_neurons = length(find(clusters==c));



		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot([1:simtime], mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		pause
	end



% show time for group and global to align



% [=================================================================]
%  Bridges and Curlies
% [=================================================================]


% [=================================================================]
%  Cluster Activity Stats
% [=================================================================]



if calculatesynchrony
	fig3 = figure;;
	for c = 1:no_clusters
		c
		clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
		clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot([1:simtime], mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		pause
	end
		% xlabel('time (ms)')
		% ylabel('mV')
end


% sync
fig3 = figure;;
	for c = 1:no_clusters
		c
		clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
		clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot([1:simtime], mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		pause
	end



