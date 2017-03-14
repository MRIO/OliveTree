% analyze_clusters_bridges.m


% [================================================]
%  what to do
% [================================================]

plotbridgeandneighbors = 1;
	plotcellscatters  = 0;
	joinedcellscatter  = 0;
plot_cluster_members = 1;
STOhistograms = 0;
calculatesynchrony = 0;


tslice = 1001:4000;

if not(exist('st_st'))
	% load('/Users/M/Public/Dropbox/simresults/clusters_curlies_bridges_26-Dec-2016.mat')
	load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
	load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
    sims = sim;
end

if not(exist('sims'))
	sims{1}.W = bridg_curlies;
	sims{2}.W = curlies;

	statevar{1} = double(sims{1}.networkHistory.V_soma(:,tslice)); % with gaps 
	statevar{2} = double(sims{2}.networkHistory.V_soma(:,tslice)); % without gaps
	statevar{3} = double(sims{3}.networkHistory.V_soma(:,tslice)); % disconnected

	% replayResults_clusters(sim{1});
	% replayResults_clusters(sim{1});

	R{1} = profile_sim(sims{1},'tslice',tslice); % bridges and curlies
	R{2} = profile_sim(sims{2},'tslice',tslice); % only curlies
	R{3} = profile_sim(sims{3},'tslice',tslice); % disconnected net

	% RR = vertcat(R{1}.allneurons, R{2}.allneurons);
	% RR(:,39) = table([zeros(1105,1) ; ones(1105,1)]);

	% sel_fields = {'g_CaL', 'g_int', 'p1',  'ampl', 'freq_each', 'meanVm', 'Var39'};
	% sel_table = RR(:,sel_fields);
	% NDscatter(sel_table, 7);

end






if plotcellscatters 
	sel_fields = {'g_CaL', 'g_int', 'p1', 'g_h',  'ampl', 'freq_each', 'meanVm'}
	sel_table = R{1}.allneurons(:,sel_fields);
	figure
	NDscatter(sel_table, 1)

	sel_table = R{3}.allneurons(:,sel_fields);
	figure	
	NDscatter(sel_table, 1)
end



if joinedcellscatter

	RR = vertcat(R{1}.allneurons, R{2}.allneurons);

	G = [ones(200,1) ;ones(200,1)*2];
	

	sel_fields = {'g_CaL', 'ampl', 'freq_each', 'g_int'};
	NDscatter(RR(:,sel_fields), G)


end


% [=================================================================]
%  set defaults
% [=================================================================]
% clusters
clusters = sims{1}.W.stats.clusters;
no_clusters = length(unique(clusters));
if ~exist('cbrewer')
		lc = jet(no_clusters);
	else
		lc = cbrewer('qual', 'Set1', no_clusters);
end
set(0,'defaultaxescolororder', linspecer(10))
set(0,'defaultfigurecolormap', linspecer(10))


% [=================================================================]
%  connectivity
% [=================================================================]

% bridge cells and neighbors:
bridgecells = bc;

% effective number of connections clusters x bridges
[connhistC bins] =  hist(curlies.stats.connections,[0:20]);
[connhistB bins] =  hist(bridges.stats.connections(find(bc)),[0:20]);
figure
bar(bins, [connhistB; connhistC]','stacked')
title('out degree')
legend({'bridges' 'curlies'})

bridges_from_cluster = single(bc .* clusters==41);
neighbors_to_bridge = find(bridges_from_cluster'*(sims{1}.W.W>0));
their_cluster = unique(clusters(neighbors_to_bridge));
targeted_cluster_cells = ismember(clusters, their_cluster).*clusters;
% plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), targeted_cluster_cells);

plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), clusters==41 | clusters==34);
view(104,56.4);




% [=================================================================]
%  bridge behavior
% [=================================================================]
bridgebehaviorplot  = 0;
if bridgebehaviorplot
	bridgeidx = find(bridgecells);
	bridge_Vm = statevar{1}(bridgeidx,:);
	[val ord] = sort(table2array(R{1}.allneurons(bridgeidx,'ampl')));

	figure
	waterfall(bridge_Vm(ord,:));
	title('bridge behavior when connected')

	figure
	bridge_Vm = statevar{3}(bridgeidx,:);
	waterfall(bridge_Vm(ord,:));
end

% [=================================================================]
%  groups
% [=================================================================]

if plotbridgeandneighbors
	bridgestoplot = find(bc.*(clusters==41))';
	for c = bridgestoplot
		figure
		thisbridge = zeros(size(bc));
		thisbridge(c) = 1;
		neighbors{c} = find(thisbridge'*(sims{1}.W.W>0));

		plot(tslice, statevar{1}(neighbors{c},:),'color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(tslice, statevar{1}(c,:),'linewidth', 2)
		

		% Q = quantile(table2array(R{1}.allneurons(neighbors{c},'freq_each')), [.25, .5, .75])
		% text(tslice(end)*.9 , max(max(statevar{1}))*.9, num2str(Q))

		% pause

		% close 

	end

end




if STOhistograms
	freqbins = [0:1:30];
	freqhistCurlies = hist(table2array(R{1}.allneurons(logical(~bridgecells),'freq_each')), freqbins )
	freqhistBridges =hist(table2array(R{1}.allneurons(logical(bridgecells),'freq_each')) ,  freqbins)
	freqhistCurlies0 = hist(table2array(R{2}.allneurons(logical(~bridgecells),'freq_each')),freqbins )
	freqhistBridges0 =hist(table2array(R{2}.allneurons(logical(bridgecells),'freq_each')) , freqbins)

	figure
	subplot(211)
	bar(freqbins, [ freqhistCurlies ; freqhistBridges]','stacked')
	xlabel('Frequency (Hz)')
	ylabel('Cells')
	title('STO freq')

	subplot(212)
	bar(freqbins, [ freqhistCurlies0 ; freqhistBridges0]','stacked')
	xlabel('Frequency (Hz)')
	ylabel('Cells')
	title('STO freq')


	ampbins = [0:1:30];
	amplitudehistCurlies  = hist(table2array(R{1}.allneurons(logical(~bridgecells),'ampl')),ampbins)
	amplitudehistBridges  = hist(table2array(R{1}.allneurons(logical(bridgecells),'ampl'))  ,ampbins)
	amplitudehistCurlies0 = hist(table2array(R{2}.allneurons(logical(~bridgecells),'ampl')),ampbins)
	amplitudehistBridges0 = hist(table2array(R{2}.allneurons(logical(bridgecells),'ampl')) ,ampbins)

	figure
	subplot(211)
	bar(ampbins, [ amplitudehistCurlies ; amplitudehistBridges]','stacked')
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO amplitude')

	subplot(212)
	bar(ampbins, [ amplitudehistCurlies0 ; amplitudehistBridges0]','stacked')
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO amplitude')
end

% relationship between cluster amplitude and oscillator amplitude

% relationship between single cell frequency and oscillator frequency

% probability of non oscillating cluster becoming oscillating cluster

% cluster synchrony distribution with and without bridges

	% mean and std of 2000s of order parameter


% activity of bridge cells and neighbors

if plot_cluster_members
	cc = 0;
	figure
	% their_cluster = 34;
	for c = their_cluster' %1:no_clusters
		cc = cc +1;
		figure
		cellsincluster = statevar{1}(find(clusters==c),:);
		cellsincluster = setdiff(cellsincluster, find(bc));
		% clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
		% clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot(tslice, mean(statevar{1}(find(clusters==c),:))+cc*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		% pause
	end
end



% show time for group and global to align



% [=================================================================]
%  Bridges and Curlies
% [=================================================================]


% [=================================================================]
%  Cluster Activity Stats
% [=================================================================]


tslice = 1:1000;

if calculatesynchrony
    
    %%
	fig3 = figure;;
	for c = 1:no_clusters
		c
		clustered{c}.sync = measureGroupSync(sims{1},'group', clusters==c,'plotme',1);
		clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot(tslice, mean(statevar{1}(find(clusters==c),tslice))+c*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		pause
	end
		% xlabel('time (ms)')
		% ylabel('mV')
        
        %%
end





