

% analyze_clusters_bridges.m

plotbridgeandneighbors = 1;
plotconnectivity  = 0;
plotscatters = 0;
plotcellscatters  = 0;
STOhistograms = 1;
calculatesynchrony = 0;
plot_cluster_members = 1;
plotreconstruction = 1;
	makevideo = 0;

tslice = 1500:5000;

if not(exist('st_st'))
	% load('/Users/M/Public/Dropbox/simresults/clusters_curlies_bridges_26-Dec-2016.mat')
	% load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_22-Jan-2017.mat');
	% load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_20-Jan-2017.mat');
	% load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_10-Feb-2017.mat')
	load('/Users/M/Projects/Experiments/Olive/model/simresults/clusters_curlies_bridges_11-Feb-2017.mat')
    sims = sim;
end

% [================================================]
%  profile simulations
% [================================================]


if not(exist('R'))
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



% [=================================================================]
%  connectivity
% [=================================================================]


clusters = sims{1}.W.stats.clusters;
no_clusters = length(unique(clusters));

if ~exist('cbrewer')
		lc = jet(no_clusters);
	else
		lc = cbrewer('qual', 'Paired', no_clusters);
end
set(0,'defaultaxescolororder', linspecer(10))
set(0,'defaultfigurecolormap', linspecer(10))

% bridge cells and neighbors:
bridgecells = bc;
if plotconnectivity
	% effective number of connections clusters x bridges
	[connhistC bins] =  hist(curlies.stats.connections,[0:20]);
	[connhistB bins] =  hist(bridges.stats.connections(find(bc)),[0:20]);
	figure


	B = bar(bins, [connhistB; connhistC]','stacked','edgecolor', 'none');

	B(1).FaceColor  = [.9 0 .2 ];
	B(2).FaceColor = [ .2 .7 .2];
	title('degree')
	xlabel('number of connections')
	legend({'bridges' 'curlies'})
end


if plotreconstruction

	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridgecells+1)
	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridg_curlies.stats.clusters)

	
	if makevideo
		clear vidObj
		vidObj = VideoWriter('this4','MPEG-4');
		vidObj.FrameRate = 24;
		vidObj.Quality = 90;
		title([])
		
		open(vidObj);
		f = 0; 
		for az = -180:180;
		 f=f+1; 
		 view(az, 45);
		 axis vis3d
		 drawnow;currframe(f) = getframe(gcf, [50 100 1000 700]);  
		 writeVideo(vidObj,currframe(f));
		end
		
		close(vidObj);
	end


end


% [=================================================================]
%  bridge behavior
% [=================================================================]
bridgeidx = find(bridgecells);
bridge_Vm = statevar{1}(bridgeidx,:);
[val ord] = sort(table2array(R{1}.allneurons(bridgeidx,'ampl')));

% [================================================]
%  cell scatters
% [================================================]


if plotcellscatters 
	% sel_fields = {'g_CaL', 'g_int', 'p1', 'g_h', 'g_ld',  'freq_each', 'meanVm','minV', 'supth'}
	% sel_fields = {'g_CaL', 'g_int', 'g_h', 'ampl', 'freq_each', 'meanVm','minV'}
	% sel_fields = {'g_CaL', 'g_ld',  'g_K_Ca', 'p1', 'g_h',  'ampl', 'freq_each', 'meanVm','minV','supth'}
	% sel_fields = {'g_CaL', 'g_ld', 'g_K_Ca',  'freq_each', 'meanVm','minV'}

	sel_fields = { 'freq_each' 'g_CaL', 'g_h', 'ampl' };

	sel_fields = { 'freq_each' 'g_CaL' 'ampl'};
	
	% sel_table = R{1}.allneurons(:,sel_fields);
	% NDscatter(sel_table, clusters)
	
	% sel_table = R{2}.allneurons(:,sel_fields);
	% NDscatter(sel_table, clusters)

	% sel_table = R{3}.allneurons(:,sel_fields);
	% NDscatter(sel_table, clusters)


	sel_table = R{1}.allneurons(:,sel_fields);
	NDscatter(sel_table, bridgecells+1)
	
	sel_table = R{2}.allneurons(:,sel_fields);
	NDscatter(sel_table, bridgecells+1)

	sel_table = R{3}.allneurons(:,sel_fields);
	NDscatter(sel_table, bridgecells+1)

end



figure
waterfall(bridge_Vm(ord,:));
title('bridge behavior when connected')

figure
bridge_Vm = statevar{3}(bridgeidx,:);
waterfall(bridge_Vm(ord,:));



% [=================================================================]
%  groups
% [=================================================================]

if plotbridgeandneighbors
	
	bc_index = find(bc)'; ind = 0;
	for c = bc_index(1:40);
		ind = ind+1;
		figure(c)
		thisbridge = zeros(size(bc));
		thisbridge(c) = 1;
		neighbors{c} = find(thisbridge'*(sims{1}.W.W>0));
		subplot(2,1,1)
		plot(statevar{1}(neighbors{c}',:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot( statevar{1}(c,:),'linewidth', 2)
		title({['bridge number: ' num2str(c) ' degree:' num2str(length(neighbors{c}))] ; ['neighbor freq (lower median upper): ' num2str(Q) 'Hz'] })

		subplot(2,1,2)
		plot(statevar{2}(neighbors{c},:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(statevar{2}(c,:),'linewidth', 2)
		title(['neighboring clusters:' num2str(clusters(neighbors{c})') ])

		Q = quantile(table2array(R{3}.allneurons(neighbors{c},'freq_each')), [.1, .5, .9])
		

		pause
		% close 

	end



end


if STOhistograms
	freqbins = [0:.5:15];
	freqhistCurlies = hist(table2array(R{1}.allneurons(logical(~bridgecells),'freq_each')), freqbins )
	freqhistBridges =hist(table2array(R{1}.allneurons(logical(bridgecells),'freq_each')) ,  freqbins)
	freqhistCurlies0 = hist(table2array(R{2}.allneurons(logical(~bridgecells),'freq_each')),freqbins )
	freqhistBridges0 =hist(table2array(R{2}.allneurons(logical(bridgecells),'freq_each')) , freqbins)

	figure
	subplot(211)
	set(gca,'colororder',[.9 0 .2 ; .2 .7 .2])
	bar(freqbins, [ freqhistCurlies ; freqhistBridges]','stacked')
	xlabel('Frequency (Hz)')
	ylabel('Cells')
	title('STO freq')

	subplot(212)
	set(gca,'colororder',[.9 0 .2 ; .2 .7 .2])
	bar(freqbins, [ freqhistCurlies0 ; freqhistBridges0]','stacked')
	xlabel('Frequency (Hz)')
	ylabel('Cells')
	title('STO freq')


	ampbins = [0:2:100];
	amplitudehistCurlies  = hist(table2array(R{1}.allneurons(logical(~bridgecells),'ampl')),ampbins)
	amplitudehistBridges  = hist(table2array(R{1}.allneurons(logical(bridgecells),'ampl'))  ,ampbins)
	amplitudehistCurlies0 = hist(table2array(R{2}.allneurons(logical(~bridgecells),'ampl')),ampbins)
	amplitudehistBridges0 = hist(table2array(R{2}.allneurons(logical(bridgecells),'ampl')) ,ampbins)

	figure
	subplot(211)
	set(gca,'colororder',[.9 0 .2 ; .2 .7 .2])
	bar(ampbins, [ amplitudehistCurlies ; amplitudehistBridges]','stacked')
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO amplitude')

	subplot(212)
	set(gca,'colororder',[.9 0 .2 ; .2 .7 .2])
	bar(ampbins, [ amplitudehistCurlies0 ; amplitudehistBridges0]','stacked')
	xlabel('Amplitude (mV) ')
	ylabel('Cells')
	title('STO amplitude')
end

boxplots = 1;
if boxplots
	BoxAmp = [  table2array(R{3}.allneurons(:,'ampl'))';
				table2array(R{2}.allneurons(:,'ampl'))';
				table2array(R{1}.allneurons(:,'ampl'))']';

	BoxFreq = [  table2array(R{3}.allneurons(:,'freq_each'))';
			 	 table2array(R{2}.allneurons(:,'freq_each'))';
			 	 table2array(R{1}.allneurons(:,'freq_each'))']';

			 % figure
			 % boxplot(BoxAmp(logical(~bridgecells),:), 'plotstyle', 'compact', 'medianstyle', 'line', 'notch', 'on',  'color', [.9 0 .2], 'boxstyle', 'filled'); hold on
			 % boxplot(BoxAmp(logical(bridgecells),:) , 'plotstyle', 'compact', 'medianstyle', 'line', 'notch', 'on', 'color', [.2 .7 .2], 'boxstyle', 'filled')
			 % set(gca,'xticklabel',{'both' 'curlies' 'disconnected'})
			 % ylabel('mV')
			 
			 figure
			 plot(BoxAmp(logical(~bridgecells),:)','color',[.9 0 .2 ],'marker', 'o');hold on
			 plot(BoxAmp(logical(bridgecells),:)' ,'color',[.2 .7 .2],'marker', 'x')
			 set(gca,'xticklabel',{'disconnected' 'curlies' 'curlies+bridges'},'xtick', [1 2 3])
			 ylabel('mV')

			 % figure
			 % boxplot(BoxFreq(logical(~bridgecells),:),'plotstyle', 'compact', 'medianstyle', 'line', 'notch', 'on','color', [.9 0 .2], 'boxstyle', 'filled' ,'plotstyle', 'compact'); hold on
			 % boxplot(BoxFreq(logical(bridgecells),:) ,'plotstyle', 'compact', 'medianstyle', 'line', 'notch', 'on','color', [.2 .7 .2], 'boxstyle', 'filled','plotstyle', 'compact')
			 % ylabel('Hz')

			 figure
			 plot(BoxFreq(logical(~bridgecells),:)','color',[.9 0 .2 ],'marker', 'o');hold on
			 plot(BoxFreq(logical(bridgecells),:)' ,'color',[.2 .7 .2],'marker', 'x')
 			 set(gca,'xticklabel',{'disconnected' 'curlies' 'curlies+bridges'},'xtick', [1 2 3])
			 ylabel('Hz')
end



% relationship between cluster amplitude and oscillator amplitude

% relationship between single cell frequency and oscilltor frequency

% probability of non oscillating cluster becoming oscillating cluster

% cluster synchrony distribution with and without bridges

	% mean and std of 2000s of order parameter


% activity of bridge cells and neighbors


% [=================================================================]
%  Cluster Activity Stats
% [=================================================================]

if plot_cluster_members
	for c = 1:max(clusters)
		c
		% clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
		% clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		plot(tslice, mean(statevar{1}(find(clusters==c),:))+c*5,'color', lc(c,:))
		hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		pause
	end
end



% [=================================================================]
%  Synchrony stats for clusters
% [=================================================================]

tslice = 1:3000;

% calculates phase coherence of clustered cells
if calculatesynchrony_clusters
    
    %%
	fig3 = figure;;
	for c = 17 %1:max(clusters)
		c = 17;
		clustered{c,1}.sync = measureGroupSync(sims{1},'group', clusters==c,'plotme',1);
		clustered{c,2}.sync = measureGroupSync(sims{2},'group', clusters==c,'plotme',1);

		clustered{c}.no_neurons = length(find(clusters==c));

		
		% plot(tslice, mean(statevar{1}(find(clusters==c),tslice))+c*5,'color', lc(c,:))
		% hold on
		

		pause
	end
		% calculate mean sync per cluster -> look at variability
end


% example cluster
if exampleclustersync
	figure
	 [hclu xclu] = hist([clustered{17,1}.sync.order_parameter, clustered{17,2}.sync.order_parameter],50);
	 BBC = bar(xclu, hclu/5000, 2, 'edgecolor', 'none');
	 BBC(1).FaceColor  = [.9 0 .2 ];
	 BBC(2).FaceColor = [ .2 .7 .2];
	title('Synchrony')
	ylabel('P(K)/ms')
	xlabel('Kuramoto Parameter')


	figure
	[hclu xclu] = hist([clustered{17,1}.sync.instantaneousFrequency(:), clustered{17,2}.sync.instantaneousFrequency(:)],5000);
	 BBC = bar(xclu, hclu/(5000*5000), 2, 'edgecolor', 'none');
	 BBC(1).FaceColor  = [.9 0 .2 ];
	 BBC(2).FaceColor = [ .2 .7 .2];
	 ylabel('P(Freq/ms)')
	 xlabel('Freq (Hz)')
	 title('Spectrum')
end




% % sync
% fig3 = figure;;
% 	for c = 1:no_clusters
% 		c
% 		clustered{c}.sync = measureGroupSync(sims{1},'group', clusters==c,'plotme',0);
% 		clustered{c}.no_neurons = length(find(clusters==c));

% 		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
% 		plot([1:simtime], mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
% 		hold on
% 		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
% 		pause
% 	end




% M{1} = measureGlobalSync(sims{1},'duration',tslice,'group',clusters==5);
% M{2} = measureGlobalSync(sims{2},'duration',tslice,'group',clusters==5);
% M{3} = measureGlobalSync(sims{3},'duration',tslice,'group',clusters==5);
