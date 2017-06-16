

% analyze_clusters_bridges.m

profile_brick_curlies = 1;
plotcellscatters  = 1;	sel_fields = { 'g_CaL', 'g_h' 'ampl' 'freq_each'};


plotreconstruction = 0; 
plotselectedclusters = 0;
	makevideo = 0
plotconnectivityhistogram  = 0; % comparison of degree between clusters and bridges
distance_histogram = 0;

plotclustermemberaverages = 0;

plotbridgeandneighbors_Vm = 0;
	plotbridgewaterfall = 0;

plotclustercellactivity = 0;

analyze_group_stim = 1;
activitydifference = 1;
spectral_clustering = 1;


STOhistograms = 0;
calculatesynchrony_clusters = 0;
exampleclustersync = 0;

bridgecond_pspace = 0;

boxplots = 0;



% [================================================]
%  default colors
% [================================================]

bridgecolor = [.2 .2 .9 ];
curlycolor = [ .2 .7 .2];


if not(exist('st_st'))
	
	load('clusters_curlies_bridges_01-Mar-2017.mat')
    sims = sim;
end
 
% [================================================]
%  profile simulations
% [================================================]

if profile_brick_curlies
	if not(exist('R'))
		sims{1}.W = bridg_curlies;
		sims{2}.W = curlies;
		sims{3}.W = curlies;
		sims{4}.W = brick;
		tslice = 500:8000;

		for nsim = 1:length(sims)
			statevar{nsim} = double(sims{nsim}.networkHistory.V_soma(:,tslice)); % with gaps 
			R{nsim} = profile_sim(sims{nsim},'tslice',tslice); % bridges and curlies
		end
		
	end
end

% [=================================================================]
%  connectivity
% [=================================================================]

bridgecells = bc;
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

if plotconnectivityhistogram
	% effective number of connections clusters x bridges
	bins = [0:20];
	
	[connhistC bins] =  hist(curlies.stats.connections,[0:20]);
	[connhistB bins] =  hist(bridges.stats.connections(find(bc)),[0:20]);
	
	figure
	B = bar(bins, [connhistB; connhistC]','stacked','edgecolor', 'none');



	B(1).FaceColor  = bridgecolor;
	B(2).FaceColor = curlycolor;
	title('degree')
	xlabel('number of connections')
	legend({'bridges' 'curlies'})

	% number of connections per cell 
	NC{1} = sum(sim_RM{1}.networkParameters.connectivityMatrix>0);
	NC{2} = sum(sim_RM{2}.networkParameters.connectivityMatrix>0);
	NC{3} = sum(sim_RM{3}.networkParameters.connectivityMatrix>0);

	% mean and std excluding bridge cells
	[mean(NC{1}(~bc)), std(NC{1}(~bc));
	 mean(NC{2}(~bc)), std(NC{2}(~bc));
	 mean(NC{3}(~bc)), std(NC{3}(~bc))]
	
	
	figure
	subplot(3,1,1)
	 hist(NC{1}(~bc) ,[0:15])
	 title('curlies')

	subplot(3,1,2)
	 hist(NC{2}(~bc),[0:15])
	 title('bridges+curlies')

	subplot(3,1,3)
	 hist(NC{3}(~bc),[0:15])
	 title('brick')

histopatch = @(x,y,c) patch( [x 0], [y 0], c);

figure
	subplot(3,1,1)
	 H1_nonBC = hist(NC{1}(~bc) ,[0:15]);
	 H1_BC    = hist(NC{1}      ,[0:15]);
	 title('curlies')

	subplot(3,1,2)
	 H2_nonBC = hist(NC{2}(~bc),[0:15]);
	 H2_BC    = hist(NC{2},[0:15]);
	 title('bridges+curlies')

	subplot(3,1,3)
	 H2_nonBC =hist(NC{3}(~bc),[0:15])
	 title('brick')



end



if distance_histogram
	[rr cc v_] = find(bridg_curlies.W);
	 %overall distance matrix
	 D = squareform(pdist(JM394_horizontal_coordinates(2:end,:)));

	 % distances only for connections between curlies
	 D_cur2cur = D(find(not(bc)),find(not(bc)));
	 Wcur = bridg_curlies.W(find(not(bc)), find(not(bc)));

	 % connections from bridges to bridges and curlies
	 Wbridge = bsxfun(@times, bridg_curlies.W, bc);
	 Wbrick  = brick.W;
	 
	 % indices of nonzero distances (in D  and ind D_cur2cur).
	 nonzero_cur2cur = find( (triu(D_cur2cur)-eye(size(D_cur2cur))).*Wcur);
	 nonzero_bri2cur = find( (triu(D)-eye(size(D))).*Wbridge);
	 nonzero_brick   = find( (triu(D)-eye(size(D))).*Wbrick);

	 dsup    = [0:10:500];
	 hcur    = hist( D_cur2cur(nonzero_cur2cur)  ,dsup);
	 hbridge = hist( D(nonzero_bri2cur) ,dsup);
	 hbrick  = hist( D(nonzero_brick)   ,dsup);
	 
	 figure
	 subplot(2,1,1)
	 patch([dsup 0],[hcur 0], curlycolor)
	 patch([dsup 0],[hbridge 0], bridgecolor)
	 alpha(.3)
	 title('curlies and bridges')
	 legend({'curlies' ; 'bridges'})
	 ylim([0 500])
	 subplot(2,1,2)
	 patch([dsup 0],[hbrick 0], [.7 .7 .7]);
	 alpha(.3)
	 ylim([0 500])
	 xlabel('distance \mum')
	 legend('all')
	 
end


if plotselectedclusters
	bridges_from_cluster = single(bc .* clusters==41);
	neighbors_to_bridge = find(bridges_from_cluster'*(sims{1}.W.W>0));
	their_cluster = unique(clusters(neighbors_to_bridge));
	targeted_cluster_cells = ismember(clusters, their_cluster).*clusters;
	% plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), targeted_cluster_cells);

	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), clusters==41 | clusters==34);
	view(22,56.4);
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


% [================================================]
%  cell scatters
% [================================================]


if plotcellscatters 
	% sel_fields = {'g_CaL', 'g_int', 'p1', 'g_h', 'g_ld',  'freq_each', 'meanVm','minV', 'supth'}
	% sel_fields = {'g_CaL', 'g_int', 'g_h', 'ampl', 'freq_each', 'meanVm','minV'}
	% sel_fields = {'g_CaL', 'g_ld',  'g_K_Ca', 'p1', 'g_h',  'ampl', 'freq_each', 'meanVm','minV','supth'}
	% sel_fields = {'g_CaL', 'g_ld', 'g_K_Ca',  'freq_each', 'meanVm','minV'}

	thesefields = {'ampl' 'freq_each'};
	RR = vertcat(R{1}.allneurons, R{2}.allneurons);
	groups = [zeros(1105,1) ; ones(1105,1)];
	sel_table = RR(:,sel_fields);
	NDscatter(sel_table, groups+1);

	sel_fields = { 'ampl' 'freq_each' 'g_CaL', 'g_h' };
	clusterswobridges = clusters;
	clusterswobridges(logical(bc))=0;
	sel_table = R{1}.allneurons(:,sel_fields);
	NDscatter(sel_table, clusters)

	sel_fields = { 'ampl' 'freq_each' 'g_CaL', 'g_h' };
	clusterswobridges = clusters;
	clusterswobridges(logical(bc))=0;
	sel_table = R{3}.allneurons(:,sel_fields);
	NDscatter(sel_table, clusters)

	
	% sel_table = R{2}.allneurons(:,sel_fields);
	% NDscatter(sel_table, clusters)

	% sel_table = R{3}.allneurons(:,sel_fields);
	% NDscatter(sel_table, clusters)

	for nsim = 1:length(sims);
		sel_table = R{nsim}.allneurons(:,sel_fields);
		NDscatter(sel_table, bridgecells+1);
	end
	

	
end




% [=================================================================]
%  bridge behavior
% [=================================================================]
bridgeidx = find(bridgecells);
bridge_Vm = statevar{1}(bridgeidx,:);
[val ord] = sort(table2array(R{1}.allneurons(bridgeidx,'ampl')));


if plotbridgewaterfall
	figure
	waterfall(bridge_Vm(ord,:));
	title('bridge behavior when connected, ordered by amplitude')

	figure
	bridge_Vm = statevar{2}(bridgeidx,:);
	waterfall(bridge_Vm(ord,:));
	title('bridge behavior when disconnected')

	figure
	bridge_Vm = statevar{4}(bridgeidx,:);
	waterfall(bridge_Vm(ord,:));
	title('homogeneous connectivity')

end


% [=================================================================]
%  groups
% [=================================================================]

if plotbridgeandneighbors_Vm
	selectedbridges = find(single(bc .* clusters==41 | bc .* clusters==34))';
	titleorder = {'bridges' 'only curlies' 'disconnected' 'brick'};

	for c = selectedbridges
		ind = ind+1;
		figure(c)
		thisbridge = zeros(size(bc));
		thisbridge(c) = 1;
		neighbors{c} = find(thisbridge'*(sims{1}.W.W>0));
				Q = quantile(table2array(R{3}.allneurons(neighbors{c},'freq_each')), [.1, .5, .9])

		subplot(4,1,1)
		plot(statevar{1}(neighbors{c}',:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot( statevar{1}(c,:),'linewidth', 2)
		title({['bridge number: ' num2str(c) ' degree:' num2str(length(neighbors{c}))] ; ['neighbor freq (lower median upper): ' num2str(Q) 'Hz'] ; titleorder{1}})

		subplot(4,1,2)
		plot(statevar{2}(neighbors{c},:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(statevar{2}(c,:),'linewidth', 2)
		title({['neighboring clusters:' num2str(clusters(neighbors{c})') ]; titleorder{2}})

		subplot(4,1,3)
		plot(statevar{3}(neighbors{c},:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(statevar{3}(c,:),'linewidth', 2)
		title({['neighboring clusters:' num2str(clusters(neighbors{c})') ]; titleorder{3}})

		subplot(4,1,4)
		plot(statevar{4}(neighbors{c},:)','color', [1 1 1]*.8)
		% plot(tslice, mean(statevar{1}(neighbors{c},:)),'color', [1 1 1]*.9)
		hold on
		plot(statevar{4}(c,:),'linewidth', 2)
		title({['neighboring clusters:' num2str(clusters(neighbors{c})') ]; titleorder{4}})


	end

end


if plotclustercellactivity
	selected_clusters = [34 41];
	% titleorder = {'bridges' 'only curlies' 'disconnected' 'brick'};
	for s= [1:4]
		for sc = selected_clusters
			figure
			plot(statevar{s}(find(clusters==sc),:)')
			title({['cluster:' num2str(sc)] ; titleorder{s}})

		end
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

if plotclustermemberaverages
	figure

	for sp = 1:length(statevar)
		subplot(1,length(statevar),sp)
		for c = 1:max(clusters)
			% clustered{c}.sync = measureGroupSync(sim{1},'group', clusters==c,'plotme',0);
			% clustered{c}.no_neurons = length(find(clusters==c));

			% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
			clusterset = setdiff(find(clusters==c), find(bc));

			plot(tslice, mean(statevar{sp}(clusterset,:))+c*5,'color', lc(c,:))
			text(tslice(1),c*5, ['ncells:' num2str(length(clusterset))])
			hold on
			axis tight
			% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		end
	end

end



% [=================================================================]
%  Synchrony stats for clusters
% [=================================================================]plotscatters


% calculates phase coherence of clustered cells
if calculatesynchrony_clusters
    
    %%
	fig3 = figure;;
	for c = 17 %1:max(clusters)
		c = 17;
		c = 10;
		clustered{c,1}.sync = measureGroupSync(sims{1},'group', clusters==c,'plotme',1);
		clustered{c,2}.sync = measureGroupSync(sims{2},'group', clusters==c,'plotme',1);

		clustered{c,1}.no_neurons = length(find(clusters==c));

		
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






if bridgecond_pspace
	load('bridge_conductance_pspace01-Mar-2017.mat');
	% sims = sim;
	clusters = sims{1}.W.stats.clusters
	g1g2 = find(clusters==41 | clusters==34);

	for ss  =1:11
		bridgepspace_sync{ss} = measureGlobalSync(sims{ss},'duration', [1000:6000],'group',g1g2,'plotme',0)
	end
	

	plot([	bridgepspace_sync{1}.stats.order_parameter_all...
			bridgepspace_sync{2}.stats.order_parameter_all...
			bridgepspace_sync{3}.stats.order_parameter_all...
	 		bridgepspace_sync{4}.stats.order_parameter_all...
			bridgepspace_sync{5}.stats.order_parameter_all...
	 		bridgepspace_sync{6}.stats.order_parameter_all...
	 		bridgepspace_sync{7}.stats.order_parameter_all...
	 		bridgepspace_sync{8}.stats.order_parameter_all...
	 		bridgepspace_sync{9}.stats.order_parameter_all...
	 		bridgepspace_sync{10}.stats.order_parameter_all...
	 		bridgepspace_sync{11}.stats.order_parameter_all...
	 		])

	hist([	bridgepspace_sync{1}.stats.order_parameter_all...
			bridgepspace_sync{2}.stats.order_parameter_all...
			bridgepspace_sync{3}.stats.order_parameter_all...
	 		bridgepspace_sync{4}.stats.order_parameter_all])

end



if  analyze_group_stim
	 load curlies_bridges_stim_pair22-Mar-2017.mat
	 sims = sim; clear sim
	 tslice = [500:8000];
	 g1 = find(clusters==34);
	 g2 = find(clusters==41);
	 g1g2 = find(clusters==41 | clusters==34);



	M{1} = phase_distribution_over_time(sims{1},'duration',tslice,'group',g1g2);
	% saveallfigs('prefix', 'withbridges_1clusterstim','style','12x6')
	close all
		
	M{2} = phase_distribution_over_time(sims{2},'duration',tslice,'group',g1g2);
	% saveallfigs('prefix', 'withbridges_2clusterstim','style','12x6')
	% close all


	M{3} = phase_distribution_over_time(sims{3},'duration',tslice,'group',g1g2);
	% saveallfigs('prefix', 'withoutbridges_1clusterstim','style','12x6')
	% close all
	
	M{4} = phase_distribution_over_time(sims{4},'duration',tslice,'group',g1g2);
	% saveallfigs('prefix', 'withoutbridges_2clusterstim','style','12x6')
	% close all

	plot(abs(M{2}.phases.orderparameter{1}));hold on
	plot(abs(M{4}.phases.orderparameter{1}))
	
	[hh xx] = hist([abs(M{2}.phases.orderparameter{2})'  abs(M{4}.phases.orderparameter{2})'],linspace(0,1,50));
	bar(xx,hh)
	legend({'with bridges' 'clusters only'})
	title('synchrony for whole network')

	figure
	subplot(2,1,1)
	plot(sims{2}.networkHistory.V_soma(g1,:)','r')
	hold on
	plot(sims{2}.networkHistory.V_soma(g2,:)','g')
	title('with bridges')

	subplot(2,1,2)
	plot(sims{4}.networkHistory.V_soma(g1,:)','r')
	hold on
	plot(sims{4}.networkHistory.V_soma(g2,:)','g')
	title('without bridges')


end


makemovies = 0;
if makemovies
	replayResults_clusters(sims{1},'savemovie',1,'time_slice', [2000:3500])
	replayResults_clusters(sims{2},'savemovie',1,'time_slice', [2000:3500])
	replayResults_clusters(sims{3},'savemovie',1,'time_slice', [2000:3500])
	replayResults_clusters(sims{4},'savemovie',1,'time_slice', [2000:3500])
end

makemoviesofstim = 0;
if makemoviesofstim
	M = single(sims{1}.perturbation.mask{1});
	sims{1}.perturbation.mask{1} = M;
	sims{1}.perturbation.mask{1}(find(g1)) = 1;
	sims{1}.perturbation.mask{1}(find(g2)) = 2;
	sims{2}.perturbation.mask{1} = M;
	sims{2}.perturbation.mask{1}(find(g1)) = 1;
	sims{2}.perturbation.mask{1}(find(g2)) = 2;
	sims{3}.perturbation.mask{1} = M;
	sims{3}.perturbation.mask{1}(find(g1)) = 1;
	sims{3}.perturbation.mask{1}(find(g2)) = 2;
	sims{4}.perturbation.mask{1} = M;
	sims{4}.perturbation.mask{1}(find(g1)) = 1;
	sims{4}.perturbation.mask{1}(find(g2)) = 2;
	

	M{1} = phase_distribution_over_time(sims{1},'duration',[2500:4000], 'animate', 1,'savemovie',1);
	M{2} = phase_distribution_over_time(sims{2},'duration',[2500:4000], 'animate', 1,'savemovie',1);
	M{3} = phase_distribution_over_time(sims{3},'duration',[2500:4000], 'animate', 1,'savemovie',1);
	M{4} = phase_distribution_over_time(sims{4},'duration',[2500:4000], 'animate', 1,'savemovie',1);
end

	

M{1} = phase_distribution_over_time(sims{1},'duration',[2900:3500], 'animate', 1);


if activitydifference
	load /Users/M/Synced/Titan/Bench4/curlies_bridges_randmaskstim22-Mar-2017.mat
	g1 = find(clusters==34);
	g2 = find(clusters==41);
	g1g2 = find(clusters==41 | clusters==34);

	figure
	sims{2}.W = bridg_curlies;
	sims{3}.W = bridg_curlies;
	replayResults_clusters(sims{2},'savemovie',0,'time_slice', [2000:3500])
	replayResults_clusters(sims{3},'savemovie',0,'time_slice', [2000:3500])
	actdif = [sim{2}.networkHistory.V_soma(:,3000:5000) - sim{3}.networkHistory.V_soma(:,3000:5000)];
	imagesc([repmat(sim{2}.perturbation.mask{1}*10-5,1,200) actdif])
	set(gca,'clim',[-5 5])
end


if render_volumetric_activity
	load /Users/M/Synced/Titan/Bench4/curlies_bridges_randmaskstim22-Mar-2017.mat
	tslice = 1000:5000;
	V1 = sim{1}.networkHistory.V_soma(:,tslice);
	V2 = sim{2}.networkHistory.V_soma(:,tslice);
	V3 = sim{3}.networkHistory.V_soma(:,tslice);
	V4 = sim{4}.networkHistory.V_soma(:,tslice); slice)
	plot_volume(V2, JM394_horizontal_coordinates(2:end,:),tslice)
	plot_volume(V3, JM394_horizontal_coordinates(2:end,:),tslice)
	plot_volume(V4, JM394_horizontal_coordinates(2:end,:),tslice)
end

if cluster_bridges_netstructs
	plotnetstruct(brick.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), (bridg_curlies.stats.clusters>-1)*3);
	saveallfigs('prefix', 'brick', 'style', '19x19')
	close all

	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridg_curlies.stats.clusters);
	saveallfigs('prefix', 'bri+cur', 'style', '19x19')
	close all

	plotnetstruct(curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), bridg_curlies.stats.clusters);
	saveallfigs('prefix', 'curlies', 'style', '19x19')
	close all

end

phase_activity_volume_snapshots = 1;
if phase_activity_volume_snapshots
	savemovie = 1;
	load('curlies_bridges_randmaskstim22-Mar-2017.mat')
	framestoprint = [2226:1:3850];
	mkdir sim1_randmask
	cd sim1_randmask/
	animate_volume_hilbert(sim_RM{1}, framestoprint,savemovie);
	close all
	cd ..

	mkdir sim2_randmask
	cd sim2_randmask/
	animate_volume_hilbert(sim_RM{2}, framestoprint,savemovie);
	close all
	cd ..

	mkdir sim3_randmask
	cd sim3_randmask/
	animate_volume_hilbert(sim_RM{3}, framestoprint,savemovie);
	close all
	cd ..

		framestoprint = [2226:1:3850];
	load('curlies_bridges_nostim22-Mar-2017.mat')

	mkdir sim1_nostim
	cd sim1_nostim/
	animate_volume_hilbert(sim_nostim{1}, framestoprint,savemovie);
	close all
	cd ..

	mkdir sim2_nostim
	cd sim2_nostim/
	animate_volume_hilbert(sim_nostim{2}, framestoprint,savemovie);
	close all
	cd ..

	mkdir sim3_nostim
	cd sim3_randmask/
	animate_volume_hilbert(sim_nostim{3}, framestoprint,savemovie);
	close all
	cd ..
end



if spectral_clustering
	% load /Users/M/Synced/Titan/Bench4/curlies_bridges_randmaskstim22-Mar-2017.mat
	% load /Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_bridges/clusters_curlies_bridges_22-Jan-2017.mat
	load /Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_bridges/clusters_curlies_bridges_20-Jan-2017.mat
	% load /Users/M/Synced/Projects/Experiments/Olive/model/simresults/clusters_bridges/curlies_bridges_01-Mar-2017.mat
	



	% for ss = 1:3
		ss = 3;
		PD = phase_distribution_over_time(sims{ss}, 'duration', [2000:5000]);
		ph = PD.phases.pop{2};
		x = sin(ph(:,2000:5000));
		cc = partialcorr(x', mean(x)'); % partial correlation w.r.t the mean
		rho = (cc + 1)/2; % normalize rho to get an affinity matrix 0<=rho<=1
		rho = (rho+rho')/2; % rho had better be symmetric
		c = spect_clust(rho, 50);
		[bic, aic] = baic(x, c);

	% end




	subplot(331)
	plot(x'), axis tight
	title('Data')
	subplot(332)
	pcolor(cc), shading flat
	title('Partial correlation')
	subplot(333)
	hist(cc(:),100), axis tight
	title('Partial correlation')


	subplot(334)
	plot(c.nc,c.gc,'-o'), axis tight
	xlabel('Number of clusters')
	ylabel('Eigenvalue gap')

	subplot(335)
	[~, ic] = sort(c.label(:,2));
	pcolor(cc(ic,ic)), shading flat
	title('Partial correlation')

	subplot(336)
	fshow(c.label(:,2))

	%% Check with information criteria

	subplot(337)
	[bic, aic] = baic(x, c);
	plot([bic aic]), axis tight % These criterion give different results.. 
	legend('BIC','AIC')
	subplot(338)
	[~,ic] = sort(c.label(:,20)); % Let's go for n = 10 as BIC suggests..
	pcolor(cc(ic,ic)), shading flat % For some reasons, it doesn't seem to work..
	subplot(339)
	fshow(c.label(:,10))
end


randstim_carpet_diffs = 1;
if randstim_carpet_diffs
	clear
	load curlies_bridges_randmaskstim22-Mar-2017.mat
	tslice = [1:10000];
	for ii = 1:3
		H{ii} = hilbert_of_membranepotential(sim_RM{ii}.networkHistory.V_soma(:,tslice));
	end

	[v ord] = sort(clusters);
	oord = ord(~bc);

	load curlies_bridges_nostim22-Mar-2017

	for ii = 1:3
		H{ii+3} = hilbert_of_membranepotential(sim_nostim{ii}.networkHistory.V_soma(:,tslice));
	end


	for s = 1:6
		for g = 1:50;
			GS{s}(g,:) = abs(mean(H{s}.hilbert(clusters==g,:))) ;
		end
	end



	f100 = figure;
			
	oord = ord(~bc(ord))
	ax(1) = subplot(2,3,1);
		imagesc(H{1}.hilbert(oord,:))
		ax(2) = subplot(2,3,2);
		imagesc(H{2}.hilbert(oord,:))
		ax(3) = subplot(2,3,3);
		imagesc(H{3}.hilbert(oord,:))

		
		
		ax(4) = subplot(2,3,4);
		imagesc(H{4}.hilbert(oord,:))
		ax(5) = subplot(2,3,5);
		imagesc(H{5}.hilbert(oord,:))
		ax(6) = subplot(2,3,6);
		imagesc(H{6}.hilbert(oord,:))


	f200 = figure;
	subplot(131), imagesc([H{4}.hilbert(oord,:) - H{1}.hilbert(oord,:)]')	
	subplot(132), imagesc([H{5}.hilbert(oord,:) - H{2}.hilbert(oord,:)]')	
	subplot(133), imagesc([H{6}.hilbert(oord,:) - H{3}.hilbert(oord,:)]')
[
end


measuregropusync = 1;
if measuregropusync

	for s = 1:3

		for g = 1:50;
			groupsync{g} = measureGroupSync(sim_nostim{s}, 'duration',[2500:6000], 'plotme', 0, 'group', find(sim_nostim{s}.W.stats.clusters==g));
			G(s,g) = groupsync{g}.sync(1);
			V(s,g) = groupsync{g}.sync(2);
		end
	end

end




	f102 = figure;

		ax(1) = subplot(3,1,1);
		imagesc(abs(GS{1}),[0 1])
		ax(2) = subplot(3,1,2);
		imagesc(abs(GS{2}),[0 1])
		ax(3) = subplot(3,1,3);
		imagesc(abs(GS{3}),[0 1])

