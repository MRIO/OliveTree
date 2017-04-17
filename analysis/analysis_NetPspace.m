% ANALYSIS noise_corr_gap_comparisons.m

% [=================================================================]
%  script parameters
% [=================================================================]
% pathtofile  = '/Users/M/Synced/Projects/Experiments/Olive/model/simresults/netpspace/';
% pathtofile  = '/Users/M/Projects/Experiments/Olive/model/simresults/netpspace/';
% pathtofile  = './';


if ~exist('transients') 
	% conntype = 'iso'; load /Users/M/SyncBox/Titan/Bench/pspace_noise_gap_24_iso_10000_17-May-2016.mat
	% load /Users/M/SyncBox/Titan/Bench/pspace_noise_gap_24_cluster_10000_18-May-2016
	% load /Users/M/Synced/Titan/Bench2/netpspace24_iso_5000_27-May-2016.mat
	% load /Users/M/SyncBox/Titan/Bench/netpspace24_iso_30000_28-May-2016.mat
	% fname = '/Users/M/SyncBox/Titan/Bench/netpspace30_iso_5000_03-Jun-2016.mat'; % pspace for gap compensation
	
	% fname =  '/Users/M/Synced/Titan/Bench/netpspace/netpspace_input_108_cluster_10000_02-Jun-2016.mat'; %pspace for noise amplitude and sigma
	% dimnames = {'noiseamp', 'noisesig', 'gap' };
	
	% fname = '/Users/M/Synced/Titan/Bench/netpspace/netpspace16_iso_3000_07-Jun-2016.mat';
	% dimnames = {'sametoall', 'gapcomp', 'gap' };
	
	% fname = '/Users/M/Synced/Titan/Bench/netpspace/netpspace16_iso_3000_07-Jun-2016_2.mat';
	% dimnames = {'sametoall', 'tau', 'gap' };

	% fname = '/Users/M/Synced/Titan/Bench/clustersize_netpspace40_cluster_5000_11-Jun-2016.mat'; 
	% dimnames = {'sametoall', 'clustersize', 'gap' };

	% fname = '/Users/M/Synced/Titan/Bench/netpspace/tauXeta_netpspace84_iso_3000_20-Jun-2016.mat';
	% dimnames = {'sametoall', 'tau', 'gap' };

	% fname = 'sigXmu_0corr_netpspace98_iso_3000_28-Jun-2016.mat'
	% dimnames = {'noisesig', 'noiseamp', 'gap' };

	% fname = 'RADIUS_X_NCORR_netpspace25_iso_3000_28-Jun-2016.mat'
	% dimnames = { 'radius', 'corr',	 'gap' };

	% fname = 'sigXmu_01corr_R2_netpspace98_iso_3000_25-Sep-2016.mat'
	% fname =	'sigXmu_01corr_R2_netpspace98_iso_5000_25-Feb-2017.mat';
	% dimnames = { 'noisesig', 'noiseamp', 'gap' };	

	fname = 'gapcomp_corr_R2_g04_netpspace24_iso_6000_06-Apr-2017.mat';
	dimnames = { 'gap', 'gapcomp', 'noisesig' };	

	
	[pth fnm] = fileparts(fname)

	% load([pathtofile fname])
	load([fname])
end


if ~exist('PTable')
	try
		load([fname(1:end-4) 'results.mat']);
		preparetables = 0;
	catch
		preparetables = 1
	end
else 
	preparetables = 0
end

plottraces = 0;
	% sims2plot = [];
	sims2plot = 1:length(transients);
	tslice = [1000:6000];


summarize = 1;
	plotsurfaces = 1;
		% measured = {'freq', 'pop_r', 'prop_f', 'FO_sync', 'SO_sync', 'All_sync','meanVm','ampl'};
		measured = {'freq', 'pop_r', 'prop_f', 'FO_sync', 'SO_sync', 'All_sync','meanVm','ampl'};
		% measured = {'pop_r', 'ampl'};

	save2svg = 0;
	save2png = 0;

computeselectedxcorr = 0;
	nwins = 1;
	plot_selected_neurons = 1;
	sims2p = [6 18];
	% sims2p = find(PTable.tau==5  & (PTable.sametoall== 0 | PTable.sametoall==.3) & PTable.meannoconn==8)';

	plotxcorrs = 1;

multiwindowxcorr = 0;

rndseed = 0;
rng(rndseed,'twister')


set(0,'defaultaxescolororder', linspecer(5));
% set(0,'defaultfigurecolormap', linspecer(length(transients{1}.Plist)));
set(0,'defaultfigurecolor', [1 1 1]);

CB = @(n) flipud(cbrewer('div', 'RdBu', n));
set(0, 'defaultfigurecolormap', CB(20))


% [=================================================================]
%  simulation parameters
% [=================================================================]

% delta = .02;
% Fs = 1000; % sampling rate =! delta

%=============================Make Parameter Table ==============================%


if preparetables 



	ci = [0.05 .5 .95];
	qt = [0.25 .5 .75];

	sims			= [1:size(Plist,1)]';
	tau 			= Plist(:,1);
	sametoall 		= Plist(:,2);
	noiseamp  		= Plist(:,3);
	gap 			= Plist(:,4);
	radius 			= Plist(:,5);
	clustersize 	= Plist(:,6);
	intraclusterP 	= Plist(:,7);
	extraclusterP	= Plist(:,8);
	meannoconn 		= Plist(:,9);
	numneurons 		= Plist(:,10);
	depth 			= Plist(:,11);
	gapcomp 		= Plist(:,12);
	noisesig 		= Plist(:,13);

	PTable = table(sims, tau, sametoall, noiseamp, noisesig, gap, gapcomp, radius, clustersize, intraclusterP, extraclusterP, meannoconn, numneurons);


	 for s = 1:size(Plist,1)
	 	s
				

		T = transients{s};
		V = T.networkHistory.V_soma;
		CAL = T.cellParameters.g_CaL;
		IH  = T.cellParameters.g_h;

		ampl_temp 		= quantile(max(V, [], 2) - min(V, [], 2), qt);
		ampl(s,:) 		= ampl_temp(2);
		meanVm(s,:) 	= quantile(mean(V, 2), ci);
		caL(s,:) 		= quantile(CAL, ci);

		if isfield(T.W.stats, 'clusters')
			firstcluster = find(T.W.stats.clusters==1);
			K 		= measureGlobalSync(transients{s}, 'duration', tslice, 'plotme', 0, 'group',firstcluster); %, 'duration', tslice
		else
			K 		= measureGlobalSync(transients{s}, 'duration', tslice, 'plotme', 0); %, 'duration', tslice
		end
		
		spikes  = spikedetect(transients{s});

		FO_sync (s, :)  = K.stats.firstordersync(1,:);
		SO_sync (s, :)  = K.stats.secondordersync(1,:);
		All_sync(s, :)  = K.stats.overallsync(1,:);

		try
			freq(s, :) = median(K.instantaneousFrequency(:));
		catch
			freq(s,:) = -1;
		end

		pop_r(s, 1)  = spikes.popfrequency;
		prop_f(s, 1) = spikes.propspkneurons;

		W = transients{s}.networkParameters.connectivityMatrix;
			W(W==0) = [];
			truegaps(s,:) = quantile(W, ci) ;

		W = transients{s}.networkParameters.connectivityMatrix;
			truenconn(s,:) = quantile(sum(W>0), ci) ;

	end

	RTable = table(truegaps, truenconn, pop_r, prop_f,  ampl, meanVm, caL, freq,  FO_sync, SO_sync, All_sync);


	PTable.Properties.UserData = X_README;
	save(['tables' fname], 'PTable', 'RTable')

end



% [=================================================================]
%  plot trace matrix
% [=================================================================]

if plottraces
	if ~exist('tslice');tslice = [1:simtime];end
	
	
		pind = 0;
		for sims = sims2plot

			thissim = 		num2str(sims);
			thistau = 		table2array(PTable(sims,'tau') );
			thiscorr = 		table2array(PTable(sims,'sametoall') );
			thisgap  = 		table2array(PTable(sims,'gap')       );
			thisnamp = 		table2array(PTable(sims,'noiseamp')  );
			thisnsig = 		table2array(PTable(sims,'noisesig')  );	
			thismeanconn =  table2array(PTable(sims,'meannoconn'));
			thisCaL =  		table2array(RTable(sims,'caL'));
			% thisIh  =  		table2array(RTable(sims,'Ih'));

			ParamVals = ...
				{['SIM:   ' num2str(thissim)]
				 ['tau:   ' num2str(thistau)];
				 ['gaps:   ' num2str(thisgap)];
				 ['ncorr:  ' num2str(thiscorr)];
				 ['nsig:   ' num2str(thisnsig)];
				 ['namp:   ' num2str(thisnamp)];
				 ['meanconn:' num2str(thismeanconn)];
				 ['CaL:' num2str(thisCaL)];
				 % ['Ih:'  num2str(thisIh)];
				 }

			thisstofreq   = median(table2array(RTable(sims,'freq')));
			thispopfreq   = table2array(RTable(sims,'pop_r'));
			thisFOsync    = table2array(RTable(sims,'FO_sync'));
			thisSOsync    = table2array(RTable(sims,'SO_sync'));

			ResVals= ...
			 {['sto freq: ' num2str(thisstofreq) ' Hz'];
			 ['pop freq: ' num2str(thispopfreq)  ' Hz'];
			 ['FO sync : ' num2str(thisFOsync )  'mu, std'];
			 ['SO sync : ' num2str(thisSOsync )  'mu, std']}
	
			% if thisnamp == inputamplitude
				pind = pind+1;

				vs = transients{sims}.networkHistory.V_soma(:,tslice);

				f = figure;

					ax(1) = subplot(2,3,[1 2])
					set(0, 'defaultaxescolororder', bone(200))
					% ax1(pind) = subplot(length(noisecorr), length(gaps), pind);
					plot(tslice, vs');
					ylabel('mV')
					xlabel('ms')
					axis tight
					ylim([-80 10])
					
					ax(2) = subplot(2,3,[4 5])
					% ax2(pind) = subplot(length(noisecorr), length(gaps), pind);
					imagesc(tslice, [1:noneurons], vs, [-80 -20])
					colormap(bone)
					set(gca,'clim',[-65 -20])
					% set(gca,'ytick', [1:56],'yticklabel', num2str(transients{ind}.Plist),'fontsize',8)


					
					linkaxes(ax,'x')

					ax(3) = subplot(2,3,3 )
					t2 = text(0,1,ResVals );
					t2.VerticalAlignment = 'top';
					t2.BackgroundColor = [1 1 1];
					t2.EdgeColor = [0 0 1];
					axis off

					ax(4) = subplot(2,3,6 )
					t2 = text(0,1,ParamVals);
					t2.VerticalAlignment = 'top';
					t2.BackgroundColor = [1 1 1];
					t2.EdgeColor = [0 0 1];
					axis off

					
					% maximize_fig
					f.Position = [0,0,560,480];

					drawnow
					if save2png
						prefix = ['netpspace_noise' num2str(pind)];
				   		fname = [prefix '_' num2str(i) '.png'];
				   		snam='12x12';
				   		s=hgexport('readstyle',snam);
					    s.Format = 'png';
					    hgexport(f,fname,s);


					% export_fig(['netpspace_noise' num2str(pind)],'-png',f,'-m2')
						% export_fig(,'-png','-r300', f)
					end

				close

			% end

	end

end




% [=================================================================]
%  parameter space grid
% [=================================================================]



if summarize
	
		% rows = PTable.tau==5  & (PTable.sametoall== 0 | PTable.sametoall==.3) & PTable.meannoconn==8;

	CTable = horzcat(PTable, RTable);
	
	cols = {'noiseamp', 'noisesig', 'sametoall', 'gapcomp' , 'clustersize', 'meannoconn', 'radius' , 'tau', 'gap',  'freq', 'pop_r', 'prop_f', 'FO_sync', 'SO_sync', 'All_sync' , 'meanVm', 'ampl'};

	

	STable = CTable(: , cols);
	 % STable.freq = median(STable.freq')';


	for m = measured

		% warning off
		R = table2grid(dimnames, m, STable);
		% warning on

		[XX YY] = meshgrid(R.ticks.X, R.ticks.Y);
		Rsz = size(R.results);
		if prod(size(Rsz))==2; Rsz(3) = 1; end;

		if plotsurfaces

			f = figure;
			nlevels = 5;
			set(f,'colormap', CB(nlevels))

			for z =1:Rsz(3)
				

				axx(z) = subplot(Rsz(3),1,z);

				[c ,h] = contourf(XX,YY,R.results(:,:,z)');

				% clabel(c,'color','w','fontweight','bold')
				 clabel([],h,'LabelSpacing',72)
				 xlabel(dimnames{1});
				 ylabel(dimnames{2});
				 title ([m ' ' R.ticks.Z(z)],'interpreter', 'none');
				 % axis equal
				 nlevels = max([nlevels length(h.LevelList)]);

			end
			try
			set(axx(:),'clim', [floor(min(R.results(:))) max(R.results(:))] )
			catch
			end

				 if save2svg
				 	plot2svg([fnm '_' m{1} '.svg'], f)
				 end

				 if save2png
				 	% maximize_fig
				 	set(f,'position', [0 0 425 800])
				 	export_fig([fname '_' m{1} '.png'],'-r300')
				 end

		end
	end

end


if computeselectedxcorr
% [12 16 19 24]
	

	c = 0;
	for simcount = sims2p
		c = c+1;
		

		thisW  = transients{simcount}.networkParameters.connectivityMatrix;
		
		if ~exist('selectedneurons')
			selectedneurons = [];
		else
			selectedneurons = X{1}.selectedneurons;
		end

		X{c} = xcorr_summa(transients{simcount},'nwins',nwins,'plotme',plotxcorrs ,'selectedneurons', selectedneurons );


		if plot_selected_neurons
	        depth   = [1:netsize(1)];
	        breadth = [1:netsize(2)];
	        height  = [1:netsize(3)];

	        % # compute adjacency matrix
	        [XX YY ZZ] = meshgrid(depth,breadth,height);
	        XX = XX(:); YY = YY(:); ZZ = ZZ(:);

	        idx = zeros(length(XX),1);
	        idx(X{c}.selectedneurons) = 1;

	        plotnetstruct(thisW,XX,YY,ZZ,idx)
	    end

	end

	
	for c = 1:length(sims2p)
		collectX(:,c) = X{c}.XcorrNoAc;
		collectAmpl(:,c) = X{c}.amplitude;
		collectDelay(:,c) = X{c}.delay;
		collectAsym(:,c) = X{c}.asymmetry;
		% hold on
	end


	figure
	plot([-400:400], collectX)
	legend(num2str(sims2p'))
	% axis tight
	axis normal
	xlabel('ms')
	ylabel('correlation (coeff)')
	xlabel('aggregate correlation')
	

try
figure; boxplot([collectAsym],'notch', 'on'); title('asym')
figure; boxplot([collectDelay],'notch', 'on'); title('delay')
figure; boxplot([collectAmpl],'notch', 'on'); title('amplitude')
catch

end

end

