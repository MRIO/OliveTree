function varargout = replayResults_clusters(varargin)

% #TODO: PROFILE_SIM AND HISTOGRAMS

% >> C = simR{2}.W.stats.clusters;
% >> hist(table2array(R.allneurons(C==1,'freq_each')),30)

% p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
% p.addRequired('time_slice')
% p.addRequired('savemovie')
% p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
% p.addParamValue('fhandle', gcf)

%TODO: accept any dt

static = 1;
plotmeanclusteractivity = 1;
plotpopactivity = 1;
calculatesynchrony = 0;
plotthreeDscatter = 0;
reconstruction = 1;

trigger = 1;


% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('time_slice',[])
p.addOptional('savemovie',0)
p.addOptional('snap_time',1)
p.addParameter('plotallfields', 0) % stdandard deviation criterion for offset threshold
p.addParameter('fhandle', gcf)
p.addParameter('plotspikes', 0)

p.parse(varargin{:});

sim = p.Results.sim;
time_slice = p.Results.time_slice;
savemovie = p.Results.savemovie;
snap_time = p.Results.snap_time
plotspikes = p.Results.plotspikes;

plotallfields  = p.Results.plotallfields;
fhandle = p.Results.fhandle;


% [=================================================================]
%  retrieve fields
% [=================================================================]
frames = [1:5000];
netsize 			= sim.networksize;
simtime             = length(sim.networkHistory.V_soma);
noneurons           = prod(netsize);
V_soma_unwrapped 	= sim.networkHistory.V_soma;
% V_soma = reshape(V_soma_unwrapped(:,1) , noneurons, 1);
pert 				= sim.perturbation;
loggedstatevars = fieldnames(sim.networkHistory)';

coords = sim.W.coords;

% if isfield(sim, 'stats.clusters')


% [=================================================================]
%  slice to plot
% [=================================================================]
 
 if ~isempty(time_slice)
 	simtime = length(time_slice);
 else
 	time_slice = [1:simtime];
 end
 


% [=================================================================]
%  % retrieve clusters and order ascendingly
% [=================================================================]

if isfield(sim.W.stats, 'clusters')
	clusters = sim.W.stats.clusters;
	[orderedclusters O] = sort(clusters);
	V_soma_ordered = V_soma_unwrapped(O,:);
	no_clusters = length(unique(clusters));
else
	disp('did not find clusters in input struct.')
	return
	O = [1:prod(netsize)];
end


if ~exist('cbrewer')
		lc = jet(no_clusters);
	else
		lc = cbrewer('qual', 'Set1', no_clusters);
end



%crop
V_soma_unwrapped = V_soma_unwrapped(:,time_slice);
V_soma_ordered   = V_soma_ordered(:,time_slice);
% [=================================================================]
%  spike detection
% [=================================================================]
sim.networkHistory.V_soma = V_soma_unwrapped;
spks = spikedetect(sim)

propspkneurons = spks.propspkneurons;
popfreq = spks.popfrequency;
f = sprintf('%.2f', popfreq);
edges = [-1 0.01:.5:max(spks.medfreq)];
histfreq = histc(spks.medfreq, edges);



% [=================================================================]
%  spike triggered count of stimulus
% [=================================================================]
% #todo

% [=================================================================]
%  spike triggered average of membrane potential
% [=================================================================]
% #todo

% [=================================================================]
%  spike triggered average of neighbor's membrane potential
% [=================================================================]
% #todo


% [=================================================================]
%  if saving movie
% [=================================================================]

if savemovie
	
	if isfield(sim,'note')
		fname = sim.note;
	else
		fname = ['sim' num2str(netsize) '_']
	end
		try
		vidObj = VideoWriter(fname,'MPEG-4');
		
		catch
			vidObj = VideoWriter(fname);
		end
		set(vidObj,'FrameRate', 24)

	open(vidObj);
end




% [=================================================================]
%  prepare figure
% [=================================================================]

if plotpopactivity
	fig1 = gcf;
	set(fig1, 'position', [1 1 1024 768],'color', [1 1 1]);	
	colormap(hot);


	% if isfield(sim, 'clusters')

	 
	 if ~isempty(time_slice)
	 	simtime = length(time_slice);
	 else
	 	time_slice = [1:simtime];
	 end
	 
	fig1 = figure('position', [  1           1        1440         900]);

	colormap(hot);

	%    _   __    __    _      __                                       
	%   (_)_/_/   / /_  (_)____/ /_____  ____ __________ _____ ___  _____
	%    _/_/    / __ \/ / ___/ __/ __ \/ __ `/ ___/ __ `/ __ `__ \/ ___/
	%  _/_/_    / / / / (__  ) /_/ /_/ / /_/ / /  / /_/ / / / / / (__  ) 
	% /_/ (_)  /_/ /_/_/____/\__/\____/\__, /_/   \__,_/_/ /_/ /_/____/  
	%                                 /____/                             

	% 1. correlation between firing and Calcium neighborhood
	% 2. correlation between gap neighborhood and phase distribution 


	% raster / vsomma
	   
	a(1) = axes('position', [0.12    0.07    0.85    0.6]);

		imagesc(V_soma_ordered,[-68 -30]); %imagesc(V_soma_unwrapped',[-68 -30]);

        
        
		% set(gca,'xtick',[1 noneurons]);
		hold on;
	    xlabel('ms');
	    ylabel('neurons');
        axis tight
        
		plotspikes1 = @(c)plot(spks.spikes{O(c)},O(c),'markersize', 3,'marker', 'o','linestyle', 'none','color', 'w');
		plotspikes2 = @(c)plot(spks.spikes{O(c)},O(c),'markersize', 10,'marker', '.','linestyle', 'none','color', 'g');
		
		for c = 1:prod(netsize)

			if not(isempty(spks.spikes{O(c)})) & plotspikes
				plotspikes1(c);
				plotspikes2(c);
			end
		end

% 		colorbar
		title('Vm')

		if isfield(sim.perturbation,'mask')
			M = sim.perturbation.mask{1}(O)*20-100;
			b(1) = axes('position', [0.07 0.07 0.05 0.6])
			imagesc(M)
		end

		if isfield(sim.W.stats, 'clusters')
			C = sim.W.stats.clusters;
			b(2) = axes('position', [0.02 0.07 0.05 0.6])
			imagesc(orderedclusters)
		end

	% spike histogram

	if sum(cell2mat(spks.spikes))
		a(2) = axes('position', [0.12    0.67    0.85    0.3]);
		[hh x] = hist(cell2mat(spks.spikes),[1:simtime]);
		K = conv(hh, gausswin(50), 'same')/sum(gausswin(50));
		K = conv(hh, ones(1,10), 'same')/sum(ones(1,10));
		xlim([time_slice(1) time_slice(end)])
		bar(x,K,'facecolor','k')
		set(a(2),'xtick',[])
	end


	if isfield(sim.perturbation,'triggers')
		% stimulus
		hold on
		a(3) = axes('position', [0.07    0.67    0.85    0.3]);
		[hh x] = hist(sim.perturbation.triggers{trigger},simtime);
		bar(x,hh,'facecolor','g')

		set(a(3),'xtick',[],'ytick',[],'color','none')
		axis off
	end


	linkaxes(a, 'x')
    linkaxes([a(1) a(3)], 'y')
    % xlim([time_slice(1) time_slice(end)])
end

% figure
% for c = 1:no_clusters

% 	MMM(c,:) = mean(V_soma_unwrapped(find(V==c),:));
% end
% keyboard
% waterfall(MMM)
% 	xlabel('time (ms)')
% 	ylabel('mV')








	if 0
		plotvolume = 1;

		if plotvolume
			
			gksz = 7;
			g3d3 = gaussKernel3d(.2,  gksz, ceil(gksz/2)); g3d3 = g3d3/sum(g3d3(:)); g3d3(g3d3<.005) = 0; g3d3(g3d3>0.005)=1; %g3d3 = g3d3/sum(g3d3(:)); 
			% figure
			% cla
			% vol3d('cdata',g3d3)
			% colorbar

			coarseness = 4;


			fig_volume = figure('color', [1 1 1]);
			ax_volume = axes;
			colormap(ax_volume, jet(32));
			colorbar

			for tt = 1:length(time_slice)
				VVVV = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), V_soma_unwrapped(:,tt));
				NNNN = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), 1);
				NNNN(NNNN==0)=1;
				VVVV = VVVV./NNNN;
				
				% CCCC = interp3(VVVV, 3);
				% CCCC = convn(VVVV, g3d3, 'same');
				% VVVV = CCCC;
				CCCC = imerode(VVVV,g3d3);

				set(0,'CurrentFigure',fig_volume);
			    % set(fig1,'CurrentAxes',a(3));


			    AAA = CCCC;
				AAA(CCCC>-35) = 1;
				AAA(AAA<=-35) = 0;

				% VVVV(VVVV>-50) = -45;
				% VVVV(VVVV<-67) = -67;
				% VVVV(1,1,1) = -67; 
				% VVVV(1,2,1) = -50;

				cla(ax_volume)
				set(ax_volume, 'clim',[-60 -40])
				
				vol3d('cdata',CCCC, 'Alpha', ~AAA*.25 , 'texture','3D');


				if tt==1;view(3) ; view(-22,-56.4);axis off; axis tight;  daspect([1 1 1]);end
				title([num2str(time_slice(tt)) 'ms'])
				drawnow
				if savemovie
					writeVideo(vidObj, getframe(fig_volume))
				end
				

			end
			if savemovie
				close(vidObj)
			end

		end
	end




if plotthreeDscatter
	fig2 = figure('colormap',lc);

	scaledV = (V_soma_unwrapped -min(min(V_soma_unwrapped)) +1      )*20;

	if reconstruction

		for tt = 1:size(V_soma_unwrapped,2)
			cla
			scatter3(coords(:,1), coords(:,2), coords(:,3), scaledV(:,tt),clusters ,'filled')
			% caxis([-80,-30])
			title(num2str(tt))
			axis equal
			colorbar
			drawnow

		end
	end


end


% figure
% neuronselection = [1:10:noneurons];
% set(0,'defaultaxescolororder', linspecer(2))
% plot(V_soma_unwrapped(neuronselection,:)','linewidth',2);

if plotmeanclusteractivity
	fig03 = figure;;
	disp('calculating cluster synchrony.')
	for c = 1:no_clusters
		try
			plot(1:simtime, mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
			% plot_mean_and_std(V_soma_unwrapped(find(clusters==c),:)+c*5,'color', lc(c,:))
		catch 
			continue
		end

		
		hold on

		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
    end
    title('mean cluster activity')
    xlabel('ms')
end



% [=================================================================]
%  Cluster Activity Stats
% [=================================================================]
if calculatesynchrony
	fig3 = figure;;
	disp('calculating cluster synchrony.')
	for c = 1:no_clusters

		c
		clustered{c}.sync = measureGroupSync(sim,'group', clusters==c,'plotme',0);
		clustered{c}.no_neurons = length(find(clusters==c));

		% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
		% plot([1:simtime], mean(V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
		% hold on
		% plot([1:simtime], V_soma_unwrapped(find(clusters==c),:))+c*5,'color', lc(c,:))
	end
		% xlabel('time (ms)')
		% ylabel('mV')
end






