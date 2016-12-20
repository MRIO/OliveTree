function replayResults_clusters(varargin)

% #TODO: PROFILE_SIM AND HISTOGRAMS

% >> C = simR{2}.W.stats.clusters;
% >> hist(table2array(R.allneurons(C==1,'freq_each')),30)

% p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
% p.addRequired('time_slice')
% p.addRequired('savemovie')
% p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
% p.addParamValue('fhandle', gcf)

%TODO: accept any dt

static = 0;

% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('time_slice',[])
p.addOptional('savemovie',0)
p.addOptional('snap_time',1)

p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
p.addParamValue('fhandle', gcf)

p.parse(varargin{:});

sim = p.Results.sim;
time_slice = p.Results.time_slice;
savemovie = p.Results.savemovie;
snap_time = p.Results.snap_time

plotallfields  = p.Results.plotallfields;
fhandle = p.Results.fhandle;


% [=================================================================]
%  retrieve fields
% [=================================================================]

netsize 			= sim.networksize;
simtime             = length(sim.networkHistory.V_soma);
noneurons           = prod(netsize);
V_soma_unwrapped 	= sim.networkHistory.V_soma;
% V_soma = reshape(V_soma_unwrapped(:,1) , noneurons, 1);
pert 				= sim.perturbation;
loggedstatevars = fieldnames(sim.networkHistory)';

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
%  % order to present the neurons
% [=================================================================]

if isfield(sim.W.stats, 'clusters')
	[V O] = sort(sim.W.stats.clusters);
	V_soma_unwrapped = V_soma_unwrapped(O,:);
	no_clusters = length(unique(V));
else
	disp('did not find clusters in input struct.')
	return
	O = [1:prod(netsize)];
end

%crop
V_soma_unwrapped = V_soma_unwrapped(:,time_slice);
% [=================================================================]
%  spike detection
% [=================================================================]
sim.networkHistory.V_soma = V_soma_unwrapped;
spks = spikedetect(sim );

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
	fname = ['sim' num2str(netsize) '_']
		
		vidObj = VideoWriter(fname,'MPEG-4');
		set(vidObj,'FrameRate', 100)

	open(vidObj);
end




% [=================================================================]
%  prepare figure
% [=================================================================]

fig1 = gcf;
set(fig1, 'position', [1 1 1024 768],'color', [1 1 1]);	
colormap(hot);


static = 0;
trigger = 2;



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
   
a(1) = axes('position', [0.07    0.07    0.85    0.6]);

	imagesc(V_soma_unwrapped,[-68 -30]); %imagesc(V_soma_unwrapped',[-68 -30]);
	% set(gca,'xtick',[1 noneurons]);
	hold on;
    xlabel('ms');
    ylabel('neurons');  

	plotspikes1 = @(c)plot(spks.spikes{c},c,'markersize', 3,'marker', 'o','linestyle', 'none','color', 'w','markersize',6);
	plotspikes2 = @(c)plot(spks.spikes{c},c,'markersize', 10,'marker', '.','linestyle', 'none','color', 'g','markersize',10);
	for c = 1:prod(netsize)
		if not(isempty(spks.spikes{c}))
			plotspikes1(c);
			plotspikes2(c);
		end
	end

	colorbar('eastoutside');
	title('Vm')

	if isfield(sim.perturbation,'mask')
		M = sim.perturbation.mask{1}(O)*20-100;
		b(1) = axes('position', [0.07 0.07 0.05 0.6])
		imagesc(M)

	end



% spike histogram

if sum(cell2mat(spks.spikes))
	a(2) = axes('position', [0.07    0.67    0.85    0.3]);
	[hh x] = hist(cell2mat(spks.spikes),simtime);
	K = conv(hh, gausswin(50), 'same')/sum(gausswin(50));
	bar(x,K,'facecolor','k')

	set(a(2),'xtick',[])
end


if isfield(sim.perturbation,'all_pulses')
	% stimulus
	hold on
	a(3) = axes('position', [0.07    0.67    0.85    0.3]);
	[hh x] = hist(sim.perturbation.all_pulses(:),simtime);
	bar(x,hh,'facecolor','g')

	set(a(3),'xtick',[],'ytick',[],'color','none')
	axis off
end
linkaxes(a, 'x')


% figure
% for c = 1:no_clusters

% 	MMM(c,:) = mean(V_soma_unwrapped(find(V==c),:));
% end
% keyboard
% waterfall(MMM)
% 	xlabel('time (ms)')
% 	ylabel('mV')



if ~exist('cbrewer')
		lc = jet(no_clusters);
	else
		lc = cbrewer('qual', 'Set1', no_clusters)
end

figure
for c = 1:no_clusters
c	
	% plot_mean_and_std([1:simtime], V_soma_unwrapped(find(V==c),:),'color', lc(c,:))
	plot([1:simtime], mean(V_soma_unwrapped(find(V==c),:))+c*5,'color', lc(c,:))
	hold on
end
	xlabel('time (ms)')
	ylabel('mV')




threeD = 1;
if threeD
		
	plotvolume = 1;
	if plotvolume
	fig2 = figure;
	[xx yy zz]  = meshgrid(1:netsize(1),1:netsize(2),1:netsize(3));
	for tt = snap_time

		V = reshape(V_soma_unwrapped(:,tt), [netsize(1) netsize(2) netsize(3)]);

		set(0,'CurrentFigure',fig2);
	    % set(fig1,'CurrentAxes',a(3));


		V(V>-30) = -20;
		V(1,1,1) = min(min(min(V))); Vq(1,1,2) = -20;

		cla
		
		vol3d('cdata',V,'texture','3D');
		view(3)
		axis off; axis tight;  daspect([1 1 1])
		pause(.1)
		drawnow

	end

end


% figure
% neuronselection = [1:10:noneurons];
% set(0,'defaultaxescolororder', linspecer(2))
% plot(V_soma_unwrapped(neuronselection,:)','linewidth',2);

end





