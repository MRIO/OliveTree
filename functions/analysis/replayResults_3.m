function [spks] = replayResults_3(varargin)
% p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
% p.addRequired('time_slice')
% p.addRequired('savemovie')
% p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
% p.addParamValue('fhandle', gcf)
% TODO: accept any dt
% TODO: time_slice _must be_ the x axis

static = 0;

meanlinewidth = 2;

% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('time_slice',[])
p.addOptional('savemovie',0)

p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
p.addParamValue('fhandle', gcf)
p.addParamValue('plot_info', 1)

p.parse(varargin{:});

sim = p.Results.sim;
time_slice = p.Results.time_slice;
savemovie = p.Results.savemovie;

plotallfields  = p.Results.plotallfields;
fhandle = p.Results.fhandle;
plot_info = p.Results.plot_info;


% [=================================================================]
%  retrieve fields
% [=================================================================]

netsize 			= sim.networksize;
simtime             = length(sim.networkHistory.V_soma);
noneurons           = prod(netsize);
pert 				= sim.perturbation;
loggedstatevars = fieldnames(sim.networkHistory)';


% [=================================================================]
%  slice to plot
% [=================================================================]
 
 if isempty(time_slice)
 	simtime = sim.duration;
 	time_slice = [1:simtime];
 
 end
 
V_soma_unwrapped = sim.networkHistory.V_soma;


% [=================================================================]
%  spike detection
% [=================================================================]

spks = spikedetect(sim,'time_slice',time_slice);




propspkneurons = spks.propspkneurons;
popfreq = spks.popfrequency;
f = sprintf('%.2f', popfreq);
edges = [-1 0.01:.5:max(spks.medfreq)];
histfreq = histc(spks.medfreq, edges);

V_soma = reshape(V_soma_unwrapped(:,1) ,prod(netsize), 1);



% order to present the neurons


% % [=================================================================]
% %  if saving movie
% % [=================================================================]

% if savemovie
% 	fname = ['sim' num2str(netsize) '_']
		
% 		vidObj = VideoWriter(fname,'MPEG-4');
% 		set(vidObj,'FrameRate', 100)

% 	open(vidObj);
% end


% [=================================================================]
%  prepare figure
% [=================================================================]


fig1 = gcf;
set(fig1, 'position', [1 1 1024 768],'color', [1 1 1]);	
try
	cmap = flipud(cbrewer('seq', 'Greys', 20));
	% cmap = flipud(cbrewer('div', 'RdBu', 30));
	colormap(cmap);	
	clf
catch
	disp(['could not find color brewer, using hot'])
end




% [=================================================================]
%  plot populations according to type of perturbation
% [=================================================================]

a(4) = axes('Position',[0.1    0.0781    0.80    0.4]);

 
nonstim_neurons = [1:noneurons];

if isfield(sim.perturbation, 'type')
	try
		pert.type
	catch
		pert.type{1} = [];
		pert.triggers{1} = [];
		pert.mask{1} = [];
	end


	nonstim_neurons = [];
	for i = 1:length(pert.type)
		nonstim_neurons = unique( [find(~pert.mask{i})'  nonstim_neurons]);
	end

	plot(V_soma_unwrapped(nonstim_neurons,time_slice)','color', [.7 .7 .7]), hold on
	



		n = 0; colors = jet(length(pert.type));
		for i = 1:length(pert.type)
  
			stn = find(pert.mask{i});
			tr = pert.triggers{i};
			if iscolumn(tr)
				tr = tr';
			end

			if not(isempty(tr))
				switch pert.type{i}
					case 'ampa'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*min(V_soma_unwrapped(:))+10],'color','r')
						plot(V_soma_unwrapped(stn,time_slice)','color', [1 .8 .8])

					case 'gaba_soma'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','b')
						plot(V_soma_unwrapped(stn,time_slice)','color', [.7 .7 1])


					case 'gaba_dend'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','b')
						plot(V_soma_unwrapped(stn,time_slice)','color', [.7 .7 1])

					case 'ou-noise'
						n = n+1;
						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','r')
						plot(V_soma_unwrapped(stn,time_slice)','color', colors(n,:))

					otherwise


				end
			end

		end
	
else
	plot(V_soma_unwrapped(:,time_slice)','color', [.7 .7 .7]), hold on
	plot(mean(V_soma_unwrapped(:,time_slice)),'color', 'k','linewidth', meanlinewidth)
	nonstim_neurons = [1:noneurons];

end


% plot averages



if isfield(sim.perturbation, 'type')
		n = 0;
		for i = 1:length(pert.type)
			stn = find(pert.mask{i});
			if not(isempty(tr))
				switch pert.type{i}
					case 'ampa'
						
						plot(mean(V_soma_unwrapped(stn,time_slice)),'r','linewidth', meanlinewidth)

					case 'gaba_soma'

						
						plot(mean(V_soma_unwrapped(stn,time_slice)),'b','linewidth', meanlinewidth)

					case 'gaba_dend'

						plot(mean(V_soma_unwrapped(stn,time_slice)),'b','linewidth', meanlinewidth)

					case 'ou_noise'
						n = n+1;
						plot(mean(V_soma_unwrapped(stn,time_slice)), 'color', colors(n,:),'linewidth', meanlinewidth)

					otherwise

				end
			end
		end
	
else

	corder = cbrewer('seq', 'Greys',50);
	set(0,'defaultaxescolororder', flip(corder))


	plot(V_soma_unwrapped'), hold on
	meantracecolor = [1 163 218]/255;
	plot(mean(V_soma_unwrapped),'color', meantracecolor,'linewidth', meanlinewidth)
	nonstim_neurons = [1:noneurons];

end




	
	xlabel('ms');ylabel('mV'); 
	xlim([time_slice(1) time_slice(end)])
	axis tight




%    _   __    __    _      __                                       
%   (_)_/_/   / /_  (_)____/ /_____  ____ __________ _____ ___  _____
%    _/_/    / __ \/ / ___/ __/ __ \/ __ `/ ___/ __ `/ __ `__ \/ ___/
%  _/_/_    / / / / (__  ) /_/ /_/ / /_/ / /  / /_/ / / / / / (__  ) 
% /_/ (_)  /_/ /_/_/____/\__/\____/\__, /_/   \__,_/_/ /_/ /_/____/  
%                                 /____/                             
% a(3) = axes('Position',[0.07 0.575 0.275 0.375]);

% bar(edges, histfreq ,'edgecolor', 'none','facecolor', [0 0 0])
% y = get(gca,'ylim');
% x = get(gca,'xlim');
% xlim([-.1 10]);
% set(gca,'xtick',[0:9]);

% xlabel('Hz'); ylabel('Number of cells');
% t4 = text(max(x)*0.95, max(y)*0.40, 50, ['Pop. freq.:' num2str(f) 'Hz']);
% t5 = text(max(x)*0.95, max(y)*0.20,['Sp. neurons.:' num2str(propspkneurons)]);
% set(t4,'backgroundcolor', [1 1 1],'horizontalalignment','right')
% set(t5,'backgroundcolor', [1 1 1],'horizontalalignment','right')





if savemovie
	fname = ['sim' num2str(netsize) '_']
		
		vidObj = VideoWriter(fname,'MPEG-4');
		set(vidObj,'FrameRate', 100)

	open(vidObj);
end


%     ____  __                        
%    / __ \/ /_  ____ _________  _____
%   / /_/ / __ \/ __ `/ ___/ _ \/ ___/
%  / ____/ / / / /_/ (__  )  __(__  ) 
% /_/   /_/ /_/\__,_/____/\___/____/  


% order to present the neurons
a(2) = axes('Position',[0.1 0.55 0.8 0.4]);
O = [1:prod(netsize)];
imagesc(V_soma_unwrapped(O,time_slice),[-66 -45]);

try
	if length(unique(sim.W.stats.clusters))>=2 & length(unique(sim.W.stats.clusters))<30
		[v O] = sort(sim.W.stats.clusters);	

		imagesc(V_soma_unwrapped(O,time_slice),[-66 -35]);
		[u_ yt] = unique(v);
		set(gca,'tickdir','out')

		set(gca,'ytick',yt,'yticklabel', num2str(u_));

		disp('ordering according to cluster identity')


		neurons

	else
		O = [1:prod(netsize)];
		disp('not ordering according to cluster identity')
	end

catch
	
end

hold on;
ylabel('neurons');  

% plotspikes1 = @(c)plot(c,spks.spikes{c},'markersize', 1,'marker', 'o','linestyle', 'none','color', 'w','markersize',6);
plotspikes2 = @(c)plot(spks.spikes{O(c)}-time_slice(1),c, 'markersize', 3,'marker', '.','linestyle', 'none','color', 'k','markersize',18);
plotspikes3 = @(c)plot(spks.spikes{O(c)}-time_slice(1),c, 'markersize', 3,'marker', '*','linestyle', 'none','color', [1 1 1],'markersize',4);
for c = [1:noneurons]
	if not(isempty(spks.spikes{O(c)}))
		% plotspikes1(c);
		plotspikes2(c);
		plotspikes3(c);
	end
end


l = line([1 prod(netsize)], [0 0],'color', 'c');

if static, simtime=1;end

colorbar('east');
title('Vm')

linkaxes(a,'x')




if isfield(sim.networkHistory, 'backgroundnoise')
	if not(isempty(sim.networkHistory.backgroundnoise))
		N = sim.networkHistory.backgroundnoise;

		a(5) = axes;
		O = [1:prod(netsize)];
		imagesc(N(O,time_slice),[-5 5]);
		colorbar('East')
		

		set(a(5),'position',[0.10  0.67  0.8  0.25],'TickLength', [0 0],'xticklabel',[])
		set(a(2),'position',[0.10  0.40  0.8  0.25],'TickLength', [0 0],'xticklabel',[])
		set(a(4),'position',[0.10  0.08  0.8  0.30],'TickLength', [0 0])

	end
end

linkaxes(a,'x')
% colormap(bone(64))


plotclusterneurons = 0;
if plotclusterneurons
	for nc = u_'
		figure
		v==nc
		plot(V_soma_unwrapped(v==nc,time_slice)');
	end

	figure, imagesc(W.W(O,O))
end






if plotallfields
	nsubplots = length(loggedstatevars);

	nsubplotperfig = 5;
	simtime = 		length(sim.networkHistory.V_soma);

	figure, set(gcf,'color', [1 1 1])
	s = 0; na = 0;
	nsubplots = length(loggedstatevars);
	for lsv = loggedstatevars

		s = s+ 1; na = na+1;
			if s > nsubplotperfig
				figure
				set(gcf,'color', [1 1 1])
				s = 1;
			end

			a(na) = subplot(nsubplotperfig, 1, s);

			eval(['plot(transpose( sim.networkHistory.' lsv{1} '));'  ]);
			eval(['legend('  ' lsv '  ');'])
			xlabel('ms')
			axis tight
			% axis off
			if s == nsubplotperfig
				axis on
			end
	end
	linkaxes(a,'x')
end






    
    
    