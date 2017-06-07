function [spks] = replayResults(varargin)
% p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
% p.addRequired('time_slice')
% p.addRequired('savemovie')
% p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
% p.addParamValue('fhandle', gcf)
% TODO: accept any dt

static = 0;

% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('time_slice',[])
p.addOptional('savemovie',0)

p.addParamValue('plotallfields', 0) % stdandard deviation criterion for offset threshold
p.addParamValue('fhandle', gcf)

p.parse(varargin{:});

sim = p.Results.sim;
time_slice = p.Results.time_slice;
savemovie = p.Results.savemovie;

plotallfields  = p.Results.plotallfields;
fhandle = p.Results.fhandle;


% [=================================================================]
%  retrieve fields
% [=================================================================]

netsize 			= sim.networksize;
simtime             = length(sim.networkHistory.V_soma);
noneurons           = prod(netsize);
V_soma_unwrapped 	= sim.networkHistory.V_soma;
V_soma = reshape(V_soma_unwrapped(:,1) , noneurons, 1);
pert 				= sim.perturbation;
loggedstatevars = fieldnames(sim.networkHistory)';


% [=================================================================]
%  slice to plot
% [=================================================================]
 
 if ~isempty(time_slice)
 	simtime = length(time_slice);
 else
 	time_slice = [1:2];
 end
 

% [=================================================================]
%  spike detection
% [=================================================================]

spks = spikedetect(sim);

propspkneurons = spks.propspkneurons;
popfreq = spks.popfrequency;

	f = sprintf('%.2f', popfreq);
	edges = [-1 0.01:.5:10];
	histfreq = histc(spks.medfreq, edges);

V_soma = reshape(V_soma_unwrapped(:,1) ,prod(netsize), 1);

% order to present the neurons
O = [1:prod(netsize)];


% [=================================================================]
%  if saving movie
% [=================================================================]

if savemovie
	fname = ['sim' num2str(netsize) '_']
		
		vidObj = VideoWriter(fname,'MPEG-4');
		set(vidObj,'FrameRate', 100)

	open(vidObj);
end

V_soma_unwrapped = V_soma_unwrapped(O,:);


% [=================================================================]
%  prepare figure
% [=================================================================]


fig1 = gcf;
set(fig1, 'position', [1 1 1024 768],'color', [1 1 1]);	
colormap(hot);
clf




% [=================================================================]
%  plot populations according to type of perturbation
% [=================================================================]

a(4) = axes('Position',[0.07 0.07 0.575 0.40]);
nonstim_neurons = [1:noneurons];

if isfield(sim, 'perturbation')
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

	plot(V_soma_unwrapped(nonstim_neurons,:)','color', [.7 1 .7]), hold on
	



		n = 0; colors = jet(length(pert.type));
		for i = 1:length(pert.type)
  
			stn = find(pert.mask{i});

			tr = pert.triggers{i}';
			if iscolumn(tr)
				tr = tr';
			end

			if not(isempty(tr))
				hold on
				switch pert.type{i}
					case 'ampa'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','r')
						plot(V_soma_unwrapped(stn,:)','color', [1 .7 .7])

					case 'gaba_soma'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','b')
						plot(V_soma_unwrapped(stn,:)','color', [.7 .7 1])


					case 'gaba_dend'

						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color','r')
						plot(V_soma_unwrapped(stn,:)','color', [.7 .7 1])

					case 'ou_noise'
						n = n+1;
						stn
						line([tr; tr] , [ones(size(tr))*min(V_soma_unwrapped(:)); ones(size(tr))*max(V_soma_unwrapped(:))],'color',colors(n,:))
						plot(V_soma_unwrapped(stn,:)','color', colors(n,:))


					otherwise


				end
			end

		end
	
else
	plot(V_soma_unwrapped','color', [.7 1 .7]), hold on
	plot(mean(V_soma_unwrapped),'g','linewidth', 3)
	nonstim_neurons = [1:noneurons];

end



% plot averages
if isfield(sim, 'perturbation')
		n = 0;
		for i = 1:length(pert.type)
			stn = find(pert.mask{i});
			if not(isempty(tr))
				switch pert.type{i}
					case 'ampa'
						
						plot(mean(V_soma_unwrapped(stn,:)),'r','linewidth', 3)

					case 'gaba_soma'

						
						plot(mean(V_soma_unwrapped(stn,:)),'b','linewidth', 3)

					case 'gaba_dend'

						plot(mean(V_soma_unwrapped(stn,:)),'b','linewidth', 3)

					case 'ou-noise'
						n = n+1;
						plot(mean(V_soma_unwrapped(stn,:)), 'color', colors(n,:),'linewidth', 3)

					case 'ou-noise_pos'
						n = n+1;
						plot(mean(V_soma_unwrapped(stn,:)), 'color', colors(n,:),'linewidth', 3)

					case 'ou-noise_neg'

						n = n+1;
						plot(mean(V_soma_unwrapped(stn,:)), 'color', colors(n,:),'linewidth', 3)

					otherwise

				end
			end
		end
	
else
	plot(V_soma_unwrapped','color', [.7 1 .7]), hold on
	plot(mean(V_soma_unwrapped),'g','linewidth', 3)
	nonstim_neurons = [1:noneurons];

end


plot(mean(V_soma_unwrapped(nonstim_neurons,:)),'g','linewidth', 3)

	
	xlabel('ms');ylabel('mV'); 
	xlim([time_slice(1) time_slice(end)])
	axis tight




%    _   __    __    _      __                                       
%   (_)_/_/   / /_  (_)____/ /_____  ____ __________ _____ ___  _____
%    _/_/    / __ \/ / ___/ __/ __ \/ __ `/ ___/ __ `/ __ `__ \/ ___/
%  _/_/_    / / / / (__  ) /_/ /_/ / /_/ / /  / /_/ / / / / / (__  ) 
% /_/ (_)  /_/ /_/_/____/\__/\____/\__, /_/   \__,_/_/ /_/ /_/____/  
%                                 /____/                             
a(3) = axes('Position',[0.07 0.575 0.275 0.375]);

% imagesc(cov(sim.networkHistory.V_soma))
bar(edges, histfreq ,'edgecolor', 'none','facecolor', [0 0 0])
y = get(gca,'ylim');
x = get(gca,'xlim');
xlim([-.1 10]);
set(gca,'xtick',[0:9]);

xlabel('Hz'); ylabel('Number of cells');
t4 = text(max(x)*0.95, max(y)*0.40, 50, ['Pop. freq.:' num2str(f) 'Hz']);
t5 = text(max(x)*0.95, max(y)*0.20,['Sp. neurons.:' num2str(propspkneurons)]);
set(t4,'backgroundcolor', [1 1 1],'horizontalalignment','right')
set(t5,'backgroundcolor', [1 1 1],'horizontalalignment','right')



% order to present the neurons
O = [1:prod(netsize)];

if savemovie
	fname = ['sim' num2str(netsize) '_']
		
		vidObj = VideoWriter(fname,'MPEG-4');
		set(vidObj,'FrameRate', 100)

	open(vidObj);
end

V_soma_unwrapped = V_soma_unwrapped(O,:);




%     ____  __                        
%    / __ \/ /_  ____ _________  _____
%   / /_/ / __ \/ __ `/ ___/ _ \/ ___/
%  / ____/ / / / /_/ (__  )  __(__  ) 
% /_/   /_/ /_/\__,_/____/\___/____/  
                                    

a(2) = axes('Position',[0.70 0.07 0.275 0.85]);
imagesc(V_soma_unwrapped',[-66 -30]);

set(gca,'xtick',[1 noneurons]);
hold on;
  ylabel('ms');
  xlabel('neurons');  



% plotspikes1 = @(c)plot(c,spks.spikes{c},'markersize', 1,'marker', 'o','linestyle', 'none','color', 'w','markersize',6);
plotspikes2 = @(c)plot(c,spks.spikes{c},'markersize', 3,'marker', '.','linestyle', 'none','color', 'g','markersize',10);
for c = 1:prod(netsize)
	if not(isempty(spks.spikes{c}))
		% plotspikes1(c);
		plotspikes2(c);
	end
end


l = line([1 prod(netsize)], [0 0],'color', 'c');

if static, simtime=1;end

colorbar;
title('Vm')

a(1) = axes('Position',[0.375 0.575 0.275 0.375]);

set(0,'CurrentFigure',fig1);
    set(fig1,'CurrentAxes',a(1));
	imagesc(V_soma_unwrapped(:,:,1),[-68 -30]);


gg = gausswin(5)*gausswin(5)'; gg = gg / sum(gg(:));

for t =time_slice;
    set(0,'CurrentFigure',fig1);
    set(fig1,'CurrentAxes',a(1));



	% state = V_soma_unwrapped(:,t);
	state = reshape(V_soma_unwrapped(:,t),netsize(1)*netsize(2), netsize(3));
	% state = reshape(V_soma_unwrapped(:,t),netsize(1), netsize(2)* netsize(3));

	imagesc(state',[-68 30])
	axis equal, axis xy, axis off

	title(num2str(t));

	[i j] = find(state>-30);
	if ~isempty(i)
		line(i,j,'markersize', 10,'marker', '.','linestyle', 'none','color', 'c')
	end
   
    set(l,'ydata', [t t])
  


   %    _   __    ____                 __                         _     __         
   %   (_)_/_/   / __ \_________  ____/ /_  __________     _   __(_)___/ /__  ____ 
   %    _/_/    / /_/ / ___/ __ \/ __  / / / / ___/ _ \   | | / / / __  / _ \/ __ \
   %  _/_/_    / ____/ /  / /_/ / /_/ / /_/ / /__/  __/   | |/ / / /_/ /  __/ /_/ /
   % /_/ (_)  /_/   /_/   \____/\__,_/\__,_/\___/\___/    |___/_/\__,_/\___/\____/ 
                                                                                 
    if savemovie
  
   		if false
   			currFrame = getframe(a(1));
   		else
       		currFrame = getframe(fig1);
       	end
       writeVideo(vidObj,currFrame);
	else
		drawnow;
		% keyboard
	end

end        







if savemovie
	% Close the file.
    close(vidObj);
    
end



%write video
   % Prepare the new file.
   



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











