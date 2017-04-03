function [out] = phase_distribution_over_time(varargin)
	% function [out] = phase_distribution_over_time(sim,duration,savemovie)
	% 
	% to plot the distributions we need the pretty style
	% 
	% for 'pretty' style we need the topology toolbox
	% 
	% it's currently commented out
	% phase_distribution_over_time(sim{1}, 'group', find(sim{1}.W.stats.clusters==5))





% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('duration',[])
p.addOptional('savemovie',[])
p.addOptional('animate',0)
p.addOptional('print2file',0)
p.addOptional('trigger',1)

p.addParameter('fhandle', gcf)
p.addParameter('group', [])
p.addParameter('fname', [])

p.parse(varargin{:});

sim = p.Results.sim;
duration = p.Results.duration;
savemovie = p.Results.savemovie;
trigger = p.Results.trigger;
frames2file = p.Results.print2file;
fhandle = p.Results.fhandle;
group = p.Results.group;
fname = p.Results.fname;
animate = p.Results.animate;



% [=================================================================]
%  retrieve fields
% [=================================================================]

netsize 			= sim.networksize;
simtime             = size(sim.networkHistory.V_soma,2);
noneurons           = prod(netsize);
VS 	= sim.networkHistory.V_soma;
pert 				= sim.perturbation;
loggedstatevars = fieldnames(sim.networkHistory)';

% retrieve fields
if ~isempty(sim.perturbation) % backward compatibility
	pert_mask			= pert.mask{trigger};
	perturbation_onsets = pert.triggers{trigger};	
else	
	pert_mask = [];
	perturbation_onsets = [];
end

VS = sim.networkHistory.V_soma;

if frames2file; animate = 1; end

	groupcolors = [.7 .7 .7 ; 0 0 1];

% [=================================================================]
%  slice to plot
% [=================================================================]
 
 if ~isempty(duration)
 	tt = duration;
 else
 	tt = [1:simtime];
 	duration = tt;
 end
 

fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');



%=============================savevideo==============================%

if savemovie
	if isempty(fname);
		% fname = ['volume' datestr(now)];
		fname = 'volume';
	end

		try
			vidObj = VideoWriter('volume','MPEG-4');
		catch
			warning('No MPEG-4 encoder in this crappy OS.')
			vidObj = VideoWriter(['phase' ]);
		end

		vidObj.FrameRate = 24;
	open(vidObj);
end



VS = sim.networkHistory.V_soma;

try
        
 	H = hilbert_of_membranepotential(VS); 
 
catch E

  warning('could not compute hilbert')

  results.stats.firstordersync = 0;
	results.stats.secondordersync = 0;
	results.stats.overallsync =     0;
	results.stats.order_parameter_g_fo = 0;
	results.stats.order_parameter_g_so = 0;
	results.hilbert = 0;

  return 

end



U = H.hilbert;

allneu = 1:noneurons;


GA = U(allneu,:);
MA = circ_mean(GA)+pi;
VA = circ_var(GA);

GR = U(group,:);
MR = circ_mean(GR)+pi;
VR = circ_var(GR);


% - select groups
% group2=[];
if ~isempty(group)
		allneu = [1:noneurons];
		group1 = group;
		group2 = setdiff(allneu, group1);
else
		allneu = (pert_mask>=0);
		group1 = find(pert_mask==1); % stimulated group
		group2 = find(pert_mask==2);
end

if isempty(group2)

		% group1 = [1:floor(noneurons/2)];
		allneu = [1:noneurons];
		group2 = setdiff(allneu, group1);
end




G1 = U(group1,:);
G2 = U(group2,:);
GA = U(allneu,:);


% circular stats
M1 = circ_mean(G1)+pi;
M2 = circ_mean(G2)+pi;
MA = circ_mean(GA)+pi;
V1 = circ_var(G1);
V2 = circ_var(G2);
VA = circ_var(GA);


% mean of exp( e^(i*(theta_k,p(t) - theta_syn(t)) ; p is group, k is neuron
	% assuming G2 is the 'non stimulated group'
	% KURAMOTO ORDER PARAMETER
	order_parameter_G1 = mean( exp(i*(G1+pi)));
	order_parameter_G2 = mean( exp(i*(G2+pi)));
	order_parameter_GA = mean( exp(i*(GA+pi)));



%            __      __      
%     ____  / /___  / /______
%    / __ \/ / __ \/ __/ ___/
%   / /_/ / / /_/ / /_(__  ) 
%  / .___/_/\____/\__/____/  
% /_/                        


fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


fig0 = figure('color',[1 1 1])
	subplot(2,2,1:2)
	plot(abs(order_parameter_G1),'r'),hold on
	plot(abs(order_parameter_G2),'b')
	plot(abs(order_parameter_GA),'g')
	title('abs(Z)')
	legend({'Group' 'Others' 'All'})
	
	subplot(2,2,3)
	scatter(real(order_parameter_G1(tt)),imag(order_parameter_G1(tt)), 'filled', 'r' ); hold on
	scatter(real(order_parameter_G2(tt)),imag(order_parameter_G2(tt)), 'filled', 'b' )
	scatter(real(order_parameter_GA(tt)),imag(order_parameter_GA(tt)), 'filled', 'g' )
	hold on
	axis equal
	axis([-1,1,-1,1])
	title('')
	legend({'Group' 'Others' 'All'})
	


fig0 = figure('color',[1 1 1])
	% circular phase plot
	a(1) = axes;

		circ_plot(G1(:,1),'pretty','',true,'linewidth',1,'color',groupcolors(1,:),'markersize',7,'marker','o');
		circ_plot(G2(:,1),'pretty','',true,'linewidth',1,'color',groupcolors(2,:),'markersize',10,'marker','.');
		hold on

		axis off

		alpha(.3)
		drawnow
		
		title({num2str(1);[' ']; [' ']})

	%static phase plot over time
fig1 = figure('color',[1 1 1])
	ttt  = tt(1):tt(end);
	 ttt1 = repmat(tt(1):tt(end), numel(group1), 1);  
	 ttt2 = repmat(tt(1):tt(end), numel(group2), 1);
	a(2) = subplot(2,1,1);

	fill_between_lines(ttt, M1(ttt)+V1(ttt), M1(ttt)-V1(ttt),[1 .8 .8])
		hold on
		fill_between_lines(ttt, M2(ttt)+V2(ttt), M2(ttt)-V2(ttt),[0 .8 1])
		hold on
		line(ttt1,M1(ttt1),'color', groupcolors(1,:)); 
		hold on
		line(ttt2,M2(ttt2),'color', groupcolors(2,:));

		title('phase average')
		xlabel('ms')
		ylabel('phase (radians)')
		axis tight



	a(3) = subplot(2,1,2)

	 % plot(sim.networkHistory.V_soma')
	 
	 

		line(ttt1',sim.networkHistory.V_soma(group1,ttt)', ...
			'color',groupcolors(1,:), 'linewidth', .5)

		line(ttt2',sim.networkHistory.V_soma(group2,ttt)', ...
			'color',groupcolors(2,:), 'linewidth', .5)
		
	   line(ttt, mean(sim.networkHistory.V_soma(group1,ttt)),'color',groupcolors(1,:),'linewidth', 3)
	   line(ttt, mean(sim.networkHistory.V_soma(group2,ttt)),'color',groupcolors(2,:),'linewidth', 3)
	   
	 
		 xlabel('time')
		 ylabel('Vm @ soma')
		 axis tight





fig2 = figure('color',[1 1 1])
	[phasedist1 tim] = hist(G1, length(tt));
	[phasedist2 tim] = hist(G2, length(tt));
	[phasedistA tim] = hist(GA, length(tt));

	ax_ph(1) = subplot(3,1,1)
	imagesc(1:length(G1),tim, phasedist1), axis xy
	colorbar
	title('phase distribution over time (group1)')
	xlabel('time (ms)')
	ylabel('Prob(phase)')

	ax_ph(2) = subplot(3,1,2)
	imagesc(1:length(G2),tim, phasedist2), axis xy
	colorbar
	title('phase distribution over time (group2)')
	
	ax_ph(3) = subplot(3,1,3)
	imagesc(1:length(GA),tim, phasedistA), axis xy
	colorbar
	title('phase distribution over time (all neurons)')

	linkaxes(ax_ph, 'x')
	xlim()



spks1 = spikedetect(sim.networkHistory.V_soma(group1,:));
spks2 = spikedetect(sim.networkHistory.V_soma(group2,:));


if animate
	for t = duration 
		% R1 = rose(U(group1,t));
		% R2 = rose(U(group2,t));

		axes(a(1))
		cla
		circ_plot(G2(:,t),'pretty','',true,'linewidth',1,'color',groupcolors(1,:),'markersize',20,'marker','.');
		hold on
		circ_plot(G1(:,t),'pretty','',true,'linewidth',1,'color',groupcolors(2,:),'markersize',7,'marker','o');
		
		% hold on
		% circ_plot(G2(:,t),'pretty','',true,'linewidth',1,'color','b','markersize',10,'marker','.');

		
		axis off

		alpha(.3)
		drawnow
		
		title({num2str(t);[' ']; [' ']})


	     if savemovie	
	  
	   		if false
	   			currFrame = getframe(a(1));
	   		else
	       		currFrame = getframe(fig0);
	       	end
	       writeVideo(vidObj,currFrame);
		else
			drawnow;
			% keyboard
		end

		if frames2file
			name = ['phase_dist@' num2str(t) ];

			frm=hgexport('readstyle','12x12');
		    frm.Format = 'pdf';
		    hgexport(fig0,name,frm);


			
		end
	end
end



if savemovie
	% Close the file.
    close(vidObj);
end


out.phases.mean = [M1 ; M2];
out.phases.var  = [circ_var(G1) ; circ_var(G2)];
out.phases.groups = {group1 ; group2};
out.phases.pop = {G1 ; G2};
out.phases.orderparameter = {order_parameter_G1 ; order_parameter_G2};
out.spikes.spks1 = spks1;
out.spikes.spks2 = spks2;


% for i = 1:1000;

% X(i,:) = hist(scaled(:,i),72)
% fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


% end

% % out.kurtosis = 
% % out.skew = 



function out = spikedetect(vsoma) 

% Detect spikes with threshold
noneurons = size(vsoma,1);

transient_to_exclude = 1;

A = zeros(noneurons,21);

for i = 1:noneurons
   abovethreshold{i} = find(vsoma(i,transient_to_exclude:end)>-30);
   b{i} = diff([-2 abovethreshold{i}]);
   spikeonset = find(b{i}>1);
   d{i} = abovethreshold{i}(spikeonset);
   ISI{i} = diff(d{i});
   freq{i} = 1./ISI{i}*1000;
   med{i} = median(freq{i});
     
        
    if ~isempty(freq{i})
        B = histc(freq{i},[0:.5:10]);  
        A(i,:) = B;
        
    end
    
    if isempty(freq{i})
        freq{i} = [0];
        med{i} = [0];
    end
    
end

popmedian = cell2mat(med);

popfreq = median(popmedian);

% Find proportion of spiking cells
spkneurons = [find(~cellfun(@isempty,d))];
propspkneurons = numel(spkneurons)/(noneurons);

 
for i = 1:noneurons;
	spkspercell(i) = numel(d{i});
end




out.spikespercell = spkspercell;
out.medfreq = popmedian;
out.popfrequency = popfreq;
out.propspkneurons = propspkneurons;
out.cellISI = ISI;



