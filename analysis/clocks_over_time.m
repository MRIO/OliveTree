function [out] = clocks_over_time(varargin)
	% only works if there are groups defined according to input perturbation!
	% function [out] = clocks_over_time(sim,duration, frames)
	% reuires: plot2svg

static = 0;

% [=================================================================]
%  Parser
% [=================================================================]

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addOptional('duration',[])
p.addOptional('frames',[])
p.addOptional('print2file',0)
p.addOptional('trigger',1)

p.addParamValue('fhandle', gcf)

p.parse(varargin{:});

sim = p.Results.sim;
duration = p.Results.duration;
frames = p.Results.frames;
trigger = p.Results.trigger;
print2file = p.Results.print2file;
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
 
 if ~isempty(duration)
 	tt = duration;
 else
 	tt = [1:simtime];
 end
 

fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


VS = sim.networkHistory.V_soma;
noneurons = prod(sim.networksize);
simtime = size(VS,2);


userdefinedgroup = [];


% retrieve fields
if ~isempty(sim.perturbation) % backward compatibility
		pert_mask			= sim.perturbation.mask{trigger};
		perturbation_onsets = sim.perturbation.triggers{trigger};	
else	
	pert_mask = [];
	perturbation_onsets = [];
end



if isempty(duration)
	duration = [1:simtime];
end
	tt = [duration(1):duration(end)];

if isempty(frames)
	frames = tt;
end


%%%%%%
H = hilbert_of_membranepotential(VS);
U = H.hilbert;
%%%%%%

% VS(VS>-30) = -30;
% N = -VS;

% deno = (max(N,[],2)-min(N,[],2))';
% scaled = bsxfun(@minus, N, min(N,[],2));
% scaled = bsxfun(@rdivide, scaled, deno');
% scaled = scaled *2 -1;


% for n = 1:noneurons
% 	U(n,:) = angle(hilbert(scaled(n,:)));
% end


if ~isempty(userdefinedgroup)
		group1 = userdefinedgroup;
		allneu = [1:noneurons];
		group2 = setdiff(group1, allneu)
	else
	allneu = (pert_mask>=0);
	group1 = (pert_mask==1); % stimulated group
	group2 = (pert_mask==0);

end



% CK1 = circ_kurtosis(U(group1,:),[],2);
% CK2 = circ_kurtosis(U(group2,:),[],2);

G1 = U(group1,:);
G2 = U(group2,:);
GA = U(allneu,:);

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

	% order_parameter_G1 = mean( exp(i*bsxfun(@minus, G1+pi, M1)));
	% order_parameter_G2 = mean( exp(i*bsxfun(@minus, G2+pi, M2)));

% find intervals with significant differences


f100 = figure(100)
	subplot(2,1,1)
	plot(abs(order_parameter_G1),'r'),hold on
	plot(abs(order_parameter_G2),'b')
	plot(abs(order_parameter_G2) - abs(order_parameter_G1),'g')
	plot(abs(order_parameter_GA),'k')

	legend({'group 1' 'group 2' 'difference' 'all neurons'})
	
	subplot(2,1,2)
	scatter(real(order_parameter_G2(tt)),imag(order_parameter_G2(tt)),'filled','b' )
	hold on
	scatter(real(order_parameter_G1(tt)),imag(order_parameter_G1(tt)) ,'filled','r')
	


f101= figure(101)

a(2) = subplot(1,2,1);

fill_between_lines(tt, M1(tt)+V1(tt), M1(tt)-V1(tt),[1 .8 .8])
hold on
fill_between_lines(tt, M2(tt)+V2(tt), M2(tt)-V2(tt),[0 .8 1])
hold on
plot(tt,M1(tt),'r'); 
hold on
plot(tt,M2(tt),'b');
title('phase average')
xlabel('ms')
ylabel('phase (radians)')
axis tight

a(3) = subplot(1,2,2)

 % plot(sim.networkHistory.V_soma')

 
 ttt = repmat(tt, numel(group2), 1);
	
	line(tt',sim.networkHistory.V_soma(group2,tt)', ...
		'color',[.8 .8 1], 'linewidth', .5)
	line(tt',sim.networkHistory.V_soma(group1,tt)', ...
		'color',[1 .8 .8], 'linewidth', .5)

  hold on
   line(tt, mean(sim.networkHistory.V_soma(group2,tt)),'color','b','linewidth', 3)
   line(tt, mean(sim.networkHistory.V_soma(group1,tt)),'color','r','linewidth', 3)
	 xlabel('time')
	 ylabel('Vm @ soma')
	 % export_fig('Vm', '-r300')
	 % export_fig('Vm2', '-r300')
	 % ylabel('simulation run')

spks1 = spikedetect(sim.networkHistory.V_soma(group1,:));
spks2 = spikedetect(sim.networkHistory.V_soma(group2,:));


	numsubplots = length(frames);
	fclock = figure
	aclock = axes;
	
	sbp = 0;
	for t = frames
		% sbp = sbp +1;
		% axx = subplot(1,numsubplots,sbp);
		
		% axes(aclock)
		% cla
		% circ_plot(G2(:,t)+pi,'oldpretty','',true,'linewidth',1,'color','b','markersize',15,'marker','+');
		% hold on
		% circ_plot(G1(:,t)+pi,'oldpretty','',true,'linewidth',1,'color','r','markersize',10,'marker','o');
		% axis off

		circ_plot(G2(:,t)+pi,'pretty','',true,'linewidth',1,'color','b','markersize',15,'marker','+');
		hold on
		circ_plot(G1(:,t)+pi,'pretty','',true,'linewidth',1,'color','r','markersize',10,'marker','o');
		axis off

		drawnow
		alpha(.3)
		
		title(num2str(t))

		if print2file
			plot2svg(['phase_dist@' num2str(t) '.svg'])
		end

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



