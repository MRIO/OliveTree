function out = network_phase_response(sims)


mph = -55;
mpd = 50;
wave = [50:170];
rows = sims{1}.rows;
cols = sims{1}.columns;
net_size = rows*cols;

bins = 15;

for s = 1:length(sims)
	[pos(s,:) amp(s,:)] = phase_response(sims{s},wave);
	length(find(amp>-10))
end

% peak distribution


figure
set(gcf, 'defaultaxescolororder', linspecer(length(sims)))
[h xx] = hist(pos'+wave(1),bins,'edgecolor','none');
plot(xx',h/net_size','linewidth',3)
set(gca,'fontsize',12)
xlabel('latency to peak')
ylabel('probability of latency')
% legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})

figure

i=1; center = [ceil(rows/2) ceil(cols/2)];
for x = 1:rows
	for y =1:cols

		D(i) = (x-center(1))^2+(y-center(2))^2;
		i = i+1;
	end
end
[V O] = sort(D);


set(gcf, 'defaultaxescolororder', linspecer(net_size))

nsp = length(sims);
ha = tight_subplot(nsp,1,[.01 .03],[.05 .01],[.05 .01]);
for s =1:nsp
	vsoma = sims{s}.networkHistory.V_soma;
	
        axes(ha(s));
        plot(vsoma(O,1:size(vsoma,2))')

		axis tight
		ylim([-65 -20])

end
	set(ha(1:length(sims)-1),'XTickLabel','');
	xlabel('ms'), ylabel('mV')


out.histogram = h-wave(1);
out.xx = xx;




function [pos amp] = phase_response(sim,wave)

networkHistory  	= sim.networkHistory;
rows 				= sim.rows;
columns				= sim.columns;
lastState 			= sim.lastState;
dt 					= sim.dt;
g_CaL 				= sim.g_CaL;
g_Gap 				= sim.g_Gap;
V_soma_unwrapped 	= sim.networkHistory.V_soma;

if isfield(sim, 'simtime')
	simtime				= sim.condition.simtime;

pert_map			= sim.condition.perturbation_map;
perturb_amplitude	= sim.condition.perturb_amplitude;
perturb_onsets      = sim.condition.perturb_onsets;
noise_level			= sim.condition.noise_level;
else
	simtime = sim.time;
end


net_size = rows*columns;

time_slice = [1:simtime];


% simtime = [1:length(time_slice)];

V_soma_unwrapped = V_soma_unwrapped(:,time_slice);

V_soma = reshape(V_soma_unwrapped,rows, columns, []);

% spks = getSpks(V_soma_unwrapped);

% if ishandle(gcf)
% 	clf
% end


i=1; center = [ceil(rows/2) ceil(columns/2)];
for x = 1:rows
	for y =1:columns

		D(i) = (x-center(1))^2+(y-center(2))^2;
		i = i+1;
	end
end
[V O] = sort(D);


% first peak distribution
mph = -55;
mpd = 50;

vsoma = sim.networkHistory.V_soma;

% keyboard
for c = 1:net_size
 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph,'minpeakdistance',mpd );
 [amp(c) pos(c) ] = max(vsoma(c,wave));
 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph);
end



[h x] = hist((pos)',10,'edgecolor','none');

% figure
% set(gcf, 'defaultaxescolororder', linspecer(5))
% plot(x',h/net_size','linewidth',3)
% set(gca,'fontsize',12)
% xlabel('latency to peak')
% ylabel('probability of latency')
% legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
% plot2svg('first_pk_hype_bars.svg')
% plot2svg(['first_pk_' simtype '_lines.svg'])
% close


% figure
% bar(x,h,'edgecolor','none')
% legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
% xlabel('latency to peak')
% ylabel('probability of latency')

% plot2svg('first_pk_hype_bars.svg')
% close
