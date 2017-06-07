function plotSelectedNeurons(sim, time_slice, neuronselection)
	% summaryResults(sim, time_slice, snap_time, neuronselection)

static = 0;

offset = 0; % for mixstim

trigger = 1;

netsize 			= sim.networksize;

noneurons           = prod(netsize);


vsoma = sim.networkHistory.V_soma;
% calcium distribution
% a(31) = axes('Position',[0.2 0.775 0.125 0.170]);
% [g_CaL_hist, g_CaL_hist_x] = hist(g_CaL, [0:0.05:2]);
% bar(g_CaL_hist_x, g_CaL_hist,'edgecolor', 'none','facecolor', [0 0 0]);
% xlim([0 1.75]);
% ylim([0 noneurons*0.275]);




% spks1 = spikedetect(sim.networkHistory.V_soma(group1,:));
% spks2 = spikedetect(sim.networkHistory.V_soma(group2,:));



figure

set(0,'defaultaxescolororder', linspecer(length(neuronselection)))
plot(vsoma(neuronselection,time_slice)','linewidth',2);






% ax2(2) = subplot(3,1,2);
% keyboard
% spkt1 = zeros(1:simtime,1); spkt1.spkspercell = 1;
% spkt2 = zeros(1:simtime,1); spkt2.spkspercell = 1;

% plot(xcorr(V_soma_unwrapped(15,:),V_soma_unwrapped(950,:)));
% set(gca,'colororder',linspecer(5))

% ax2(2) = subplot(3,1,3);
% plot(xcorr(V_soma_unwrapped(15,:),V_soma_unwrapped(950,:),'unbiased'));




