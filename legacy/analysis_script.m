
mph = -30;
mpd = 125;
vsoma = simIO_gpu_transients.networkHistory.V_soma;
wave = [1:simIO_gpu_transients.duration];

s = 1;
netsize = simIO_gpu{s}.rows*simIO_gpu{s}.columns;
for c = 1:netsize
	putativePK = findpeaks(vsoma(c,wave), 'minpeakheight',mph,'minpeakdistance',mpd );
	if ~isempty(putativePK)
 		[amp(c) pos(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph,'minpeakdistance',mpd );
 	else
 		amp(c)=0;
 		pos(c)=0;
 	end
 	c
end



§

subplot(2,1,1)
hist(amp,50)
skew_amp = skewness(amp);
ylabel('neurons')
xlabel('STO Peak Amplitude')
subplot(2,1,2)
hist(g_CaL,50)
skew_CaL = skewness(g_CaL);
xlabel('Low Threshold Calcium Conductance')



i=1;
for x = 1:50
	for y =1:50

		D(i) = (x-26)^2+(y-25)^2;
		i = i+1;
	end
end
[V O] = sort(D);



	%     ____  ___________________  ________ 
	%    / __ \/  _/ ___/_  __/ __ \/  _/ __ )
	%   / / / // / \__ \ / / / /_/ // // __  |
	%  / /_/ // / ___/ // / / _, _// // /_/ / 
	% /_____/___//____//_/ /_/ |_/___/_____/  
	                                        

% noisy_transients.mat sim_depo_perturb.mat sim_hyper_pert.mat   sims.mat             sims_depol_pert.mat  sims_mixed.mat       sims_transients.mat


% load sim_depo_perturb
% simtype = 'depo';

load sim_hyper_pert.mat %sim_hyper_pert
simtype = 'hype';

% first peak distribution
mph = -55;
mpd = 50;
wave = [200:280]	
for i =1:length(simIO_gpu)
	vsoma = simIO_gpu{i}.networkHistory.V_soma;

	for c = 1:2500
	 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph,'minpeakdistance',mpd );
	 [amp(i,c) pos(i,c) ] = max(vsoma(c,wave));
	 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph);
	end

end

set(0, 'defaultaxescolororder', linspecer(250))

[h x] = hist((pos)',10,'edgecolor','none');

figure
plot(x',h/2500','linewidth',3)
set(gca,'fontsize',12)
xlabel('latency to peak')
ylabel('probability of latency')
legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
% plot2svg('first_pk_hype_bars.svg')
plot2svg(['first_pk_' simtype '_lines.svg'])
close


% figure
% bar(x,h,'edgecolor','none')
% legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
% xlabel('latency to peak')
% ylabel('probability of latency')

% plot2svg('first_pk_hype_bars.svg')
% close


% second peak distribution
mph = -55;
mpd = 50;
wave = [201:350]
for i =1:length(simIO_gpu)
	vsoma = simIO_gpu{i}.networkHistory.V_soma;

	for c = 1:2500
	 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph,'minpeakdistance',mpd );
	 [amp(i,c) pos(i,c) ] = max(vsoma(c,wave));
	 % [amp{i}(c) pos{i}(c) ] = findpeaks(vsoma(c,wave), 'minpeakheight',mph);
	end

end
% figure
[h x] = hist((pos-52)',10,'edgecolor','none');
plot(x',h/2500','linewidth',3)
set(gca,'fontsize',12)
xlabel('latency to peak')
ylabel('probability of latency')
legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
plot2svg(['second_pk_' simtype '_lines.svg'])
% close
% plot2svg('second_pk_hype_bars.svg')


% figure
% bar(x,h,'edgecolor','none')
% xlabel('latency to peak')
% ylabel('probability of latency')
% legend({'post peak', 'early trough', 'trough', 'rising' ,'peak'})
% plot2svg('second_pk_hype_bars.svg')
% close





	%  _       _____ _    _____________
	% | |     / /   | |  / / ____/ ___/
	% | | /| / / /| | | / / __/  \__ \ 
	% | |/ |/ / ___ | |/ / /___ ___/ / 
	% |__/|__/_/  |_|___/_____//____/  
	                                 



set(0, 'defaultaxescolororder', linspecer(2500))
% depoloarization waves
clf
nsp = length(simIO_gpu);
ha = tight_subplot(nsp,1,[.01 .03],[.05 .01],[.05 .01])
for s =1:nsp
	vsoma = simIO_gpu{s}.networkHistory.V_soma;
	
        axes(ha(s));
        plot(vsoma(O,1:250)')
	
		axis tight
		ylim([-65 -20])

end
	set(ha(1:length(simIO_gpu)),'XTickLabel','');
	xlabel('ms'), ylabel('mV')




% hypoerpolarizations
clf
ha = tight_subplot(5,1,[.01 .03],[.05 .01],[.05 .01])
for s =1:length(simIO_gpu)
	vsoma = simIO_gpu{s}.networkHistory.V_soma;
	
        axes(ha(s));
        plot(vsoma(O,201:400)')
	
		axis tight
		ylim([-75 -30])

end
	set(ha(1:4),'XTickLabel','');
	xlabel('ms'), ylabel('mV')


% MOVIES

for sim_ = 1:4
	figure(1),clf	
	replayResults(simIO_gpu{sim_},[],1)
	eval(['!mv sim.avi sim_rhytmic_hy_pert' num2str(sim_)])
end









clf
ha = tight_subplot(1,4,[.01 .03],[.05 .01],[.05 .01])
for s =1:4
	vsoma = simIO_mix_depo{s}.networkHistory.V_soma;
	
        axes(ha(s));
        imagesc(vsoma(O,1:300)',[-75 -30])
		

end
	set(ha(1:4),'XTickLabel','');
	xlabel('neurons'), ylabel('ms')
	colorbar





clf
ha = tight_subplot(1,4,[.01 .03],[.05 .01],[.05 .01])
for s =1:4
	vsoma = simIO_mix_depo{s}.networkHistory.V_soma;
	
        axes(ha(s));
        imagesc(vsoma(O,1:300)',[-75 -30])
		

end
	set(ha(1:4),'XTickLabel','');
	xlabel('neurons'), ylabel('ms')
	colorbar



