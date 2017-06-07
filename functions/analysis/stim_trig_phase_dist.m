% stim_trig_phase_dist.m

function out = stim_trig_phase_dist(sim)

trigger = 1;
spontaneous = 0;

% window around trigger
win = 1000;
netsize = sim.networksize;
noneurons = prod(netsize);

duration = [1:100000];
triggers = sim.perturbation.triggers{trigger};
triggers(triggers>duration(end)-win-1) = [];
triggers = triggers(1:end-1)';
notriggers = length(triggers);


fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


W = sim.networkParameters.connectivityMatrix;
gapneighborhood = full(sum(triu(W)))+eps;


spks = sim.spikes;
Q = quantile(spks.spikespercell,[.25 , .5, .75]);


	mask1 = double(sim.perturbation.mask{trigger});
	W1 = double(sim.networkParameters.connectivityMatrix>0);
	Vm1 =  sim.networkHistory.V_soma;

	
	fo_neighbors1 = setdiff(find(W1 * mask1), find(mask1));
	neither_nor1  = setdiff([1:prod(netsize)], union(fo_neighbors1, find(mask1)));

	

	GS1 = measureGlobalSync(sim, 'duration', duration, 'plotme', 0);
	H1 = GS1.hilbert.hilbert;
	 
	[hH1 x] = hist(H1,100);

	
	Z1 = zeros(100,2*win+1);
	
	 for t = triggers
	 	Z1 = Z1 + hH1(:, t-win:t+win);

	 end
	 Z1 = Z1/notriggers;

	% sim1
	n =0; 
	for t = triggers
		n = n+1;
		interv = t-win:t+win;
		
		Km = H1(find(mask1), interv );
		KKm(n,:)  = abs(mean(exp(i*bsxfun(@minus, Km, circ_mean(Km)))));

		Kfo = H1(fo_neighbors1, interv );
		KKfo(n,:) = abs(mean(exp(i*bsxfun(@minus, Kfo, circ_mean(Km)))));

		Ko =  H1(neither_nor1, interv);
		KKo(n,:)  = abs(mean(exp(i*bsxfun(@minus, Ko, circ_mean(Km)))));

		% figure(10000)
		% subplot(3,1,1)
		% imagesc(Km)
		% subplot(3,1,2)
		% imagesc(Kfo)
		% subplot(3,1,3)
		% imagesc(Ko)
		% pause(1)

	end


	R1.kp_mask = KKm;
	R1.kp_neig = KKfo;
	R1.kp_othr  = KKo;


	n =0; 
	for t = triggers
		n = n+1;
		interv = t-win:t+win;
		Vneig = Vm1(fo_neighbors1, interv);
		VVneig(n,:) = mean(Vneig);

		Vothr = Vm1(neither_nor1, interv);
		VVothr(n,:) = mean(Vothr);

		Vmask  = Vm1(find(mask1), interv);
		VVmask(n,:) = mean(Vmask);

		Vnm = Vm1(find(not(mask1)), interv);
		VVnm(n,:) = mean(Vnm);	

		% imagesc(Vneig'),hold on
		% plot(mean(Vneig),'r', 'linewidth',3)
		% pause(1)
		% clf
	end

	R1.meanVm_mask = mean(VVmask);
	R1.meanVm_neig = mean(VVneig);
	R1.meanVm_othr = mean(VVothr);
	R1.meanVm_notm = mean(VVnm);

	R1.Vm_mask = VVmask;
	R1.Vm_neig = VVneig;
	R1.Vm_othr = VVothr;
	R1.Vm_notm = VVnm;



xax = linspace(-win,win,win*2+1);




figure
	imagesc(Z1,[0,20]) 
	title('stimulus triggered phase concentration')
	xlabel('ms')

	% triggered kuramoto parameter distribution 
	% #todo, make it dependent on phase
	


figure

	pl = plot(xax, abs(R1.kp_mask)');
	colorlines_gradient(pl,'r'), hold on
	plot(xax, abs(R1.kp_mask)','r','linewidth',2)
	title('kuramoto parameter (mask)')


figure
	pl = plot(xax,abs(R1.kp_othr)','color',[1 .9 .9]); hold on
	colorlines_gradient(pl,'r')
	plot(xax,abs(R1.kp_othr)','r','linewidth',2)
	axis tight
	ylim([0 1])
	title('kuramoto parameter (other)')

	
figure

	plot(xax, R1.meanVm_mask','b','linewidth',1), hold on
	plot(xax, R1.meanVm_othr','b:','linewidth',1)
	plot(xax, R1.meanVm_neig','b-.','linewidth',1)
	plot(xax, R1.meanVm_notm','b--','linewidth',1)

	legend({ 'mask WT' 'other WT' 'neigh WT' 'notmask WT'})


figure

	
	var1 = R1.Vm_mask;
	var2 = R1.Vm_neig;
	var3 = R1.Vm_othr;
	
	Q1  = quantile(var1 , [0.05 .5 .95],1);
	Q2  = quantile(var2 , [0.05 .5 .95],1);
	Q3  = quantile(var3 , [0.05 .5 .95],1);

	fill_between_lines(xax,Q1(1,:),Q1(3,:), 'r'), alpha(.3), hold on
	fill_between_lines(xax,Q2(1,:),Q2(3,:), 'b'), alpha(.3)
	fill_between_lines(xax,Q3(1,:),Q3(3,:), 'g'), alpha(.3)

figure
	plot_mean_and_std(VVmask)
	plot_mean_and_std(VVneig)
	plot_mean_and_std(VVothr)

	% pl = plot(xax,var1','color',[1 .9 .9]); hold on
	% colorlines_gradient(pl,'r')
	% pl = plot(xax,var2','color',[.9 .9 1]); hold on
	% colorlines_gradient(pl,'b')


figure
	plot(xax,mean(var1)','r','linewidth',2)
	plot(xax,mean(var2)','b','linewidth',2)
	plot(xax,mean(var3)','g','linewidth',2)
	legend({ 'mask WT' 'neigh WT' 'other WT' })
	




out.R1 = R1;
out.hilbert = H1;