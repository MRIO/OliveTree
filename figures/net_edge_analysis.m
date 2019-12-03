% net_edge_analysis.m

	netsize = [2 10 10];
		noneurons = prod(netsize);

	plotthis  = 0;
	rd = 2;
	meannoconn = 8;
	normleak  = 1;
	randomize = 1;
	scaling   = 0.04;
	maxiter	  = 1;
	somatapositions = [];
	randomize = 1;
	symmetrize = 1;

	W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);


	plotnetstruct(W.W,W.coords(:,1),W.coords(:,2),W.coords(:,3),W.stats.clustercoeff.bu)

	borders = W.coords(:,2)==1 | W.coords(:,2)==10 | W.coords(:,3)==1 | W.coords(:,3)==10 

	figure
	BU = W.stats.clustercoeff.bu;
	[h1 x] = hist(BU(borders));
	h2 = hist(BU(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('cluster coefficient')
	legend({'border' 'center'})
	[sig p] = kstest2(BU(borders), BU(~borders))
	title({'cluster coefficient (binary undirected)' ;  ['sig: ' num2str(sig) ' ; p-val: ' num2str(p)]})

	figure
	gapleaks = sum(W.W);
	[h1 x] = hist(gapleaks(borders));
	h2 = hist(gapleaks(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('gap leaks')
	legend({'border' 'center'})
	[sig p] = kstest2(gapleaks(borders), gapleaks(~borders))
	title({'gap leaks'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})
	
	figure
	conns = sum(W.W>0);
	[h1 x] = hist(conns(borders));
	h2 = hist(conns(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('number of connections')
	[sig p] = kstest2(conns(borders), conns(~borders))
	legend({'border' 'center' })
	title({'connections'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})
	
	figure
	medfr = simresults{1}.spikes.medfreq;
	[h1 x] = hist(medfr(borders));
	h2 = hist(medfr(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('median frequency')
	[sig p] = kstest2(medfr(borders), medfr(~borders))
	legend({'border' 'center' })
	title({'median frequency (Hz)'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})

	figure
	firingrate = simresults{1}.spikes.spikespercell/(50*4); % (4 runs of 50s each )
	[h1 x] = hist(firingrate(borders));
	h2 = hist(firingrate(~borders),x);
	bar(x,[h1 ; h2]','stacked')
	xlabel('Firing rate (Hz)')
	[sig p] = kstest2(firingrate(borders), firingrate(~borders))
	legend({'border' 'center' })
	title({'firing rate (Hz)'; ['sig: ', num2str(sig) ' ; p-val: ' num2str(p)]})



