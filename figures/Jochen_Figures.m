%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Figures for XCorr paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure:

% 	- Summary of effect of coupling on network of heterogeneous components
% 	- Summary of simulation and XCorrs for coupled and uncoupled networks under stimulation
% 	- Summary of Xcorr for conditions

% Suppl: 
% 
% 	- Gap compensation parameter space


% Compares the cross correlations of coupled and uncoupled networks
% for different cell selections, including 'stimulated' 'neighbors';

plotstyle = '5x5';
saveallfigs = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% load 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/Users/M/Projects/Experiments/Olive/model/simresults/xcorr')
addpath('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/xcorr')

load periodic_ampa_xcorr_stim_tau_30_WT_4_iso_1Hz_50000_4_.mat
 
JS{1} = joinsim(simresults, [1:4])

load periodic_ampa_xcorr_stim_tau_30_MT_4_iso_1Hz_50000_4_.mat 

JS{2} = joinsim(simresults, [1:4])

load periodic_ampa_test_MT_1_iso_spont_50000_1_.mat
SPONT_MT = simresults{1};

load periodic_ampa_test_WT_1_iso_spont_50000_1_.mat
SPONT_WT = simresults{1};

neighbors 		= retrieveNeuronsByClass(JS{1}, 'neighbors');
nextneighbors 	= retrieveNeuronsByClass(JS{1}, 'nextneighbors');
stimulated 		= retrieveNeuronsByClass(JS{1}, 'stimulated');

noneu = [length(neighbors) length(nextneighbors) length(stimulated)]

nonstimulated 	= setdiff([1:200], stimulated);
stimulated 		= stimulated(1:10);
nextneighbors 	= nextneighbors(1:10);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% xcorrs and plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% compare neighbor cells with and without gap junctions
	XC_stim_WT = xcorr_summa(JS{1}, 'selectedneurons', stimulated);
	XC_stim_MT = xcorr_summa(JS{2}, 'selectedneurons', stimulated);
	
	% a static selection of non stimulated cells
	nonstim = [2    10    19    76    96   108   127   131   140];
	XC_nostim_WT = xcorr_summa(JS{1}, 'selectedneurons', nonstim);
	XC_nostim_MT = xcorr_summa(JS{2}, 'selectedneurons', nonstimulated);

	XC_NEIG_WT = xcorr_summa(SPONT_WT, 'selectedneurons', neighbors);
	XC_NEIG_MT = xcorr_summa(SPONT_MT, 'selectedneurons', neighbors);




cmap = flipud(cbrewer('qual', 'Paired',6));
set(0,'defaultaxescolororder', cmap)

XX(1,:) = mean(XC_stim_WT.XC{1});
XX(2,:) = mean(XC_nostim_WT.XC{1});
XX(3,:) = mean(XC_stim_MT.XC{1});
XX(4,:) = mean(XC_nostim_MT.XC{1});
XX(5,:) = mean(XC_NEIG_WT.XC{1});
XX(6,:) = mean(XC_NEIG_MT.XC{1});


plot([-400:400], XX', 'linewidth',3)
axis tight
xlabel('lag (ms)')
ylabel('mean of correlations for selected cells (N=10)')
legend({'stim WT; corr = .15' 'no stim WT; corr = .15' 'stim MT; corr = .15' 'no stim MT; corr = .15' 'NEIGH WT; corr = 0' 'NEIGH MT; corr = 0'})

