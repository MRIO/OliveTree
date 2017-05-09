%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% % Jochen_Xcorr_comparisons.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotstyle = '5x5';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% load 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
nonstimulated 	= setdiff(randperm(200,10), stimulated);
stimulated 		= stimulated(1:10);
nextneighbors 	= nextneighbors(1:10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% xcorrs and plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
% compare neighbor cells with and without gap junctions
XC_stim_WT = xcorr_summa(JS{1}, 'selectedneurons', stimulated);
saveallfigs('prefix', 'xcorr_stim_WT', 'style', plotstyle)
close all

XC_stim_MT = xcorr_summa(JS{2}, 'selectedneurons', stimulated);
saveallfigs('prefix', 'xcorr_stim_MT', 'style', plotstyle)
close all
end

if 1
nonstim = [2    10    19    76    96   108   127   131   140];
XC_nostim_WT = xcorr_summa(JS{1}, 'selectedneurons', nonstim);
% saveallfigs('prefix', 'xcorr_nonstim_WT', 'style', plotstyle)
close all

nonstim = [2    10    19    76    96   108   127   131   140];
XC_nostim_MT = xcorr_summa(JS{2}, 'selectedneurons', nonstim);

% saveallfigs('prefix', 'xcorr_nonstim_MT', 'style', plotstyle)
close all


XC_NEIG_WT = xcorr_summa(SPONT_WT, 'selectedneurons', neighbors);
% saveallfigs('prefix', 'xcorr_neighbors_spont_WT', 'style', plotstyle)
close all

XC_NEIG_MT = xcorr_summa(SPONT_MT, 'selectedneurons', neighbors);
% saveallfigs('prefix', 'xcorr_neighbors_spont_MT', 'style', plotstyle)
close all

end

if 0

% X(1,:) = zscore(.XcorrNoAc);
% X(2,:) = zscore(.XcorrNoAc);
% plot(X')


XZ2(1,:) = zscore(XC{3}.XcorrNoAc);
XZ2(2,:) = zscore(XC{4}.XcorrNoAc);
XZ2(3,:) = zscore(XC{5}.XcorrNoAc);
plot(X2')
legend({'stim MT' 'stim WT' 'spont WT'})

end

% gap = 0.05;  noisesig = .3; noiseamp = -.3 ; tau = 30; sametoall = 0.0; simtype = 'spont'; conntype = 'iso' ;  gapcomp = 0;

cmap = flipud(cbrewer('qual', 'Paired',6));
set(0,'defaultaxescolororder', cmap)

XX(1,:) = sum(XC_stim_WT  .XC{1});
XX(2,:) = sum(XC_nostim_WT.XC{1});
XX(3,:) = sum(XC_stim_MT .XC{1});
XX(4,:) = sum(XC_nostim_MT.XC{1});
XX(5,:) = sum(XC_NEIG_WT .XC{1});
XX(6,:) = sum(XC_NEIG_MT .XC{1});


plot([-400:400], XX', 'linewidth',3)
axis tight
xlabel('lag (ms)')
ylabel('sum of correlations for selected cells (N=10)')
legend({'stim WT; corr = .15' 'no stim WT; corr = .15' 'stim MT; corr = .15' 'no stim MT; corr = .15' 'NEIGH WT; corr = 0' 'NEIGH MT; corr = 0'})

