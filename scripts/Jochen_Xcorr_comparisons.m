% Jochen_Xcorr_comparisons.m


load periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_spont_50000_4_17-Jan-2017.mat

js{1} = joinsim(simresults, [1:4])
js{2} = joinsim(simresults, [5:8])

load periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017.mat
js{3} = joinsim(simresults, [1:4])
js{4} = joinsim(simresults, [5:8])


neighbors = retrieveNeuronsByClass(js{2}, 'neighbors');
nextneighbors = retrieveNeuronsByClass(js{2}, 'nextneighbors');
stimulated= retrieveNeuronsByClass(js{2}, 'stimulated');
stimulated= stimulated(1:10);
nextneighbors= nextneighbors(1:10);

% compare neighbor cells with and without gap junctions
XC{1} = xcorr_summa(js{1}, 'selectedneurons', neighbors)
saveallfigs('prefix', 'xcorr_neighbors_spont_MT', 'style', '29x9')
close all

XC{2} = xcorr_summa(js{2}, 'selectedneurons', neighbors)
saveallfigs('prefix', 'xcorr_neighbors_spont_WT', 'style', '12x12')
close all

X(1,:) = zscore(XC{1}.XcorrNoAc);
X(2,:) = zscore(XC{2}.XcorrNoAc);
plot(X')


% compare stimulated neurons before and after stimulation -- also with gaps
XC{3} = xcorr_summa(js{3}, 'selectedneurons', stimulated)
saveallfigs('prefix', 'xcorr_stimulated_stim_MT', 'style', '12x12')
close all

XC{4} = xcorr_summa(js{4}, 'selectedneurons', stimulated)
saveallfigs('prefix', 'xcorr_stimulated_stim_WT', 'style', '29x9')
close all

XC{5} = xcorr_summa(js{2}, 'selectedneurons', stimulated)
close all
saveallfigs('prefix', 'xcorr_stimulated_spont', 'style', '29x9')


XZ2(1,:) = zscore(XC{3}.XcorrNoAc);
XZ2(2,:) = zscore(XC{4}.XcorrNoAc);
XZ2(3,:) = zscore(XC{5}.XcorrNoAc);
plot(X2')
legend({'stim MT' 'stim WT' 'spont WT'})

xcorr_summa(js{3}, 'selectedneurons', neighbors)
xcorr_summa(js{3}, 'selectedneurons', stimulated)
xcorr_summa(js{3}, 'selectedneurons', nextneighbors)

xcorr_summa(js{4}, 'selectedneurons', neighbors)
xcorr_summa(js{4}, 'selectedneurons', stimulated)
xcorr_summa(js{4}, 'selectedneurons', nextneighbors)

