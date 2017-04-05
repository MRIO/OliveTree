% Jochen_Xcorr_comparisons.m


load periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_spont_50000_4_17-Jan-2017.mat
js{1} = joinsim(simresults, [1:4])
js{2} = joinsim(simresults, [5:8])

load periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017.mat
js{3} = joinsim(simresults, [1:4])
js{4} = joinsim(simresults, [5:8])


neighbors = retrieveNeuronsByClass(js{2}, 'neighbors');
nextneighbors = retrieveNeuronsByClass(js{2}, 'nextneighbors');
stimulated= retrieveNeuronsByClass(js{3}, 'stimulated');
stimulated= stimulated(1:10);
nextneighbors= nextneighbors(1:10);


xcorr_summa(js{1}, 'selectedneurons', neighbors)
xcorr_summa(js{1}, 'selectedneurons', stimulated)
xcorr_summa(js{1}, 'selectedneurons', nextneighbors)

xcorr_summa(js{2}, 'selectedneurons', neighbors)
xcorr_summa(js{2}, 'selectedneurons', stimulated)
xcorr_summa(js{2}, 'selectedneurons', nextneighbors)

xcorr_summa(js{3}, 'selectedneurons', neighbors)
xcorr_summa(js{3}, 'selectedneurons', stimulated)
xcorr_summa(js{3}, 'selectedneurons', nextneighbors)

xcorr_summa(js{4}, 'selectedneurons', neighbors)
xcorr_summa(js{4}, 'selectedneurons', stimulated)
xcorr_summa(js{4}, 'selectedneurons', nextneighbors)

