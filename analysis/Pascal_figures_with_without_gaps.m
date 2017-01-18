% analyze_clusters_bridges.m

calculatesynchrony = 1;

tslice = 1001:4000;

if not(exist('sim'))
	load('periodic_ampa_replay_06_12_16_with_spont_gaptest2_iso_spont_5000_1_')
end


statevar{1} = simresults{1}.networkHistory.V_soma(:,tslice);
statevar{2} = simresults{2}.networkHistory.V_soma(:,tslice);

% replayResults_clusters(sim{1});
% replayResults_clusters(sim{1});

R{1} = profile_sim(simresults{1},'tslice',tslice);
R{2} = profile_sim(simresults{2},'tslice',tslice);

set(0,'defaultaxescolororder', linspecer(10))
set(0,'defaultfigurecolormap', linspecer(10))



freqbins = [0:1:30];
freqhistwithout  = hist(table2array(R{1}.allneurons(:,'freq_each')), freqbins)
freqhistwith  = hist(table2array(R{2}.allneurons(:,'freq_each')), freqbins)


figure
bar(freqbins, [ freqhistwithout ; freqhistwith]')
xlabel('Amplitude (mV) ')
ylabel('Cells')
title('STO freq')

ampbins = [0:1:30];
amplitudehistwithout = hist(table2array(R{1}.allneurons(:,'ampl')),ampbins)
amplitudehistwith 	 = hist(table2array(R{2}.allneurons(:,'ampl'))  ,ampbins)

figure
bar(ampbins, [amplitudehistwithout ; amplitudehistwith]')
xlabel('Amplitude (mV) ')
ylabel('Cells')
title('STO amplitude')

% relationship between cluster amplitude and oscillator amplitude

% relationship between single cell frequency and oscilltor frequency

% probability of non oscillating cluster becoming oscillating cluster

% cluster synchrony distribution with and without bridges

	% mean and std of 2000s of order parameter


% activity of bridge cells and neighbors


% show time for group and global to align



% [=================================================================]
%  Bridges and Curlies
% [=================================================================]


% [=================================================================]
%  Cluster Activity Stats
% [=================================================================]



if calculatesynchrony
	fig3 = figure;;
	
		sim1_sync = measureGroupSync(simresults{1},'plotme',1);
		sim2_sync = measureGroupSync(simresults{2},'plotme',1);
end



