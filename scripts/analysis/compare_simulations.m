% compare_simulations.m

% stim_triggered_spikes.m

% F1 = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_12-Jun-2016.mat';
% F2 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';





% periodic_ampa_8_cluster_gallop_1000_4_03-Jun-2016.mat
% periodic_ampa_8_iso_gallop_50000_4_31-May-2016.mat
% periodic_ampa_8_iso_spont_50000_4_01-Jun-2016
% periodic_ampa_8_cluster_1Hz_50000_4_03-Jun-2016.mat
% periodic_ampa_8_cluster_gallop_50000_4_03-Jun-2016.mat
% periodic_ampa_8_cluster_spont_50000_4_04-Jun-2016.mat
% F = 'periodic_ampa_8_iso_1Hz_50000_4_02-Jun-2016.mat';

% F1 = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_12-Jun-2016.mat';
% F2 = 'periodic_ampa_2_iso_0.04_spont_50000_2_12-Jun-2016.mat';
% F3 = 'periodic_ampa_2_iso_0.04_gallop_50000_2_12-Jun-2016.mat';
% % F = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_20-Jun-2016.mat';
% % F = 'periodic_ampa_2_iso_0.04_gallop_50000_2_20-Jun-2016.mat';
% % F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_1Hz_+w(2)00_2_28-Jun-2016.mat'
% % F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_gallop_50000_2_29-Jun-2016.mat'

addpath('/Users/M/Projects/Experiements/Olive/model/simresults')
F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
F2 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_gallop_50000_4_25-Sep-2016.mat';
F3 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_spont_50000_4_25-Sep-2016.mat';






% cellselection = [35 105 115 170];
cellselection  = [7 35 105 118 115 170]; %[7 35 55 115]

cellofinterest = cellselection(1);

plotstuff = 1;


%=============================gather data==============================%

if ~exist('Joinedsim')
	load (F1)
	numruns = 4;
	Joinedsim{1}  = joinsim(simresults,[1:numruns]); 

	load (F2)
	numruns = 4;
	Joinedsim{2}  = joinsim(simresults,[1:numruns]); 

	load (F3)
	numruns = 4;
	Joinedsim{3}  = joinsim(simresults,[1:numruns]); 
end


% [=================================================================]
%  compare spikes added
% [=================================================================]
duration_in_s = (Joinedsim{1}.duration/1e3);
hist([Joinedsim{1}.spikes.spikespercell; Joinedsim{2}.spikes.spikespercell]'/duration_in_s,20)
xlabel('cell frequency (Hz)')
ylabel('neurons')

% [=================================================================]
%  compare isi histograms
% [=================================================================]
M1 = cell2mat(Joinedsim{1}.spikes.cellISI);
M2 = cell2mat(Joinedsim{2}.spikes.cellISI);
M3 = cell2mat(Joinedsim{3}.spikes.cellISI);
histsup = linspace(0,5000,500);
H1 = hist(M1,histsup);
H2 = hist(M2,histsup);
H3 = hist(M3,histsup);
plot(histsup, [H1; H2 ; H3]');

figure

M1 = Joinedsim{1}.spikes.cellISI{cellofinterest};
M2 = Joinedsim{2}.spikes.cellISI{cellofinterest};
M3 = Joinedsim{3}.spikes.cellISI{cellofinterest};
histsup = linspace(0,5000,100);
H1 = hist(M1,histsup);
H2 = hist(M2,histsup);
H3 = hist(M3,histsup);
plot(histsup, [H1; H2 ; H3]');
xlabel('ms')
ylabel('count')
title('ISI stim x spont')

netfiring(1) = sum(Joinedsim{1}.spikes.spikespercell)/((Joinedsim{1}.duration/1e3)*prod(Joinedsim{1}.networksize));
netfiring(2) = sum(Joinedsim{2}.spikes.spikespercell)/((Joinedsim{2}.duration/1e3)*prod(Joinedsim{2}.networksize));
netfiring(3) = sum(Joinedsim{3}.spikes.spikespercell)/((Joinedsim{3}.duration/1e3)*prod(Joinedsim{2}.networksize));

sum(Joinedsim{1}.spikes.spikespercell)-sum(Joinedsim{2}.spikes.spikespercell);

Joinedsim{1}.spikes.spikespercell(cellofinterest)/(Joinedsim{1}.duration/1e3)
Joinedsim{2}.spikes.spikespercell(cellofinterest)/(Joinedsim{2}.duration/1e3)
Joinedsim{3}.spikes.spikespercell(cellofinterest)/(Joinedsim{3}.duration/1e3)


w = [-1000 1000];

findearlyspike  = @(c) c(find(c>=0 & c<=50, 1, 'first'));
isemptycell     = @(c) ~ isempty(c);
firsttrigoffset = 15;
trigger = 1;


for i = cellselection;

	figure
	triggers1 = Joinedsim{1}.perturbation.triggers{trigger}-firsttrigoffset; % Subtract offset to shift triger to first ampa pulse	
	allspikes1 = Joinedsim{1}.spikes.spikes;
	trigrast1 = ETR(triggers1, allspikes1{i} , 'waves', Joinedsim{1}.networkHistory.V_soma(i,:),'bin', 10, 'span', w(2),'plotQ',plotstuff);
	title([num2str(i) ' 1Hz'])

	early = cellfun(findearlyspike,trigrast1.eventTriggeredRaster,'uniformoutput',0);
	early = find(cellfun(isemptycell, early));
	late  = setdiff([1:length(triggers1)], early)';

	figure
	plot_mean_and_std([w(1):w(2)],trigrast1.eventTriggeredWaveforms(:,early)', 'color', [1 0 0]); 
	hold on
	plot_mean_and_std([w(1):w(2)],trigrast1.eventTriggeredWaveforms(:,late)', 'color', [0 0 1]);
	alpha(.7)
	title('1Hz')

	
	triggers2 = Joinedsim{2}.perturbation.triggers{trigger}-firsttrigoffset; % Subtract offset to shift triger to first ampa pulse	
	allspikes2 = Joinedsim{2}.spikes.spikes;

	figure
	trigrast2 = ETR(triggers2, allspikes2{i} , 'waves', Joinedsim{2}.networkHistory.V_soma(i,:),'bin', 10, 'span', w(2),'plotQ',plotstuff);
	title([num2str(i) ' gallop'])

	early = cellfun(findearlyspike,trigrast2.eventTriggeredRaster,'uniformoutput',0);
	early = find(cellfun(isemptycell, early));
	late  = setdiff([1:length(triggers1)], early)';

	figure
	plot_mean_and_std([w(1):w(2)],trigrast2.eventTriggeredWaveforms(:,early)', 'color', [1 0 0]); 
	hold on
	plot_mean_and_std([w(1):w(2)],trigrast2.eventTriggeredWaveforms(:,late)', 'color', [0 0 1]);
	alpha(.7)
	title('gallop')
	legend({'early' 'late'})
	


	
	triggers = Joinedsim{2}.perturbation.triggers{1}-firsttrigoffset; % Subtract offset to shift triger to first ampa pulse	
	triggers3 = triggers(1:2:end);
	triggers4 = triggers(2:2:end);

	figure
	trigrast3 = ETR(triggers3, allspikes2{i} , 'waves', Joinedsim{2}.networkHistory.V_soma(i,:),'bin', 10, 'span', w(2),'plotQ',plotstuff);
	title([num2str(i) ' gallop (short)'])

	figure
	trigrast4 = ETR(triggers4, allspikes2{i} , 'waves', Joinedsim{2}.networkHistory.V_soma(i,:),'bin', 10, 'span', w(2),'plotQ',plotstuff);
	title([num2str(i) ' gallop (long)'])


	early = cellfun(findearlyspike,trigrast3.eventTriggeredRaster,'uniformoutput',0);
	early = find(cellfun(isemptycell, early));
	late  = setdiff([1:length(triggers1)], early)';

	figure
	plot_mean_and_std([w(1):w(2)],trigrast3.eventTriggeredWaveforms(:,early)', 'color', [1 0 0]); 
	hold on
	plot_mean_and_std([w(1):w(2)],trigrast3.eventTriggeredWaveforms(:,late)', 'color', [0 0 1]);
	alpha(.7)
	title(['gallop (short), cell: ' num2str(i) ])


	early = cellfun(findearlyspike,trigrast3.eventTriggeredRaster,'uniformoutput',0);
	early = find(cellfun(isemptycell, early));
	late  = setdiff([1:length(triggers1)], early)';

	figure
	plot_mean_and_std([w(1):w(2)],trigrast4.eventTriggeredWaveforms(:,early)', 'color', [1 0 0]); 
	hold on
	plot_mean_and_std([w(1):w(2)],trigrast4.eventTriggeredWaveforms(:,late)', 'color', [0 0 1]);
	alpha(.7)
	title(['gallop (long), cell: ' num2str(i) ])


	saveallfigs('prefix', ['trials_with_and_without_spikes_cell_' num2str(i)])
	close all

end



