		% under construction
% Jochen_stats.m



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






% [================================================]
%  stats
% [================================================]
simtime = 200;



freq_stim_cells_WT = mean(JS{1}.spikes.spikespercell(stimulated));
freq_neigh_cells_WT = mean(JS{1}.spikes.spikespercell(neighbors));
freq_nostim_cells_WT = mean(JS{1}.spikes.spikespercell(setdiff([1:200], stimulated)));

freq_stim_cells_MT = mean(JS{2}.spikes.spikespercell(stimulated));
freq_neigh_cells_MT = mean(JS{2}.spikes.spikespercell(neighbors));
freq_nostim_cells_MT = mean(JS{2}.spikes.spikespercell(setdiff([1:200], stimulated)));

[freq_stim_cells_WT freq_neigh_cells_WT freq_nostim_cells_WT ;
freq_stim_cells_MT freq_neigh_cells_MT freq_nostim_cells_MT]/simtime

% must compute:
% - probability of response for stimulated/neighbors/nostim cells WT MT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% xcorrs and plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





files = {'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0_MT_4_iso_1Hz_50000_4_0_29-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.15_MT_4_iso_1Hz_50000_4_0_30-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.3_MT_4_iso_1Hz_50000_4_0_30-Aug-2017..mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.45_MT_4_iso_1Hz_50000_4_0_30-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0_WT_4_iso_1Hz_50000_4_0_29-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.15_WT_4_iso_1Hz_50000_4_0_30-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.3_WT_4_iso_1Hz_50000_4_0_30-Aug-2017.mat';
		'periodic_ampa_xcorr_stim_ncorr_pspace_ncorr_0.45_WT_4_iso_1Hz_50000_4_0_31-Aug-2017.mat'};

		ff = 0;
		for f= files'
			ff = ff+1;
			load(f{1});
			JS{ff} = joinsim(simresults, [1:4]);
			stimulated 		= retrieveNeuronsByClass(JS{1}, 'stimulated');
			nonstimulated   = setdiff([1:200], stimulated);
			XC_stim{ff}       = xcorr_summa(JS{ff}, 'selectedneurons', stimulated);
			XC_nostim{ff}     = xcorr_summa(JS{ff}, 'selectedneurons', nonstimulated);





		end

		neighbors 		= retrieveNeuronsByClass(JS{1}, 'neighbors');
		nextneighbors 	= retrieveNeuronsByClass(JS{1}, 'nextneighbors');
		stimulated 		= retrieveNeuronsByClass(JS{1}, 'stimulated');




	% compare neighbor cells with and without gap junctions
	XC_stim_WT = xcorr_summa(JS{1}, 'selectedneurons', stimulated);
	XC_stim_MT = xcorr_summa(JS{2}, 'selectedneurons', stimulated);
	
	% a static selection of non stimulated cells
	nonstim = [2    10    19    76    96   108   127   131   140];
	XC_nostim_WT = xcorr_summa(JS{1}, 'selectedneurons', nonstim);
	XC_nostim_MT = xcorr_summa(JS{2}, 'selectedneurons', nonstim);

	XC_NEIG_WT = xcorr_summa(SPONT_WT, 'selectedneurons', neighbors);
	XC_NEIG_MT = xcorr_summa(SPONT_MT, 'selectedneurons', neighbors);


% repeated measures anova
%  term: experimenter adjusts factors measures responses in an attempt to determine an effect.

% balanced design (all cells, same number of observations)
% fixed-effect model ()
% multiple-factors
% 
% effects: coupling, membership, background correlation


% measures: xcorr peak, baseline firing rate, added spikes


% distribution of median instantaneous firing rates across different correlation conditions



