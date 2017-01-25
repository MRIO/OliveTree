% sanity_check_periodic_ampa.m



addpath('/Users/M/Synced/Titan/Bench2/periodic_ampa/')
addpath('/Users/M/Synced/Titan/Bench2/')
addpath('/Users/M/Synced/Titan/Bench/')


% F1 = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_12-Jun-2016.mat';
% F2 = 'periodic_ampa_2_iso_0.04_spont_50000_2_12-Jun-2016.mat';
% F3 = 'periodic_ampa_2_iso_0.04_gallop_50000_2_12-Jun-2016.mat';
% F = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_20-Jun-2016.mat';
% F = 'periodic_ampa_2_iso_0.04_gallop_50000_2_20-Jun-2016.mat';
% F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_1Hz_50000_2_28-Jun-2016.mat'
% F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_gallop_50000_2_29-Jun-2016.mat'

% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_gallop_50000_4_25-Sep-2016.mat';
F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0_1Hz_50000_4_25-Sep-2016.mat';


load('periodic_ampa_replay_06_12_16_with_spont_gaptest2_iso_spont_5000_1_');
NEW = simresults;
profile_sim(NEW{1})
profile_sim(NEW{2})
profile_sim(NEW{3})

load('periodic_ampa_2_iso_0.04_spont_50000_2_12-Jun-2016.mat');
OLD = simresults;
