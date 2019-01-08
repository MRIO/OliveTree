% figure_model_Pascal.m


netsize = [2 10 10];
% netsize = [3 30 30];
	noneurons = prod(netsize);

plotthis  = 0;
rd = 2;
meannoconn = 10;

normleak  = 1;
randomize = 1;
scaling   = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;

seed = 90005;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% create network structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak,'seed',seed);

[X Y Z] = meshgrid(1:netsize(1),1:netsize(2),1:netsize(3));
 X = X(:); Y = Y(:); Z = Z(:);


plotnetstruct(W.W,X,Y,Z,sum(W.W) )
axis off
% plotnetstruct(W,X,Y,Z,idx, varargin)

    % if isempty(varargin)
    %     plotconnections = 1;
    %     plotneurons = 1;
    %     onlynetstruct = 1;
    %     plotsingleneuronconnections = 0;

create_input_mask(netsize, 'dist_to_point', 'radius', 3,'cell_coordinates', W.coords,'projection_center', (netsize+1)/2,'synapseprobability',.5,'plotme',1);



% [================================================]
%  datafiles
% [================================================]
	addpath('/Users/M/Synced/Titan/Bench2/periodic_ampa/')
	addpath('/Users/M/Synced/Titan/Bench2/')
	addpath('/Users/M/Synced/Titan/Bench/')

% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
F2 = 'periodic_ampa_nonoise.mat';

F1 = 'periodic_ampa_replay_06_12_16_with_spont1_iso_0.04_spont_50000_1_.mat'


%=============================gather data==============================%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% sample activity in response to noise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load (F1)
s1 = simresults{1};
load (F2)
s2 = simresults{1};


s1.networkHistory.backgroundnoise = [];
s2.networkHistory.backgroundnoise = [];

figure
replayResults_3(s1, [4800:6300])
% colormap(cmap)
axis off

figure
replayResults_3(s2, [4800:6300])
% colormap(cmap)
axis off

