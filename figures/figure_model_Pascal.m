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


% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat';
F2 = 'periodic_ampa_nonoise.mat';

F1 = 'periodic_ampa_replay_06_12_16_with_spont1_iso_0.04_spont_50000_1_.mat'

addpath('/Users/M/Synced/Titan/Bench2/periodic_ampa/')
addpath('/Users/M/Synced/Titan/Bench2/')
addpath('/Users/M/Synced/Titan/Bench/')


%=============================gather data==============================%

load (F1)
s1 = simresults{1};
load (F2)
s2 = simresults{1};


s1.networkHistory.backgroundnoise = [];
s2.networkHistory.backgroundnoise = [];

figure
replayResults_3(s1, [4800:6300])
colormap(cmap)
axis off

figure
replayResults_3(s2, [4800:6300])
colormap(cmap)
axis off



cmap = [...
    0.3137    0.3137    0.3137;
    0.3307    0.3307    0.3307;
    0.3478    0.3478    0.3478;
    0.3648    0.3648    0.3648;
    0.3818    0.3818    0.3818;
    0.3988    0.3988    0.3988;
    0.4158    0.4158    0.4158;
    0.4329    0.4329    0.4329;
    0.4499    0.4499    0.4499;
    0.4669    0.4669    0.4669;
    0.4839    0.4839    0.4839;
    0.5010    0.5010    0.5010;
    0.5180    0.5180    0.5180;
    0.5350    0.5350    0.5350;
    0.5520    0.5520    0.5520;
    0.5690    0.5690    0.5690;
    0.5861    0.5861    0.5861;
    0.6031    0.6031    0.6031;
    0.6201    0.6201    0.6201;
    0.6371    0.6371    0.6371;
    0.6541    0.6541    0.6541;
    0.6712    0.6712    0.6712;
    0.6882    0.6882    0.6882;
    0.7052    0.7052    0.7052;
    0.7222    0.7222    0.7222;
    0.7392    0.7392    0.7392;
    0.7563    0.7563    0.7563;
    0.7733    0.7733    0.7733;
    0.7903    0.7903    0.7903;
    0.8073    0.8073    0.8073;
    0.8243    0.8243    0.8243;
    0.8414    0.8414    0.8414;
    0.8584    0.8584    0.8584;
    0.8630    0.8630    0.8630;
    0.8675    0.8675    0.8675;
    0.8721    0.8721    0.8721;
    0.8767    0.8767    0.8767;
    0.8812    0.8812    0.8812;
    0.8858    0.8858    0.8858;
    0.8904    0.8904    0.8904;
    0.8949    0.8949    0.8949;
    0.8995    0.8995    0.8995;
    0.9041    0.9041    0.9041;
    0.9086    0.9086    0.9086;
    0.9132    0.9132    0.9132;
    0.9178    0.9178    0.9178;
    0.9223    0.9223    0.9223;
    0.9269    0.9269    0.9269;
    0.9315    0.9315    0.9315;
    0.9360    0.9360    0.9360;
    0.9406    0.9406    0.9406;
    0.9452    0.9452    0.9452;
    0.9498    0.9498    0.9498;
    0.9543    0.9543    0.9543;
    0.9589    0.9589    0.9589;
    0.9635    0.9635    0.9635;
    0.9680    0.9680    0.9680;
    0.9726    0.9726    0.9726;
    0.9772    0.9772    0.9772;
    0.9817    0.9817    0.9817;
    0.9863    0.9863    0.9863;
    0.9909    0.9909    0.9909;
    0.9954    0.9954    0.9954;
    1.0000    1.0000    1.0000];