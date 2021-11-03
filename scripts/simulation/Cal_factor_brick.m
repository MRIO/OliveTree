% Sweep cal factor brick
% 
% set(0,'DefaultFigureColor', [1 1 1])
set(0, 'DefaultAxesColormap', cbrewer('div', 'Spectral',30))

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

simtime  = 2000; %ms
delta = .025;
gpu = 1;

if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)

% [================================================]
% variables to report
% [================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

% variables to store
to_report = vsoma;


% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)

% nconns_curlies = 8;
% nconns_bridges = 6;

% gap_curlies = .1;
% gap_bridges = .1;
% plotconn = 1;
% normalize = 1;

nconns = 6;
gap = .05;
plotconn = 0;
normalize = 1;

load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);


%     __    __
%    / /_  / /___ _
%   / __ \/ / __ `/
%  / /_/ / / /_/ /
% /_.___/_/\__,_/


	median_soma_distance = 20;
	rad_bri = median_soma_distance * 6;

cell_function = 'vanilla'; % 'devel'


%%================================================]
% 		 compute transients/steadystate
%=================================================]

nets= 2;
calfactors = [-.3:0.1 :-.2];
cal_sim = cell(length(calfactors),nets)
for n = 1:nets
    
    brick = createW('3d_reconstruction', [], rad_bri, 1, 1, 0, [], nconns, somatapositions,1,[0 0 0 0]);
    W = brick.W*gap;
    thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
    parfor s = 1:length(calfactors)
        neurons = struct();
        thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
        neurons = createDefaultNeurons(noneurons,'celltypes','randomized', 'rng', thisseed);
        neurons.g_CaL = neurons.g_CaL + calfactors(s);
        neurons.seed = thisseed;

         cal_sim{s,n} = IOnet( 'cell_parameters', neurons, ...
                'perturbation', [], ...
                'networksize', [1 1 noneurons] ,'time',simtime ,'W', W  , ...
                'to_report', to_report ,'gpu', gpu , ...
                'cell_function', cell_function ,'delta',delta,'sametoall', 0);
         cal_sim{s,n}.note = ['calciumL factor:' num2str(calfactors) , ' net:' num2str(n)] ;
    end
    seed = seed +1;
end

%
%% oscillating cells in network

figure

c = 0;
for n = 1:nets
    for i = 1:length(calfactors)
    c = c+1;
	subplot(nets,length(calfactors), c)
    osc_cells{i,n} = count_oscillating_cells(cal_sim{i,n},[500:1000], -1)
    osc_cells{i,n}.CaL_factor  = calfactors(i);
    title({['Ca L factor: ' num2str(calfactors(i))] ; ['prop osc:' num2str(osc_cells{i,n}.proportion)]})
    
    
    collected_histograms(c,:) = osc_cells{i,n}.histogram.Values;
    end
end
ylabel('cells')
xlabel('log amplitude (log(mV))')

figure
waterfall(osc_cells{i,n}.histogram.BinEdges(1:end-1), calfactors, collected_histograms(1:11,:))




figure
cmap = flipud(cbrewer('div', 'Spectral',30));
c = 0;
for n = 1:nets
    for i = 1:length(calfactors)
    c = c+1;
    subplot(nets,length(calfactors), c)
    imagesc(cal_sim{i,n}.networkHistory.V_soma, [-65 -38]); colormap(cmap)
    title({['Ca L factor: ' num2str(calfactors(i))] ; ['prop osc:' num2str(osc_cells{i,n}.proportion)]})
    
    end
xlabel('log amplitude (log(mV))')
end
ylabel('cells')
    

%%
figure
for n = 1:nets
for i = 1:length(calfactors)

     osc_pspace(i,n) = osc_cells{i,n}.proportion;
     xtick(i) = osc_cells{i,n}.CaL_factor;
     
end
end

plot(xtick, osc_pspace,'o-')
title('prop osc cells')
xlabel('Ca L')
ylabel('proportion')
legend({'net 1' ; 'net 2'})