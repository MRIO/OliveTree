% Exc Inh Duel.m
% 
% set(0,'DefaultFigureColor', [1 1 1])
set(0, 'DefaultAxesColormap', cbrewer('div', 'Spectral',30))

% [================================================]
% 		 simulation parameters
% [================================================]
% clear

steadystate_time = 1000; %ms
simtime  = 3000; %ms
delta = .025;
gpu = 1;

if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
thisseed.Seed

% [================================================]
% simulations to perform
% [================================================]

conjuctive_stimulation = 1;


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

nconns_curlies = 5;
nconns_bridges = 5;
cells_in_cluster = 20;  % upperbound estimate from: % Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.
gap_curlies = .05;
gap_bridges = .05;
plotconn = 1;
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


if not(exist('curlies'))
	% curlies:
	% create a network with distance based connectivity for close by connections
	% this network is clusterized with about 20cells per cluster, according to a k-means algo.

	% Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.

	% radius within which we allow connections
	% between curlies:
	median_soma_distance = 20;
	rad_cur = median_soma_distance * 6;
	% between bridges:
	rad_bri = median_soma_distance * 12;


% 	
% 
% 	curlies = createW('3d_reconstruction', [], rad_cur, 1, 0, plotconn, [], nconns_curlies, somatapositions,1,[1 cells_in_cluster 1 0]);
% 
% 	% create a network with distance based connectivity for further apart cells: bridges
% 	% these cells are not bound to specific clusters.
% 	bridges = createW('3d_reconstruction', [], rad_bri, 1, 0, plotconn, [], nconns_bridges, somatapositions,1,[1 cells_in_cluster 0 1]);
% 
% 
% 	% define the indices of 10% of the cells, these will be bridges
% 	bc = randperm(noneurons); % randomly permute cell indices
% 	z = zeros(noneurons,1) ; % initialize index vector
% 	z(bc(1:round(.1*noneurons))) = 1; % make 10% of the cells == bridges
% 	bc =z;
% 	bridge_idx = find(bc);
% 
% 	% remove from curlie adjacency matrix all of those that will become bridges
% 	curlies.W = bsxfun(@times, curlies.W, ~z);
% 	curlies.W = bsxfun(@times, curlies.W, ~(z')); % multiply by the 'unitary' conductance
% 	curlies_bu = curlies; %_bu -> binary undirected
% 	curlies.W = curlies.W*gap_curlies;
% 	% curlies.stats.clusters(bridge_idx) = 0;
% 	
% 	cstats = connectivity_statistics(bridges);
% 	curlies.stats = cstats.stats ;
% 
% 	% remove connections from curlies to bridges from bridge adjacency matrix
% 	bridges.W = bsxfun(@times, bridges.W, z);
% 	% create bridge cells connectivity 
% 	bridges.W = (bridges.W+bridges.W');
% 	bridges_bu = bridges; % _bu -> binary undirected
% 	bridges.W = bridges.W*gap_bridges;
% 
% 	bstats = connectivity_statistics(bridges);
% 	bridges.stats = bstats.stats ;
% 
% 	bridg_curlies.coords = curlies.coords;
% 
% 	bridg_curlies.W = curlies.W + bridges.W;
% 	bridg_curlies.stats = connectivity_statistics(bridg_curlies);
% 	bridg_curlies.stats.clusters = curlies.stats.clusters;
% 	% bridg_curlies.stats.clusters(bridge_idx) = 0;
% 
% 	clusteridx = bridg_curlies.stats.clusters;
% 	clusteridx(logical(bc)) = 70;
% 	plotnetstruct(bridg_curlies.W, bridg_curlies.coords(:,1), bridg_curlies.coords(:,2), bridg_curlies.coords(:,3), clusteridx)

	brick = createW('3d_reconstruction', [], rad_bri, 1, 1, plotconn, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
	brick_bu = brick;
	brick.W = brick.W*gap_curlies;

end



% Wcluster150 = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], plotconn, [1 150 .9 .01],1);


% [=================================================================]
%  create cells
% [=================================================================]

cell_function = 'vanilla'; % 'devel'

def_neurons = createDefaultNeurons(noneurons,'celltypes','random_norm', 'rng', thisseed) 

def_neurons.gbar_gaba_dend = def_neurons.gbar_gaba_dend + 0.75; % subthreshold
def_neurons.g_CaL = def_neurons.g_CaL + .1;

% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% ounoise_params = [.2 .3 0 5];
ounoise_params = [0 0 0 double(thisseed.Seed)];
sametoall = 0.05;


% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]


pert.mask     {1} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-50, -50, 0], 'synapseprobability', 1,'plotme',0)
pert.amplitude{1} = 1;
pert.duration {1} = 1;
pert.type	  {1} = 'ampa_dend';

pert.mask     {2} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-50, 50, 0], 'synapseprobability', .85,'plotme',0) % probability adjusted to make number of cells in cluster match
pert.amplitude{2} = 1
pert.duration {2} = 1;
pert.type	  {2} = 'gaba_dend';



scatter3(somatapositions(:,1), somatapositions(:,2), somatapositions(:,3), 200, pert.mask{1} + pert.mask{2}*2 -1,'filled'), axis equal


cm = [150 147 130 ; 250 35 29 ; 60 35 250 ; 200 35 190]/255;
colormap(cm)




% apply some current to check the behavior of the cells
I_app = [];




% [===========================================================================================================]
 
%%

%%================================================]
% 		 compute transients/steadystate
%=================================================]
s = 0;
neurons = def_neurons;
calfactors = [-.1 0 .1];
for calfact = calfactors
    s = s+1;
	neurons.g_CaL = def_neurons.g_CaL + calfact;
    
	 cal_sim{s} = IOnet( 'cell_parameters', neurons, ...
	 		'perturbation', [], ...
		   	'networksize', [1 1 noneurons] ,'time',steadystate_time ,'W', brick.W ,'ou_noise', ounoise_params , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 cal_sim{s}.note = 'calciumL factors';
end

%
%% oscillating cells in network

figure
for i = 1:length(calfactors)

subplot(3,1,i)
osc_cells{i} = count_oscillating_cells(cal_sim{i},[500:1000], -1)
osc_cells{i}.CaL_factor  = calfactors(i);
title({['Ca L factor: ' num2str(calfactors(i))] ; ['prop osc:' num2str(osc_cells{i}.proportion)]})
ylabel('cells')
end
xlabel('log amplitude (log(mV))')

figure
cmap = flipud(cbrewer('div', 'Spectral',30));
for i = 1:length(calfactors)

subplot(3,1,i)
imagesc(cal_sim{i}.networkHistory.V_soma, [-65 -38]); colorbar; colormap(cmap)
title({['Ca L factor: ' num2str(calfactors(i))] ; ['prop osc:' num2str(osc_cells{i}.proportion)]})
ylabel('cells')
end
xlabel('log amplitude (log(mV))')


%%


%%================================================]
% 		 compute transients/steadystate
%=================================================]
 if ~exist('st_st','var')
	disp('calculating transients')

	 st_st = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], ...
		   	'networksize', [1 1 noneurons] ,'time',steadystate_time ,'W', brick.W ,'ou_noise', ounoise_params , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);
	 st_st.note = 'brick steady state' 
end

% 	% st_st.Plist = Plist;
% end

%%
if conjuctive_stimulation

    i = 0;
    for interval = [-150:25:150]
        i = i + 1;
        
        onset_of_exc = 1000;
        onset_of_inh = onset_of_exc + interval;

        pert.triggers {1} = onset_of_exc;
        pert.triggers {2} = onset_of_inh;

        sim{i} = IOnet( 'cell_parameters', def_neurons, ...
                'perturbation', pert, 'tempState', st_st.lastState, ...
                'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', ounoise_params  ,...
                'to_report', to_report ,'gpu', gpu , ...
                'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

        sim{i}.networkHistory.V_soma = single(sim{i}.networkHistory.V_soma);
        sim{i}.networkHistory.I_cx36 = [];
        sim{i}.networkHistory.backgroundnoise = [];
        sim{i}.note = ['exc vs inh ' num2str(interval)]

    end
   
%%
    
ff = figure
savemovie = 0;
animate = 0;
for f = 1:5
    combgroup = find(pert.mask{1}&pert.mask{2});
    ph_dist{f} = phase_distribution_over_time(sim{f},'duration', [500:1500],'animate',animate, 'fname', ['phasedist_' num2str(f)], 'savemovie',savemovie, 'group', combgroup')
    ph_dist{f}.pert = sim{f}.perturbation;
end

%%
for f = 1:5
    plot_volume(sim{f}.networkHistory.V_soma,somatapositions,[500:2000])
    close all 
end


%%

		eval(['save exc_inh_net' num2str(seed) '_'  date ' -v7.3'])

end
%%
    
for f = 1:5
    figure
    subplot(2,1,1)
    m = pert.mask{1}+pert.mask{2}*2;
    [v s] = sort(m)
    imagesc(sim{f}.networkHistory.V_soma(s,:),[-70 -40])

    subplot(2,1,2)
    plot_mean_and_std(sim{f}.networkHistory.V_soma(pert.mask{1},:)), hold on
    plot_mean_and_std(sim{f}.networkHistory.V_soma(pert.mask{2},:),'color', [0 0 1])
    plot_mean_and_std(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},:),'color', [0 1 0])
    title(sim{f}.note)
    legend({'exc' 'excmean' 'inh' 'inhmean'  'comb' 'combmean'})
    alpha(.5)

end

%% TRIGGERED RESPONSES
figure
cmap = jet(5)
for f=1:5
    subplot(1,5,f)
    plot_mean_and_std(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},800:1500),'color', cmap(f,:)), hold on
    title(num2str(sim{f}.perturbation.triggers{1}-sim{f}.perturbation.triggers{2}))
    ylim([-75,-40])
end


%%
figure
for f = 1:5
    collected(f,:) = mean(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},:));
end
waterfall(collected)

%%

imagesc(collected)


%% 
combmask = pert.mask{1}&pert.mask{2};
for f = 1:5
    
    last_stim  = max([sim{f}.perturbation.triggers{1},sim{f}.perturbation.triggers{2}]);
    
    osc_cells_pre_stim = count_oscillating_cells(sim{f}, [500:600], -1);
    osc_cells_pos_stim = count_oscillating_cells(sim{f}, [last_stim+50:last_stim+150], -1);
    osc_cells_late_stim = count_oscillating_cells(sim{f}, [last_stim+700:last_stim+900], -1);
    
    pre_stim(f)  = osc_cells_pre_stim.proportion;
    post_stim(f) = osc_cells_pos_stim.proportion;
    late_stim(f) = osc_cells_late_stim.proportion;
    
    pre_stim_amp(f,:)  = max(sim{f}.networkHistory.V_soma(:,200:300),[],2);
    post_stim_amp(f,:) = max(sim{f}.networkHistory.V_soma(:,last_stim:last_stim+100),[],2);
    late_stim_amp(f,:) = max(sim{f}.networkHistory.V_soma(:,1800:1900),[], 2);
    
 
end

f_prepost = figure;
ax = axes;

data = [pre_stim', post_stim', late_stim']';
plot(ax,[pre_stim', post_stim', late_stim']', '-o')
legend(num2str([-125 -68 0 68 125]'))
title('proportion of oscillating cells')
ax.XTick= [1 2 3];
ax.XTickLabel = {'Pre', 'Post', 'Late'}
ax.ColorOrder = cbrewer('qual', 'Set1', 5);

%%

f_prepost_amp = figure;
ax2 = subplot(1,2,1);

amplitudes = [mean(pre_stim_amp,2), mean(post_stim_amp,2), mean(late_stim_amp,2)]';
plot(ax2,amplitudes, '-o')
legend(num2str([-125 -68 0 68 125]'))
title('amplitude of rebound')
ax2.XTick= [1 2 3];
ax2.XTickLabel = {'Pre', 'Post', 'Late'}
ax2.ColorOrder = cbrewer('qual', 'Set1', 5);
ax2.YLim = [-60 -40];

ax3 = subplot(1,2,2);
amplitudes = [mean(pre_stim_amp(:,combmask),2), mean(post_stim_amp(:,combmask),2), mean(late_stim_amp(:,combmask),2)]';
plot(ax3,amplitudes, '-o')
legend(num2str([-125 -68 0 68 125]'))
title('amplitude of rebound (stimulated cells)')
ax3.XTick= [1 2 3];
ax3.XTickLabel = {'Pre', 'Post', 'Late'}
ax3.ColorOrder = cbrewer('qual', 'Set1', 5);
ax3.YLim = [-60 -40];

    