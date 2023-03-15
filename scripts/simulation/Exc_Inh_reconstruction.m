% Produces simulations and plots for figure 8 of Loyola et al.
% Exc Inh reconstruction.m
% 
% set(0,'DefaultFigureColor', [1 1 1])
set(0, 'DefaultAxesColormap', cbrewer('div', 'Spectral',30))

% [================================================]
% 		 simulation parameters
% [================================================]
%% clear

intervals = [-180:10:180];
% intervals = [-180 100 120 180];
onset_of_inh = 500;

steadystate_time = 1000; %ms
simtime  = 1500; %ms
delta = .025;
gpu = 1;

if exist('seed') ; seed = seed +1 ; else ; seed = 0; end
thisseed = rng(int8(seed),'twister') % random seed only for simulations (not for network cells)
thisseed.Seed

%% [================================================]
% steps to perform
% [================================================]

conjuctive_stimulation = 1;

produce_plots = 1;

% [================================================]
% variables to report
% [================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

selection = {'V_soma', 'I_cx36', 'Hcurrent_q' ,'Calcium_l', 'Calcium_r'};

% variables to store
to_report = vsoma;


% [================================================]
% 		 input
% [================================================]

% currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2
% ounoise_params = [.2 .3 0 5];
% ounoise_params = [1/50 .1 0 double(thisseed.Seed)];
ounoise_params = [0 0 0 double(thisseed.Seed)];
sametoall = 0.05;

% [================================================]
% 		connectivity
% [================================================]
 
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize,normalize)


% Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.
nconns_curlies = 10;
nconns_bridges = 10;
% cells_in_cluster = 20;  % upperbound estimate from: 
gap_curlies = .03;
gap_bridges = .03;
plotconn = 0;
normalize = 1;


load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
noneurons = length(somatapositions);


%    __________  _   ___   _________________________    ______________  __
%   / ____/ __ \/ | / / | / / ____/ ____/_  __/  _/ |  / /  _/_  __/\ \/ /
%  / /   / / / /  |/ /  |/ / __/ / /     / /  / / | | / // /  / /    \  /
% / /___/ /_/ / /|  / /|  / /___/ /___  / / _/ /  | |/ // /  / /     / /
% \____/\____/_/ |_/_/ |_/_____/\____/ /_/ /___/  |___/___/ /_/     /_/



if not(exist('curlies'))
	% curlies:
	% create a network with distance based connectivity for close by connections
	
	% Parameters from: N. Vrieler, S. Loyola, Y. Yarden-Rabinowitz, J. Hoogendorp, N. Medvedev, T. M. Hoogland, C. I. D. Zeeuw, E. d. Schutter, Y. Yarom, M. Negrello, B. Torben-Nielsen, and M. Y. Uusisaari. Variability and directionality of inferior olive neuron dendrites revealed by detailed 3D characterization of an extensive morphological library. Brain structure & function, 92(4):e52068 – 19, 2019.

	% radius within which we allow connections
	median_soma_distance = 20;
	rad_bri = median_soma_distance * 6;

	brick = createW('3d_reconstruction', [], rad_bri, 1, 1, plotconn, [], nconns_curlies, somatapositions,1,[0 0 0 0]);
	brick_bu = brick;
	brick.W = brick.W*gap_bridges;

end


%% [================================================================]
%            create cells and set conductance parameters
%  [================================================================]
cal_factor = -0.05; %60% intrinsically oscillating cells without noise g=0.03


cell_function = 'vanilla'; % 'devel'

def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized2', 'rng', thisseed); 

% AMPA adjusted to reach 5mV amplitude of hyperpolarization

def_neurons.g_CaL = def_neurons.g_CaL - cal_factor;

def_neurons.gbar_gaba_soma = ones(noneurons,1)*.1;
def_neurons.gbar_gaba_dend = ones(noneurons,1)*1;



%%
% [================================================]
%  Projection Fields of AMPA and GABA
% [================================================]

pert.mask     {1} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-60, 0, 0], 'synapseprobability', .85,'plotme',0);
pert.amplitude{1} = 1;
pert.duration {1} = 1;
pert.type	  {1} = 'ampa_dend';

pert.mask     {2} =  create_input_mask(somatapositions, 'reconstruction','radius',100, 'offset', [-40, 0, 0], 'synapseprobability', .85,'plotme',0); % probability adjusted to make number of cells in cluster match
pert.amplitude{2} = 1;
pert.duration {2} = 5; %ms was the duration of the stimulus in Tycho's paper.
pert.type	  {2} = 'gaba_dend';

% adding somatic gaba for longer hyperpolarization effect
pert.mask     {3} = pert.mask{2};
pert.amplitude{3} = 1;
pert.duration {3} = 5; %ms was the duration of the stimulus in Tycho's paper.
pert.type     {3} = 'gaba_soma';



% Plot the overlapping arborinzation masks
% 
% figure
% scatter3(somatapositions(:,1), somatapositions(:,2), somatapositions(:,3), 200, pert.mask{1} + pert.mask{2}*2 -1,'filled'), axis equal
% 
% cm = [150 147 130 ; 60 35 250  ; 250 35 29 ; 200 35 190]/255;
% colormap(cm)
% alpha(.3)
% 
% % apply some current to check the behavior of the cells
% I_app = [];


%%================================================]
% 		 compute transients/steadystate
%=================================================]

if conjuctive_stimulation
 if ~exist('st_st','var')
	disp('calculating transients')
%%   
     st_st.note = 'brick steady state';
	 st_st = IOnet( 'cell_parameters', def_neurons, ...
	 		'perturbation', [], ...
		   	'networksize', [1 1 noneurons] ,'time',steadystate_time ,'W', brick.W ,'ou_noise', ounoise_params , ...
		   	'to_report', to_report ,'gpu', gpu , ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall, 'displaytext', st_st.note);
	 
end

% to list cell parameters uncomment
% st_st.Plist = Plist;


%%


    i = 0;
    for interval = intervals
        i = i + 1;
        

        onset_of_exc = onset_of_inh + interval;
    
        pert.triggers {1} = onset_of_exc;
        pert.triggers {2} = onset_of_inh:onset_of_inh+60;
        pert.triggers {3} = onset_of_inh:onset_of_inh+60;

        sim{i}.note = ['exc vs inh ' num2str(interval)];
        sim{i} = IOnet( 'cell_parameters', def_neurons, ...
                'perturbation', pert, 'tempState', st_st.lastState, ...
                'networksize', [1 1 noneurons] ,'time',simtime ,'W', brick.W ,'ou_noise', ounoise_params  ,...
                'to_report', to_report ,'gpu', gpu , ...
                'cell_function', cell_function ,'delta',delta,'sametoall', sametoall, ...
                'displaytext', sim{i}.note);

        sim{i}.networkHistory.V_soma = single(sim{i}.networkHistory.V_soma);
        sim{i}.note = ['exc vs inh ' num2str(interval)]

        
    end
   		eval(['save exc_inh_net' num2str(seed) '_'  date ' -v7.3'])
end


if produce_plots


   %% oscillation metrics

    combmask = pert.mask{1} & pert.mask{2};
    mean_ampl = @(x) mean(max(x,[],2)-min(x,[],2));
    mean_max_ampl = @(x) mean(max(x,[],2));



    %% phase distributions

    savemovie = 0;
    animate = 0;
    IOI = onset_of_inh-400:onset_of_inh+400;
    IOI = 200:800;
    IOI = 1:simtime;

    %% 

    for f = 1:length(intervals)
        
        ph_dist{f} = phase_distribution_over_time(sim{f},'duration', IOI,...
            'animate',animate, 'fname', ['phasedist_' num2str(f)], 'savemovie',savemovie, 'group', combgroup, ...
            'frames2print', [100 after2(1) late(end-50)]);
        
        ph_dist{f}.pert = sim{f}.perturbation; 


    end

    %% synchrony and proportion of oscillating cells for group and all cells

    for f = 1:length(intervals)
        
        trig = sort([sim{f}.perturbation.triggers{1}, sim{f}.perturbation.triggers{2}]);
        early  = [1:100];
        after1 = [trig(1):trig(1)+100];
        after2 = [trig(end):trig(end)+100];
        late   = [simtime-200:simtime];
        

        sync_estimate(f,1) = mean(abs(ph_dist{f}.phases.orderparameter{2}(early)));
        sync_estimate(f,2) = mean(abs(ph_dist{f}.phases.orderparameter{2}(after1)));
        sync_estimate(f,3) = mean(abs(ph_dist{f}.phases.orderparameter{2}(after2)));
        sync_estimate(f,4) = mean(abs(ph_dist{f}.phases.orderparameter{2}(late)));
        
        sync_estimate_g(f,1) = mean(abs(ph_dist{f}.phases.orderparameter{1}(early)));
        sync_estimate_g(f,2) = mean(abs(ph_dist{f}.phases.orderparameter{1}(after1)));
        sync_estimate_g(f,3) = mean(abs(ph_dist{f}.phases.orderparameter{1}(after2)));
        sync_estimate_g(f,4) = mean(abs(ph_dist{f}.phases.orderparameter{1}(late)));

        osc_cells{f,1} = count_oscillating_cells(sim{f}, early , .1);
        osc_cells{f,2} = count_oscillating_cells(sim{f}, after1, .1);
        osc_cells{f,3} = count_oscillating_cells(sim{f}, after2, .1);
        osc_cells{f,4} = count_oscillating_cells(sim{f}, late,   .1);

        prop_osc_cells(f,1) = osc_cells{f,1}.proportion;
        prop_osc_cells(f,2) = osc_cells{f,2}.proportion;
        prop_osc_cells(f,3) = osc_cells{f,3}.proportion;
        prop_osc_cells(f,4) = osc_cells{f,4}.proportion;

        prop_osc_cells_g(f,1) = length(intersect(osc_cells{f,1}.oscillating, find(combmask)))/noneurons;
        prop_osc_cells_g(f,2) = length(intersect(osc_cells{f,2}.oscillating, find(combmask)))/noneurons;
        prop_osc_cells_g(f,3) = length(intersect(osc_cells{f,3}.oscillating, find(combmask)))/noneurons;
        prop_osc_cells_g(f,4) = length(intersect(osc_cells{f,4}.oscillating, find(combmask)))/noneurons;

        amplitude_osc(f,1) = mean_ampl(sim{f}.networkHistory.V_soma(:,early));
        amplitude_osc(f,2) = mean_ampl(sim{f}.networkHistory.V_soma(:,after1));
        amplitude_osc(f,3) = mean_ampl(sim{f}.networkHistory.V_soma(:,after2));
        amplitude_osc(f,4) = mean_ampl(sim{f}.networkHistory.V_soma(:,late));

        amplitude_osc_g(f,1) = mean_ampl(sim{f}.networkHistory.V_soma(combmask,early));
        amplitude_osc_g(f,2) = mean_ampl(sim{f}.networkHistory.V_soma(combmask,after1));
        amplitude_osc_g(f,3) = mean_ampl(sim{f}.networkHistory.V_soma(combmask,after2));
        amplitude_osc_g(f,4) = mean_ampl(sim{f}.networkHistory.V_soma(combmask,late));

        amplitude_osc(f,1) = mean_max_ampl(sim{f}.networkHistory.V_soma(:,early));
        amplitude_osc(f,2) = mean_max_ampl(sim{f}.networkHistory.V_soma(:,after1));
        amplitude_osc(f,3) = mean_max_ampl(sim{f}.networkHistory.V_soma(:,after2));
        amplitude_osc(f,4) = mean_max_ampl(sim{f}.networkHistory.V_soma(:,late));

        amplitude_osc_g(f,1) = mean_max_ampl(sim{f}.networkHistory.V_soma(combmask,early));
        amplitude_osc_g(f,2) = mean_max_ampl(sim{f}.networkHistory.V_soma(combmask,after1));
        amplitude_osc_g(f,3) = mean_max_ampl(sim{f}.networkHistory.V_soma(combmask,after2));
        amplitude_osc_g(f,4) = mean_max_ampl(sim{f}.networkHistory.V_soma(combmask,late));
     
    end


    


    %% TRIGGERED RESPONSES OVERVIEW
    figure

    % cmap = jet(length(intervals))
    cmap = cbrewer('div', 'RdYlBu', length(intervals))
    % cmap = flipud(cmap)

    IOI = onset_of_inh-400:onset_of_inh+400;

    for f=1:length(intervals)
        subplot(1,length(intervals),f)
        plot(mean(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},IOI)),'color', cmap(f,:),'linewidth',2), hold on
        title(num2str(sim{f}.perturbation.triggers{1}-onset_of_inh))
        axis off
        xlim([IOI(1), IOI(end)])
        ylim([-75,-40])
        alpha(.3)
    end


    %%
    figure
    for f=1:length(intervals)
        collected(f,:) = mean(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},:));
    end
    waterfall(collected)

    %%
    figure
    imagesc(sim{f}.networkHistory.V_soma)

 
    %% statistics summaries
    f_prepost = figure;
    ax0 = axes;

    plot(ax0,prop_osc_cells(:,[1 3 4])', '-o')
    legend(num2str(intervals'))
    title('proportion of oscillating cells')
    ax0.XTick= [1 2 3];
    ax0.XTickLabel = {'Pre', 'after','late'};
    ax0.ColorOrder = cbrewer('div', 'RdYlBu', length(intervals));


    f_prepost_sync = figure;
    ax1 = axes;
    plot(ax1,sync_estimate', '-o')
    legend(num2str(intervals'))
    title('Synchrony (network)')
    ax1.XTick= [1 2 3];
    ax1.XTickLabel = {'Pre', 'after1', 'after2','late'};
    ax1.ColorOrder = cbrewer('div', 'RdYlBu', length(intervals));
    ylim([0.5,1])


    f_prepost_sync = figure;
    ax11 = axes;
    plot(ax11,sync_estimate_g', '-o')
    legend(num2str(intervals'))
    title('Synchrony (group)')
    ax11.XTick= [1 2 3];
    ax11.XTickLabel = {'Pre', 'after1', 'after2','late'};
    ax11.ColorOrder = flipud(cbrewer('div', 'RdYlBu', length(intervals)));


    f_prepost_amp = figure;
    ax2 = subplot(1,2,1);
    plot(ax2,amplitude_osc', '-o')
    legend(num2str(intervals'))
    title('amplitude of rebound')
    ax2.XTick= [1 2 3];
    ax1.XTickLabel = {'Pre', 'after1', 'after2','late'};
    ax2.ColorOrder = flipud(cbrewer('div', 'RdYlBu', length(intervals)));
    % ax2.YLim = [0 20];


    ax2 = subplot(1,2,2);
    plot(ax2,amplitude_osc_g', '-o')
    legend(num2str(intervals'))
    title('amplitude of rebound (group)')
    ax2.XTick= [1 2 3];
    ax1.XTickLabel = {'Pre', 'after1', 'after2','late'};
    ax2.ColorOrder = flipud(cbrewer('div', 'RdYlBu', length(intervals)));
    % ax2.YLim = [0 20];

  

%%  activity of network
   
    plotvolume = 0;
    savemovie = 0;

    if plotvolume

        for f = [19 31 ]
            sim{f}.networkParameters.coords = somatapositions;

            % plot_volume(V,somatapositions,[755:1:1555])
            % animate_volume_hilbert(V,somatapositions,IOI)
            animate_volume_hilbert(sim{f},[1:1000], savemovie, somatapositions)
            % plot_volume(ph_dist{f}.networkHistory.V_soma,somatapositions,[905:50:1255])

        end
    end


    %% responses to stimulation (detail)
        
    for  f =  1:length(intervals)
        figure
        subplot(2,1,1)
        m = pert.mask{1}+pert.mask{2}*2;
        [v s] = sort(m)
        imagesc(sim{f}.networkHistory.V_soma(s,:),[-65 -50])

        subplot(2,1,2)
        plot_mean_and_std(sim{f}.networkHistory.V_soma(pert.mask{1}&pert.mask{2},IOI),'color', [0 1 0])
        title(sim{f}.note)
        legend({'exc' 'excmean' 'inh' 'inhmean'  'comb' 'combmean'})
        alpha(.2)

    end



end