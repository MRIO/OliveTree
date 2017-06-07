% [results] = ...
% IOnet('parameter', 'value', ...)
% 
% [=================================================================]
%  DESCRIPTION
% [=================================================================]
% 
% simulates an inferior olive !
% 
% Matlab implementation of the model used in:
% Negrello, Sokol extensions (GABA, AMPA) to the model in De Gruijl et al., 2012, PLoS Computational Biology
%
% authors: mnegrello@gmail.com, piotr sokol, jornt de gruijl
%
% [=================================================================]
%  INPUTS
% [=================================================================]
% 
% 'parameter', default   => explanation
% ------------------------------------------------------------------
% 
% 'delta', .05                          => forward euler integration step
% 'time', 1000                          => simulation time in ms
% 'networksize', [10 10]                => in the absence of an adjacency matrix (W), creates a uniform grid up to 3D
% 'g_CaL',  0.7                         => low-threshold calcium conductance (default: 0.7 mS/cm^2) -- columns matrix containing for every cell the
% 'W',  []                              => adjacency matrix weighted by gap junctions
% 'connections','2d neighborhood' (8 neighbors)
% 'tempState',          []              => statevariables from field 'lastState' -- for continuing the simulation where we stopped
% 'gpu',          true   => if we are using the gpu (only works with CUDA)
% 'ou_noise',    [0 0 0 randomseed])  => parameters for the ohrstein uhlenbeck noise generator: [ theta, sigma, mu, seed ]
% 'cell_parameters', []                 => see createDefaultNeurons.m
% 'perturbation',  [])                  => perturbation structure with fields 
%       pert.mask                               == column vector (neurons x 1) 
%       pert.amplitude                          == scalar
%       pert.duration                           == scalar (duration of pulse in ms)
%       pert.onsets                             == integer - stimulation times in ms == [t1 ... tn];
% 
% 'to_report',{'V_soma'}                => fields to report (cell array with strings from below)
%    available fields for reporting: 
%       {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r', ...
%       'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a' , 'g_GABA_soma','g_GABA_dend', 'g_AMPA', 'Ca2_soma'
% 
% 'sametoall', 0                        => value between 0 and 1 representing how much of the noise is shared between neurons
% 
% [=================================================================]
%  OUTPUTS
% [=================================================================]
% 
% results.networkHistory = netHist (v_soma, plus any other fields in 'to_report')
% results.networksize = networksize;
% results.lastState = tempState;
% results.initialState = gather(tempState);
% results.duration = t*delta; ( brain time)
% results.timeElapsed = toc(tstart); (simulation time)
% results.adjacencyMatrix = sparse(gather(W));
% results.dt = delta;
% 
% == see notes in code for TODO's
% TODO: correct the external application of CURRENT
% TODO: write cell parameters in the output 
% TODO: write coords and W to ouptut
% TODO: change parameter name ou_noise to noise_generator
% TODO: write W with coordinates in the output structure
% 
% TODO: (done) try attempt_gpu catch gpu = 0 end
% TODO: (done) introduce parameter for noise noisemixing
% TODO: reintroduce V_app in cell function
% TODO: (done): createDefaultNeurons is not 100% working.
% TODO: (done): Save all state variables in LastState
% TODO: (done): pass gbarAMPA and gbargaba to the cell
% TODO: (done) REFACTORY: change from rows x columns to arbitrary netsize in 3d
% TODO: SOON: (done)introduce excitatory and inhibitory synaptic noise!
% TODO: (done): remove useless parameters
% TODO: 4 MATHWORKS: allow the passing of parameters in a structure to the gpuarray!
% TODO: make delta a global -- tried, doesn't work
% TODO: (done): allow for non gpu
% TODO: (done): introduce probes to_report (e.g. Calcium)
% TODO: (done): recheck tempState variables
% TODO: (done): pass text


function [results] = IOnet(varargin)

p = inputParser;

p.addParamValue('delta',         .05)  % dt
p.addParamValue('time',          100) % ms
p.addParamValue('networksize',   [10 10 1])
p.addParamValue('W',             [])
p.addParamValue('connections',   [])   % 'nearest4'...'yosi'
p.addParamValue('tempState',     [])
p.addParamValue('appCurrent',    [])     %uA/cm^2 -> .1 nA for a cell with 10000um^2
p.addParamValue('appVoltage',    [])  % legacy
p.addParamValue('gpu',          true)  
p.addParamValue('ou_noise',    [0 0 0 0]);  %  ohrstein-uhlenbeck current noise applied to soma [ theta, sigma, mu, seed  ]
p.addParamValue('perturbation',  []);  % [ amplitude, shift, seed  ]
p.addParamValue('cell_parameters', [])
p.addParamValue('cell_function', 'original') % devel, ode, betta
p.addParamValue('to_report',{'V_soma'}) 
p.addParamValue('sametoall',0) % value between 0 and 1 representing how much of the noise is shared between neurons
p.addParamValue('saveappliednoise',1) 
p.addParamValue('displaytext',[]) 
p.addParamValue('debug',[])

p.parse(varargin{:});

delta = p.Results.delta;
time                = p.Results.time;
networksize         = p.Results.networksize;
tempState           = p.Results.tempState;
W                   = p.Results.W; % REPRESENTS GAP JUNCTION WEIGHTS
connections         = p.Results.connections;
ou_noise          = p.Results.ou_noise;
I_app               = p.Results.appCurrent;
V_app               = p.Results.appVoltage;
use_gpu             = p.Results.gpu;
% compensate_for_gapleak = p.Results.compensate_for_gapleak;
pert                = p.Results.perturbation; % includes input masks, a dirac train, pulse duration and amplitude
% pert.mask ,pert.amplitude, pert.triggers, pert.pulseduration
cell_parameters     = p.Results.cell_parameters;
to_report           = p.Results.to_report;
cell_function       = p.Results.cell_function;
sametoall           = p.Results.sametoall;
saveappliednoise    = p.Results.saveappliednoise;
displaytext         = p.Results.displaytext;
debugging           = p.Results.debug;

interrupt_when_fail = 1;
report_all_dt = false;

if saveappliednoise & sum(ou_noise(1:3))
    to_report = [ to_report 'backgroundnoise'];
end

% [================================================]
%          gpu requirements
% [================================================]

if strcmp(cell_function, 'ode'); use_gpu = 0; disp('ode45 not compatible with gpu'); end
if strcmp(cell_function, 'devel'); disp('warning: debugging not compatible with gpu'); end
try 
    if gpuDeviceCount == 0; use_gpu = 0; disp('no gpu found in this machine'); end     
catch E
    use_gpu = 0; disp('no gpu found in this machine');  
        % keyboard
end


%     _       _ __  _       ___             __  _           
%    (_)___  (_) /_(_)___ _/ (_)___  ____ _/ /_(_)___  ____ 
%   / / __ \/ / __/ / __ `/ / /_  / / __ `/ __/ / __ \/ __ \
%  / / / / / / /_/ / /_/ / / / / /_/ /_/ / /_/ / /_/ / / / /
% /_/_/ /_/_/\__/_/\__,_/_/_/ /___/\__,_/\__/_/\____/_/ /_/ 
   disp('Initializing...')                                                       

% clear the persistent function
clear apply_perturbation


noNeurons = prod(networksize);
% produce default weight matrix (with no connections)
if isempty(W)
    if isempty(connections)
        W = zeros(noNeurons);
        warning('no Weight matrix given, we are creating a matrix with no connections.')
    end

    if isstruct(W)
        % Wstruct = W;
        W = W.W;

    end

end
if use_gpu
    W = gpuArray(full(W));
    wgpu = 'withgpu';
else
    wgpu = [];
end

if noNeurons <= 20
    wgpu = 'few neurons, no gpu'
end

Vdiffs = 0;

if use_gpu; datatyp = 'gpuArray'; else; datatyp = 'double'; end
% Init with or without tempstate (according to gpu flag)

state = initNetState(noNeurons, use_gpu, tempState);

simSteps       = ceil(time/delta);
clock_freq = round(1/delta); % this may lead to errors (when 1/delta had periodic decimals)
state.delta_vector = delta*ones(noNeurons,1, datatyp);

g_cx36 = zeros(noNeurons,1,datatyp); % alloc

if isempty(cell_parameters);
   cell_parameters = createDefaultNeurons(noNeurons);
end


if length(to_report)>1
    if report_all_dt
        for ftr = 1:length(to_report)
            eval(['netHist.' to_report{ftr} '= zeros(noNeurons, time*1/delta);']);
        end
    else
        for ftr = 1:length(to_report)
            eval(['netHist.' to_report{ftr} '= zeros(noNeurons, time*1);']);
        end
    end

end

% If we are doing a Voltage clamp (#BROKEN?)
if isempty(V_app) 
    vclamp = zeros(noNeurons,1, datatyp);
elseif numel(V_app) == 1 & ~isempty(V_app);
    vclamp = ones(noNeurons,1, datatyp)*V_app;
else
    vclamp = V_app.*ones(noNeurons,simSteps, datatyp);
end


% [=================================================================]
%  % USE SPECIFIED RANDOM SEED FOR NOISE
% [=================================================================]
rng(ou_noise(4),'twister');

totalcurr_noise = zeros(noNeurons,1); % for logging purposes
tstart  = tic; 
results.failed = 0; 

th = ou_noise(1);  sig = ou_noise(2); mu = ou_noise(3);
mixalpha = sametoall;


% Schweighofer's gap non linearity (2004)
fgap = @(DeltaV) (0.8 .* exp(-1.*DeltaV.*DeltaV/100) + 0.2);
gaps = logical(sum(sum(W)));


% logicals are not gpu compatible in matlab2013a, we cast.
Wlogical = double(W>0);
if use_gpu
    Wlogical = gpuArray(Wlogical);
end



%    _______  ______     _____(_)___ ___ 
%   / ___/ / / / __ \   / ___/ / __ `__ \
%  / /  / /_/ / / / /  (__  ) / / / / / /
% /_/   \__,_/_/ /_/  /____/_/_/ /_/ /_/ 
           

for t = 1:simSteps

    % [================================================]
    %    compute gap currents
    % [================================================]
    % SCHWEIGHOFER 2004 VERSION    
    % update dendritic compartment according to potential differences

    Vdiffs = bsxfun(@minus, state.V_dend', state.V_dend);
    g_cx36 = fgap(Vdiffs).*W ; % gap gain factor per dendritic pair
    state.I_cx36 = sum(g_cx36 .*Vdiffs)'; % gap currents



    % [================================================]
    %    synaptic currents (pert.)
    % [================================================]
    % dirac tstamps ('pert.triggers') trigger release of neurotransmitter according to specified masks ('pert.masks').

    if ~isempty(pert)
        [state.g_ampa_soma state.g_ampa_dend state.g_gaba_dend state.g_gaba_soma state.Ca2_soma state.curr_noise_pert] = ...
            apply_perturbation(pert,t,clock_freq, state.g_ampa_soma, state.g_ampa_dend, state.g_gaba_dend, state.g_gaba_soma, state.Ca2_soma);
    end


    % [================================================]
    %   generate noisy current for all the neurons
    % [================================================]

    if ~ou_noise(2) | sum(ou_noise)>0
        % ohrstein uhlenbeck noise
        % see: (https://math.stackexchange.com/questions/1287634/implementing-ornstein-uhlenbeck-in-matlab)
        % x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;

        % state.ou_noise =  state.ou_noise +th*(mu-state.ou_noise)*delta + ...
        %                     (1-mixalpha)*sig*sqrt(delta)*randn(noNeurons,1) + ...
        %                      mixalpha*sig*sqrt(de§lta)*randn*ones(noNeurons,1);


        state.ou_noise =  state.ou_noise +th*(mu-state.ou_noise)*delta + ...
                            (1-mixalpha)*sig*sqrt(delta)*randn(noNeurons,1) + ...
                             mixalpha*sig*sqrt(delta)*randn*ones(noNeurons,1);


    end


    
    % [================================================]
    %          Apply current  / voltage / noise
    % [================================================]
    if ~isempty(I_app)
    
        if prod(size(I_app))==1
            state.current = state.curr_noise_pert + I_app(:) + state.ou_noise;
            % disp('injecting 1')
        elseif sum(I_app(:,t)~=0)
            state.current = state.curr_noise_pert + I_app(:,t) + state.ou_noise;
            % disp('injecting 2')
        else
            state.current = state.curr_noise_pert + state.ou_noise;
            % disp('injecting 3')
        end
    else
            state.current = state.curr_noise_pert + state.ou_noise;
            % disp('injecting ounoise')

    end

    
    state.backgroundnoise = state.ou_noise + state.curr_noise_pert;

    totalcurr_noise = totalcurr_noise + state.ou_noise + state.curr_noise_pert; % logging total amount of noise applied



    if numel(V_app) == noNeurons*simSteps  % if we received a V_app for every t
            state.vclamp = V_app(:,t);
    end

    % [================================================]
    %          compute io cell
    % [================================================]
    %     ________            ____
    %    /  _/ __ \________  / / /
    %    / // / / / ___/ _ \/ / / 
    %  _/ // /_/ / /__/  __/ / /  
    % /___/\____/\___/\___/_/_/ 


            try
            
                 state = IOcell_wrapper(state, cell_parameters, cell_function); 

            catch E
                
                state
                cell_parameters
                E.stack
                disp('problem computing the selected cell function!')
                keyboard
                
            end

    
    % [================================================]
    %          collect results
    % [================================================]

    
    % log history according fields in 'to_report' cell array
    if report_all_dt

        if length(to_report)>0
            for ftr = 1:length(to_report)
                eval(['netHist.' to_report{ftr} '(:,t) = gather(state.' to_report{ftr} ');']);
            end
        end

    
    elseif ~mod(t,clock_freq) % report every milliseconds
        if length(to_report)>0
            for ftr = 1:length(to_report)
                try
                    eval(['netHist.' to_report{ftr} '(:,t/clock_freq) = gather(state.' to_report{ftr} ');']);
                    % eval(['netHist.' to_report{ftr} '(:,t/clock_freq) = state.' to_report{ftr} ';']);
                catch
                    eval(['netHist.' to_report{ftr} '(:,t/clock_freq) = state.' to_report{ftr} ';']);
                    
                end

            end
        end

        
        % imagesc(state.curr_noise)
    end


    % display progress
    if ~mod(t,1/delta*10)
        
        if find(isnan(gather(state.V_soma)))
            results.failed = [find(isnan(gather(state.V_soma)))];
            if debugging
                keyboard
            end
            if interrupt_when_fail

                disp('[%===========================================================%]');
                disp('    the simulation returned NaNs : in most cases means that');
                disp('    means that dt is too large. Try reducing it, or check the ');
                disp('    activation variables for inconsistencies.');
                disp('[%===========================================================%]');

                break
            end
        end

        clc

        disp('[================================================]')
        disp(['Running simulation...' wgpu])
        disp(['cell_function: ' cell_function])
        disp(['Complete:' num2str(100*t/simSteps) '%'])
        disp(['Time elapsed:' num2str(toc(tstart)) 's'])
        disp(displaytext)
        disp('[================================================]')
    end


end
disp('Done!')

%% save all state variables of last state
allfields = {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', ...
                        'Calcium_l', 'V_dend', 'Calcium_r', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q_s', ...
                        'Ca2Plus', 'I_CaH', ' V_axon', 'Sodium_h_a' ,  ...
                        'Potassium_x_a', 'I_cx36', 'curr_noise', 'vclamp', ...
                        ' g_gaba_soma', 'g_ampa_dend', ' g_ampa_soma', ' g_gaba_dend', 'Ca2_soma','current', 'backgroundnoise'};

for fts = fields(state)'
    try
        eval( ['tempState.' fts{1} '= gather(state.' fts{1} ');'] ); 
    catch
        eval( ['tempState.' fts{1} '= state.' fts{1} ';'] );
    end

end


% tempState =  [V_soma, Sodium_h, Potassium_n, Potassium_x_s, Calcium_k, Calcium_l, V_dend, Calcium_r, Potassium_s,...
%    Hcurrent_q, Ca2Plus, V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a , g_GABA_soma,g_GABA_dend, g_AMPA, Ca2_soma];


% [=================================================================]
% 
%                      __        __
%   ____  __  ______  / /___  __/ /_
%  / __ \/ / / / __ \/ __/ / / / __/
% / /_/ / /_/ / /_/ / /_/ /_/ / /_
% \____/\__,_/ .___/\__/\__,_/\__/
%           /_/
% 
% [=================================================================]

results.networkHistory = netHist;

results.networksize = networksize;
results.networkParameters.connectivityMatrix = sparse(gather(W));
% results.networkParameters.connectivity = Wstuct;
results.cellParameters = cell_parameters;
results.cellFunction = cell_function;

results.simulationParameters = p.Results;

results.lastState = tempState;
results.simulationParameters.initialState = tempState;
results.simulationParameters.time = time;

results.duration = t*delta;
results.timeElapsed = toc(tstart);

results.perturbation = pert;
results.perturbation.noise = ou_noise;




