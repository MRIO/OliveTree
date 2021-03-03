% Matlab implementation of the model used in:
% De Gruijl et al., 2012, PLoS Computational Biology
%
% [networkHistory,theSteps] = ...
% IOtoplevelCa(delta,time,rows,columns,conductance,randInit,appCurrent,appVoltage,g_CaL)
%
% OUTPUT
% networkHistory = the structure containing the state of every cell for
% every simulated point in time
% 
% theSteps       = the number of time increments that every cell was 
% simulated prior to the actual simulation if random initialization was 
% selected (if not, theSteps == [])
%
% INPUT
% delta       = the time increment (ms)
% time        = the time to be simulated (ms)
% rows        = the number of rows of cells in the grid
% columns     = the number of columns of cells in the grid
% conductance = maximum conductance value for the gap junctions
% randInit    = boolean to toggle initialization with random cell states
% (heterogeneous) or identical cell states (homogeneous)
% appCurrent  = rows by columns by (time/delta) matrix containing values
% for current applied to the dendritic compartment (0 for no current)
% appVoltage  = rows by columns by (time/delta) matrix containing set 
% values for the somatic membrane potential (NaNs for no applied voltage)
% g_CaL       = rows by columns matrix containing for every cell the
% maximum low-threshold calcium conductance (default: 0.7 mS/cm^2)

function [results] = IOtoplevelCa_gpu(varargin)


% [hist,steps] = IOtoplevelCa(0.05,1000,3,3,0.04,1,zeros(3,3,1000/0.05),nan(3,3,1000/0.05),0.7*ones(3,3));


p = inputParser;

p.addParamValue('delta',         .05) % dt
p.addParamValue('time',          1000) % ms
p.addParamValue('rows',          10)
p.addParamValue('columns',       10)
p.addParamValue('g_Gap',   .04) %gap g_Gap
p.addParamValue('randInit',      false)
p.addParamValue('appCurrent',    0)
p.addParamValue('appVoltage',    0)
p.addParamValue('g_CaL',         []) % 0.45 (oscillating) - 1 (spiking)
p.addParamValue('tempState', [])
p.addParamValue('gpu', true) % to be implemented
p.addParamValue('W',[])
p.addParamValue('connections',[]) % 'nearest4'...'yosi'

    p.parse(varargin{:});

delta     = p.Results.delta;
time      = p.Results.time;
rows      = p.Results.rows;
columns   = p.Results.columns;
randInit  = p.Results.randInit;
appCurrent = p.Results.appCurrent;
appVoltage = p.Results.appVoltage;
g_CaL     = p.Results.g_CaL;
g_Gap     = p.Results.g_Gap;
tempState = p.Results.tempState;
W         = p.Results.W;
connections= p.Results.connections;
gpu         = p.Results.gpu;

noNeurons = rows*columns;

if isempty(W)
    if isempty(connections)
        connections = 'nearest8';
        W = createW(connections,rows,columns);
    else
        W = createW(connections,rows,columns);
    end
end


if isempty(g_CaL)
    g_CaL = 0.7*ones(rows*columns,1);
end


tic

disp('Initializing...')

simSteps       = ceil(time/delta);

if ~isempty(tempState)

    for i = 1:23
        eval(['S' num2str(i) '= tempState(:,' num2str(i) ');' ])
        eval(['S' num2str(i) '= gpuArray(S' num2str(i) ');' ])
    end

    S24 =  delta   ;        
    S25 =  0       ;        % I_c   = S25;
    S26 =  0       ;        % I_app = S26;
    S27 =  0       ;        % V_app = S27;
    S28 =  0       ;        % g_CaL = S28;


else

    %% Initial somatic parameters
    tempState.V_soma = -60;
    tempState.PotassiumSlow_soma = 0;
    tempState.PotassiumFast_soma = 0;

    tempState.Calcium_soma       = 0;
    tempState.Leak_soma_dend     = 0;

    % Sodium
    tempState.Sodium_m = 1.0127807;
    tempState.Sodium_h = 0.3596066;

    % Potassium
    tempState.Potassium_n   = 0.2369847;
    tempState.Potassium_x_s = 0.1;

    % Low-threshold calcium
    tempState.Calcium_k = 0.7423159;
    tempState.Calcium_l = 0.0321349;

    %% Initial dendritic parameters
    tempState.V_dend = -60;

    % High-threshold calcium
    tempState.Calcium_r = 0.0112788;

    % Calcium-dependent potassium
    tempState.Potassium_s = 0.0049291;

    % H current
    tempState.Hcurrent_q = 0.0337836;

    % Calcium concentration
    tempState.Ca2Plus = 3.7152;

    % High-threshold calcium current
    tempState.I_CaH   = 0.5;
    tempState.I_K_Ca  = 0;
    tempState.I_cx36  = 0;

    %% Initial axonal parameters
    tempState.V_axon = -60;

    % Sodium
    tempState.Sodium_m_a    = 0.003596066;
    tempState.Sodium_h_a    = 0.9;

    % Potassium
    tempState.Potassium_x_a  = 0.2369847;

    S1=   tempState.V_soma; %V_soma;    
    S2=   tempState.PotassiumSlow_soma ; %I_K_s;     
    S3=   tempState.PotassiumFast_soma ; %I_Kdr_s;   
    S4=   tempState.Calcium_soma       ; %I_CaL;     
    S5=   tempState.Leak_soma_dend     ; %I_ds;      
    S6=   tempState.Sodium_m           ; %m;         
    S7=   tempState.Sodium_h           ; %h;         
    S8=   tempState.Potassium_n        ; %n;         
    S9=   tempState.Potassium_x_s      ; %x_s;       
    S10 =  tempState.Calcium_k          ; %k;        
    S11 =  tempState.Calcium_l          ; %l;        
    S12 =  tempState.V_dend             ; %V_dend;   
    S13 =  tempState.Calcium_r          ; %r;        
    S14 =  tempState.Potassium_s        ; %s;        
    S15 =  tempState.Hcurrent_q         ; %q;        
    S16 =  tempState.Ca2Plus            ; %Ca2Plus;  
    S17 =  tempState.I_CaH              ; %I_CaH;    
    S18 =  tempState.I_K_Ca             ; %I_K_Ca;   
    S19 =  tempState.I_cx36             ; %I_c;      
    S20 =  tempState.V_axon             ; %V_axon;   
    S21 =  tempState.Sodium_m_a         ; %m_a;      
    S22 =  tempState.Sodium_h_a         ; %h_a;      
    S23 =  tempState.Potassium_x_a      ; %x_a;      


    S24 =  delta   ;        
    S25 =  0       ;        % I_c   = S25;
    S26 =  0       ;        % I_app = S26;
    S27 =  0       ;        % V_app = S27;
    S28 =  0       ;        % g_CaL = S28;

    for i = 1:28
        eval(['S' num2str(i) '= S' num2str(i) '*gpuArray(ones(noNeurons,1));' ])
        eval(['S' num2str(i) '=gpuArray(S' num2str(i) ');' ])
    end


end


netHist.V_soma = zeros(noNeurons, time);
netHist.V_dend = zeros(noNeurons, time);
netHist.V_axon = zeros(noNeurons, time);

theSteps = zeros(noNeurons,1);


clock_freq = 1/delta;

V_dend = gpuArray(zeros(noNeurons,1));
I_c = gpuArray(zeros(noNeurons,1));
V_app=NaN(noNeurons,1);
W = gpuArray(W);

tstart  = tic;
for t = 2:simSteps
    
    V_dend = S12;
    
    V = V_dend - W*V_dend;

    f = 0.8 .* exp(-1.*V.*V/100) + 0.2;     % SCHWEIGHOFER 2004 VERSION
    I_c = I_c + (g_Gap .* f .* V);

    S25= I_c;
    S26= appCurrent(:,t);
    S27= V_app;
    S28= g_CaL;


    % IOcellFunctionCa_gpu(delta,I_c,I_app,V_app, State, g_CaL)
    
    [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23] =...
         arrayfun(@IOcellFunctionCa_gpu, ...
            S1(:),S2(:),S3(:),S4(:),S5(:),S6(:),S7(:),S8(:),S9(:),S10(:), ...
            S11(:),S12(:),S13(:),S1 4(:),S15(:),S16(:),S17(:),S18(:),S19(:),S20(:), ...
            S21(:),S22(:),S23(:),S24(:),S25(:),S26(:),S27(:),S28(:));

    if ~mod(t,clock_freq)
        netHist.V_soma(:,t/clock_freq) = gather(S1 );
        netHist.V_dend(:,t/clock_freq) = gather(S12);
        netHist.V_axon(:,t/clock_freq) = gather(S20);
    end

    if ~mod(t,1/delta*10)
        clc
        disp('Running simulation...')
        disp(['Complete:' num2str(100*t/simSteps) '%'])
        disp(['Time elapsed:' num2str(toc(tstart)) 'ms'])
    end


end
disp('Done!')
toc


lastState = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23];

results.networkHistory = netHist;
results.rows = rows;
results.columns = columns;
results.lastState = lastState;
results.duration = t*delta;
results.dt = delta;
results.g_CaL = g_CaL;
results.g_Gap = g_Gap;
results.time = time;
results.timeElapsed = toc(tstart);
results.initialState = tempState;
results.adjacencyMatrix = W;




function W = createW(connections, rows, columns )


switch connections
        case 'all to all'
            W = ones(noNeurons)-eye(noNeurons);
            W = W/(noNeurons-1);
        case 'nearest4'
            E = eye(noNeurons);
            W = circshift(E,[0 -1]) + circshift(E,[0 1]) + circshift(E,[0 -rows]) + circshift(E,[0 rows]);
            W(1,:) = W(1,:)*0; W(end,:) = W(end,:)*0;
            W = bsxfun(@rdivide, W, sum(W,2));
            
        case 'nearest8'
            E = eye(noNeurons);
            W = circshift(E,[0 -1]) + circshift(E,[0 1]) + circshift(E,[0 -rows]) + circshift(E,[0 rows]);
            W = bsxfun(@rdivide, W, sum(W,2));

        case 'no_torus'
            CONNECTED = 8;
            %# which distance function
            if CONNECTED == 4,     distFunc = 'cityblock';
            elseif CONNECTED == 8, distFunc = 'chebychev'; 
            end

            %# compute adjacency matrix
            [X Y] = meshgrid(1:rows,1:columns);
            X = X(:); Y = Y(:);
            W = squareform( pdist([X Y], distFunc) == 1 );
            W = bsxfun(@rdivide, W, sum(W,2));
        case 'random'
            CONNECTED = 8;
            %# which distance function
            if CONNECTED == 4,     distFunc = 'cityblock';
            elseif CONNECTED == 8, distFunc = 'chebychev'; end

            %# compute adjacency matrix
            [X Y] = meshgrid(1:rows,1:columns);
            X = X(:); Y = Y(:);
            W = squareform( pdist([X Y], distFunc) == 1 );


            W = bsxfun(@rdivide, W, sum(W,2));

        case 'yosi'
            CONNECTED = 8;
            %# which distance function
            if CONNECTED == 4,     distFunc = 'cityblock';
            elseif CONNECTED == 8, distFunc = 'chebychev'; 
            end

            %# compute adjacency matrix
            
            [X Y] = meshgrid(1:rows,1:columns);[]
            X = X(:); Y = Y(:);
            W = squareform( pdist([X Y], distFunc) == 1 );

            mx = ceil(rows/2);
            my = ceil(columns/2);

            % break a division
            boundary_left = sub2ind([rows*columns,rows*columns] ,[mx:rows:rows*columns] , [mx-1:rows:rows*columns]);
            boundary_right = sub2ind([rows*columns,rows*columns],[mx:rows:rows*columns] , [mx+1:rows:rows*columns]);
            
            W(boundary_right) = .001;
            W(boundary_left) = .001;

            W(boundary_right) = .000;
            W(boundary_left) = .000;


            % W(mx:rows:rows*columns,:) = 0;
            W(:,mx:rows:rows*columns) = 0;
            
            

            % normalize weights
            W = bsxfun(@rdivide, W, sum(W,2));

    case 'yosi'
        CONNECTED = 8;
            %# which distance function
            if CONNECTED == 4,     distFunc = 'cityblock';
            elseif CONNECTED == 8, distFunc = 'chebychev'; 
            elseif CONNECTED == 100, distFunc = 'minkowski'; 
            end

            %# compute adjacency matrix
            
            [X Y] = meshgrid(1:rows,1:columns);
            X = X(:); Y = Y(:);
            W = squareform( pdist([X Y], distFunc) == 1 );

            D = squareform(pdist([X Y], distFunc));


            mx = ceil(rows/2);
            my = ceil(columns/2);

            % break a division
            boundary_left = sub2ind([rows*columns,rows*columns] ,[mx:rows:rows*columns] , [mx-1:rows:rows*columns]);
            boundary_right = sub2ind([rows*columns,rows*columns],[mx:rows:rows*columns] , [mx+1:rows:rows*columns]);
            
            W(boundary_right) = .001;
            W(boundary_left) = .001;

            W(boundary_right) = .000;
            W(boundary_left) = .000;


            % W(mx:rows:rows*columns,:) = 0;
            W(:,mx:rows:rows*columns) = 0;

            
            

            % normalize weights
            W = bsxfun(@rdivide, W, sum(W,2));



end

[a b] = ind2sub([rows, columns], [1:rows*columns]);
gplot(double(W),[a' b']);
title(connections)
            



