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

function [results] = IOtoplevelCa_M(varargin)


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
conductance = p.Results.g_Gap;
randInit  = p.Results.randInit;
appCurrent = p.Results.appCurrent;
appVoltage = p.Results.appVoltage;
g_CaL     = p.Results.g_CaL;
tempState = p.Results.tempState;
W         = p.Results.W;

simSteps        = ceil(time/delta);
noNeurons       = rows*columns;


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

disp('Initializing...')


    if ~isempty(tempState)
        networkHistory(1:noNeurons,1:simSteps) = repmat(tempState,[1,simSteps]);
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

    networkHistory(1:noNeurons,1:simSteps) = repmat(tempState,[noNeurons,simSteps]);
end


%% fill starting state slice

%todo: check for network size compatibility



theSteps = zeros(rows*columns,1);

if randInit == true
    disp('Running prelude (randInit == 1) ...')
    for nrn = 1:length(W)
        
        
        theSteps(nrn) = randi(round(100/delta),1);
        
        if appCurrent == 0
            I_app = 0;
        else
            I_app = appCurrent(nrn,1);
        end

        if appVoltage == 0
            V_app = NaN;
        else
            V_app = appVoltage(nrn,1);
        end

            

        for t = 1:theSteps(nrn)
    
            networkHistory(nrn, 1) = IOcellFunctionCa(delta, 0, I_app, V_app, networkHistory(nrn,1), g_CaL(nrn));
    
        end

   
    end
else
    disp('No prelude (randInit == 0) ...')
    theSteps = [];
end





% % cellState = IOcellFunctionCa(delta,I_c,I_app,V_app,previousState, g_CaL)
% % rd = @(i,j,k) networkHistory(i,j,k).V_dend;
% % W  = rand(10,10);
% keyboard



V_dend = zeros(noNeurons,1);
I_c = zeros(noNeurons,1);


tstart  = tic;
for t = 2:simSteps
    
    % this is unfortunatelly fairly slow
    % [c_i c_j] = find(ones(rows,columns));
    % tic
    % V_dend = reshape(arrayfun(@(i,j) networkHistory(i,j,n).V_dend, c_i,c_j),5,5);
    % toc

    % we do this, instead
    tic
            for nrn = 1:noNeurons
                    V_dend(nrn,1) = networkHistory(nrn,t-1).V_dend;
            end

                    V = V_dend - W*V_dend;
            

                    f = 0.8 .* exp(-1.*V.*V/100) + 0.2;     % SCHWEIGHOFER 2004 VERSION
                    I_c = I_c + (conductance .* f .* V);
                    
                    


                    if appCurrent == 0
                        I_app = 0;
                    else
                        I_app = appCurrent(nrn,t);
                    end

                    % if appVoltage == 0
                        V_app = NaN;
                    % else
                    %     V_app = appVoltage(nrn,t);
                    % end


                % end
            
            for nrn = 1:noNeurons
                networkHistory(nrn, t) = IOcellFunctionCa(delta, I_c(nrn), I_app, V_app, networkHistory(nrn,t-1), g_CaL(nrn));
            
            end
    
        

            
    
    if ~mod(t,100)

    clc
    disp('Running simulation...')
    disp(['Complete:' num2str(t/simSteps) '%'])
    disp(['Time elapsed:' num2str(toc(tstart)) 'ms'])
    end

    if t==simSteps
        try
        lastState = networkHistory(:,t)
        catch
            keyboard
        end
    end

end
disp('Done!')
toc





results.networkHistory = networkHistory;
results.rows = rows;
results.columns = columns;
results.lastState = lastState;
results.duration = t*delta;
results.dt = delta;
results.g_CaL = g_CaL;
results.g_Gap = conductance;
results.time = time;
results.timeElapsed = toc(tstart);






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
            





conductance;