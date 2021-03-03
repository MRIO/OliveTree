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

function [networkHistory,theSteps] = IOtoplevelCa(delta, time, rows, columns, conductance, randInit, appCurrent,appVoltage,g_CaL)

tic

disp('Initializing...')

simSteps       = ceil(time/delta);

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


%% fill starting state slice
         
networkHistory(1:rows, 1:columns,1:simSteps) = repmat(tempState,[rows,columns,simSteps]);
I_app = 0;        % applied current
V_app = NaN;      % applied voltage

theSteps = zeros(rows,columns);

if randInit == 1
    disp('Running prelude (randInit == 1) ...')
    for n = 1:rows    
        for m = 1:columns
            
            % how many steps of init per cell
            theSteps(n,m) = randi(round(100/delta),1);
            
            for t = 1:theSteps(n,m)
        
                networkHistory(n, m, 1) = IOcellFunctionCa(delta, 0, I_app, V_app,networkHistory(n,m,1),g_CaL(n,m));
        
            end            
        end    
    end
else
    disp('No prelude (randInit == 0) ...')
    theSteps = [];
end

disp('Running simulation...')

for n = 2:simSteps
    
    for i = 1:rows        
        for j = 1:columns
            
            I_c = 0;
            
            for k = -1:1                
                for m = -1:1                    
                    if (((i + k > 0) && (j + m > 0)) && ((i + k <= rows) && (j + m <= columns)))
                        V = (networkHistory(i, j, n - 1).V_dend - networkHistory(i + k, j + m, n - 1).V_dend);
                        
                        f = 0.8 .* exp(-1.*V.*V/100) + 0.2;     % SCHWEIGHOFER 2004 VERSION
                        I_c = I_c + (conductance .* f .* V);
                        
                    end                    
                end                
            end
            
            I_app = appCurrent(i,j,n);
            V_app = appVoltage(i,j,n);
            networkHistory(i, j, n) = IOcellFunctionCa(delta, I_c, I_app, V_app, networkHistory(i,j,n-1),g_CaL(i,j));
            
        end        
    end    
end
disp('Done!')
toc