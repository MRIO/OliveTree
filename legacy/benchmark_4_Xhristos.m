% benchmark_4_Xhristos.m



dt = .05; fs = 1/dt;
rows = 20; % 50*15
columns = 20;
gap_c = 0.01; % mS/cm^2


alpha  = [32.5 32.5]; %edges for alpha parameter
beta = [12 12];       %edges for beta parameter
res = 1;
scale = 1.35; % 1.35

GCAL = BetaDistributions('alpha', alpha, 'beta', beta, 'rows', rows, 'columns', columns,'resolution', res, 'scale', scale); 
g_CaL = GCAL.sampleDraws{1}';

simtime = 10; %ms 
% noise = randn(rows*columns,fs*simtime)*noise_level;
I_app = 0;

disp('==============================================')
disp('computing transients with decoupled neurons...')

%        ______  ____  _   _________ _____    ____  ____  _______________   _____    __ 
%       / / __ \/ __ \/ | / /_  __( ) ___/   / __ \/ __ \/  _/ ____/  _/ | / /   |  / / 
%  __  / / / / / /_/ /  |/ / / /  |/\__ \   / / / / /_/ // // / __ / //  |/ / /| | / /  
% / /_/ / /_/ / _, _/ /|  / / /    ___/ /  / /_/ / _, _// // /_/ // // /|  / ___ |/ /___
% \____/\____/_/ |_/_/ |_/ /_/    /____/   \____/_/ |_/___/\____/___/_/ |_/_/  |_/_____/
                                                                                      

original_time = tic;


[networkHistory,theSteps] = IOtoplevelCa(dt, simtime, rows, columns, .01, 0, I_APP,0 ,g_CaL)

toc

%    __________  __  __   __  ______   ________________
%   / ____/ __ \/ / / /  /  |/  /   | / ____/  _/ ____/
%  / / __/ /_/ / / / /  / /|_/ / /| |/ / __ / // /     
% / /_/ / ____/ /_/ /  / /  / / ___ / /_/ // // /___   
% \____/_/    \____/  /_/  /_/_/  |_\____/___/\____/   
%                                                      		

gpu_time = tic;
[transients] = IOtoplevelCa_gpu('rows',rows,'columns', columns,'appCurrent',I_app,'time',simtime,'delta',dt,'g_CaL',g_CaL,'g_Gap',0.01 ,'connections','no_torus');

toc

