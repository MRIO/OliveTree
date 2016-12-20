% noise_corr_gap_comparisons.m


% conntype = 'cluster'; HPCGPU_netpspace_noise_corr_gap;
% conntype = 'iso'    ; HPCGPU_netpspace_noise_corr_gap;

% [=================================================================]
%  script parameters
% [=================================================================]


compsim = 1; gpu = 1;
	savehist = 1;
summarize = 0;
	calcfreq = 0;
parametertable = 0;

rndseed = 0;
rng(rndseed,'twister')

fname = ['pspace_conn_rad__' date ];

% [=================================================================]
%  parameter space grid
% [=================================================================]

gaps = [0.05];
noiseamps = [1.2]; %uA/cm^2 -> x .1 nA for a cell with 10000um^2


if ~exist('noisecorr');noisecorr = .1;end
if ~exist('simtype')  ; simtype  = '1Hz' ;end
if ~exist('conntype') ; conntype = 'iso' ;end
if ~exist('numruns')  ; numruns  = 2     ;end
if ~exist('sametoall'); sametoall = 0.1  ;end
if ~exist('tau')	  ; tau = 5  		 ;end
if ~exist('noiseamp') ; noiseamp = 1.5	 ;end

displaytext = [simtype '_' conntype '_' num2str(numruns) '_' num2str(sametoall)];


% [=================================================================]
%  simulation parameters
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
steadystate_time = 300;
simtime  = 50;

delta = .02; 
Fs = 1000; % sampling rate =! delta

% [================================================]
%  net parameters
% [================================================]

netsize = [4 10 10];
	noneurons = prod(netsize);

rds = [3 4 5 6];
meannoconns = [6 8 10 12];

plotconn  = 1;
normleak  = 1;
randomize = 0; %!check
scaling   = 1;
plotthis  = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;

c = 0;
for rd = rds;
	for meannoconn = meannoconns
	c = c+1'
		W{c} = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
		W{c}.radius = rd;
		W{c}.meannoconn = meannoconn;
	end
end

% [=================================================================]
%  Neurons
% [=================================================================]

rng(rndseed,'twister')
def_neurons = createDefaultNeurons(noneurons,'celltype', 'randomized');
rng(rndseed,'twister')
gapless_neurons = createDefaultNeurons(noneurons,'celltype', 'randomized','nogapcompensation', 1);
Plist = def_neurons.Plist;


disp('[=================================================================]')
disp( ['using:' conntype])
disp( ['netsize:' num2str(netsize)])
disp( ['gap conditions:' num2str(gaps   )])
disp( ['noisecorr:' num2str(noisecorr)])
disp( ['noiseamps:' num2str(noiseamps)])

disp('[=================================================================]')


% [=================================================================]
%  readouts
% [=================================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};

to_report = vsoma;


% [=================================================================]
%  OU noise parameters
% [===========================================================q======]

tau = 5;
theta = 1/tau;








%===========================================================%
%     ____  __  ___   __   _____ ______  ________
%    / __ \/ / / / | / /  / ___//  _/  |/  / ___/
%   / /_/ / / / /  |/ /   \__ \ / // /|_/ /\__ \ 
%  / _, _/ /_/ / /|  /   ___/ // // /  / /___/ / 
% /_/ |_|\____/_/ |_/   /____/___/_/  /_//____/  
% 
%===========================================================%            

if compsim

	%%================================================]
	% 		 compute transients/steadystate
	%=================================================]
	if ~exist('st_st','var')
		disp('calculating transients')
		st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
		st_st.Plist = def_neurons.Plist;
	end


	I_app = [];



% [===========================================================================================================]
   simcount= 0;
   for simcount = 1:length(W)
   	
   	 	displaytext = [ num2str(W{simcount}.rd) '_' '_' num2str(ncorr)];


   	 	rng(0,'twister')
   	 	gnoise = [theta namp 0 5];
   	 	simcount = simcount+1;

   	 	if gs > 0.01
		   neurons = def_neurons;
		else
		   neurons = gapless_neurons;
		end

		transients{simcount} = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons , ...
		   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W{simcount}.W*gs ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu ,  ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', ncorr,'saveappliednoise',0,...
		   	'displaytext', displaytext);

		   transients{simcount}.scriptParameters.noiseAmplitude   = namp;
		   transients{simcount}.scriptParameters.noiseCorrelation = ncorr;
		   transients{simcount}.scriptParameters.gapAmplitudes    = gs;
		   transients{simcount}.spikes = spikedetect(transients{simcount});
		   transients{simcount}.Plist = Plist;
		   transients{simcount}.W = W{simcount};

		   if not(savehist)
		   		transients{simcount}.networkHistory = [];
		   end

			eval(['save ' fname ' -v7.3'])
				
end

disp(['computed ' num2str(simcount) ' simulations.'])
   

% [===========================================================================================================]

% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'), legend(num2str(Plist)),title('V vs Calcium\_r')

%=============================Make Parameter Table ==============================%

