% pspace tests

% from Schweighofer 1998:
% Using the following parameter values (conductances in mS/
% cm2), soma: g 􏰅 70, g 􏰅 18, g 􏰅 1.0, and g 􏰅 1.5; Na K_dr CA_l h
% dendrite: gCa_h 􏰅 4.0 and gK_Ca 􏰅 35; leak: gls 􏰅 gld 􏰅 0.015 and vl 􏰅 􏰄10 mV; cell morphology:
% gint 􏰅 0.13 and p 􏰅 0.20, data from our standard cell model were consistent with many reported experimental findings from actual IO neurons.
% STEADY-STATE PROPERTIES. The membrane potential with no
% input current (I 􏰅 0 􏰈A/cm2) was 􏰄57 mV, whereas the app
% input resistance derived from the voltage-current curve was 36 Mohm􏰐. When I 􏰅= -􏰄5 􏰈A/cm2, the potential was 􏰄80.3 mV app
% and the input resistance was 14 M􏰐ohm, reflecting the effect of the h current.
% When Iapp 􏰅 􏰃5 􏰈A/cm2, the membrane potential was 􏰄46 mV and the input resistance was only 10 M􏰐 because of the effect of the delayed rectifier current.
% These steady-state values are in agreement with published experimental data (see for instance Table 1 in Llina ́s and Yarom 1981a;
% Fig. 1A in Yarom and Llina ́s 1987; Manor 1995). Note that 1 􏰈A/cm2 corresponds to 0.1 nA for a total cell surface of 10,000 u􏰈m2.


clear

% [================================================]
% 		 globals
% [================================================]

gpu = 1;

% [================================================]
% 		 simulation parameters
% [================================================]

delta = .025;
activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};

to_report = currents;
cell_function = 'vanilla';
% cell_function = 'devel';
steadystate_time = 500;
simtime  = 3000;


% [================================================]
% 		 cell parameters 
% [================================================]


def_neurons = createDefaultNeurons(0 ,'celltypes', 'psweep_gh_gcal');
noneurons = length(def_neurons.g_CaL);
netsize = [noneurons 1 1]

W = zeros(noneurons);
% W = createW()


% [================================================]
% 		 perturbations
% [================================================]

currentstep = 10; %uA/cm^2 -> x .1 nA for a cell with 10000um^2

I_app = zeros(noneurons, simtime*(1/delta));
I_app(:,(110:120)*(1/delta)) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(1500:1600)*(1/delta)) = -currentstep; % nAmpere 20/dt [nA/s.cm^2] 

gnoise = [0 0 0 0];

pert.mask{1}  	  = ones(noneurons,1);
pert.amplitude{1} = 2;
pert.type{1}	  = 'ampa';
pert.duration{1}  = 1;
pert.triggers{1}  = 100;



%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta );
	
end

simcount= 0;

simcount = simcount+1;

% [===========================================================================================================]
   [transients] = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, 'perturbation', pert , ...
   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise , ...
   	'to_report', to_report ,'gpu', gpu ,  ...
   	'cell_function', cell_function ,'delta',delta);
% [===========================================================================================================]
	
simcount = simcount+1;

% [===========================================================================================================]
   [transients] = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, 'perturbation', pert , ...
   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise , ...
   	'to_report', to_report ,'gpu', gpu ,  ...
   	'cell_function', cell_function ,'delta',delta);
% [===========================================================================================================]




R = profile_sim(transients);

figure
imagesc(transients.networkHistory.V_soma,[-150 20]), colorbar
set(gca,'ytick', [1:noneurons],'yticklabel', num2str(def_neurons.Plist),'fontsize',8)

figure
ca = axis;
set(0,'defaultaxescolororder', linspecer(20))
p = plot(transients.networkHistory.V_soma')
legend(num2str(def_neurons.Plist))


figure
imagesc(st_st.networkHistory.V_soma,[-150 20]), colorbar
set(gca,'ytick', [1:noneurons],'yticklabel', num2str(def_neurons.Plist),'fontsize',8)

