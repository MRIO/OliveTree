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
% These steady-state values are in agreement with published experimental data (see for instance Table 1 in Llinas and Yarom 1981a;
% Fig. 1A in Yarom and Llina ́s 1987; Manor 1995). Note that 1 􏰈A/cm2 corresponds to 0.1 nA for a total cell surface of 10,000 u􏰈m2.


clear
rng(0,'twister')

gpu = 1;

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36'};

to_report = gapcur;
cell_function = 'vanilla'; % 'devel'
steadystate_time = 3;
simtime  = 1000;

delta = .025;
currentstep = 9; %uA/cm^2 -> x .1 nA for a cell with 10000um^2

% 8 Dimensional GRID: parameter ranges
p1 = [1]; 		% CalciumL - conductance range
p2 = [.25];     	    % C_m
p3 = [.13]; 		% g_int
p4 = [.12];      	% g_h
p5 = [-38];       	% V_h
p6 = [45];			% Ca act Potassium: not voltage dependent 
p7 = [4.5];
p8 = [0.013];    % leak
[p{1} p{2} p{3} p{4} p{5} p{6} p{7} p{8}] = ndgrid(p1,p2,p3,p4,p5,p6,p7,p8);

Plist = [p{1}(:) p{2}(:) p{3}(:) p{4}(:) p{5}(:) p{6}(:) p{7}(:) p{8}(:)]; 

noneurons = length(p{1}(:));
netsize = [1 noneurons 1];
noneurons = prod(netsize);

gnoise = [3 1 0 0];
% gnoise = [0 0 0 0];

def_neurons = createDefaultNeurons(noneurons);
twins = createDefaultNeurons(noneurons);
def_neurons.g_CaL    = p{1}(:);
def_neurons.C_m    = p{2}(:);
def_neurons.g_int 	 = p{3}(:);
def_neurons.g_h 	 = p{4}(:);
def_neurons.V_h 	 = p{5}(:);
def_neurons.g_K_Ca   = p{6}(:);       
def_neurons.g_CaH    = p{7}(:);     % High-threshold calcium
% def_neurons.arbitrary= p{8}(:);
def_neurons.g_ld = p{8}(:);
def_neurons.g_ls = p{8}(:);
def_neurons.g_la = p{8}(:);


 W = zeros(noneurons);
 % W = createW(noneurons);

gap = 0.05;
W = createW('all to all', [1 noneurons 1], [], gap, 0, 0, 0, 10, []);

%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
	st_st.Plist = Plist;
end

simcount= 0;

simcount = simcount+1;

% apply some current to check the behavior of the cells
I_app = [];
% I_app = zeros(noneurons, simtime*(1/delta));
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 



% [===========================================================================================================]
   [simresults] = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, ...
   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise , ...
   	'to_report', to_report ,'gpu', gpu ,  ...
   	'cell_function', cell_function ,'delta',delta,'sametoall', .0);

   simresults.Plist = Plist;
% [===========================================================================================================]
	
spks = spikedetect(simresults, 0,0);


figure
imagesc(st_st.networkHistory.V_soma,[-80 -20]), colorbar
set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
legend(num2str(Plist))


figure
ca = axis;
set(0,'defaultaxescolororder', linspecer(length(Plist)))
p = plot([1:steadystate_time],   st_st.networkHistory.V_soma');
legend(num2str(Plist))


figure
imagesc(simresults.networkHistory.V_soma,[-80 -20]), colorbar
set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
title([num2str(spks.popfrequency) ' Hz'])

figure
ca = axis;
set(0,'defaultaxescolororder', linspecer(length(Plist)));
p = plot([1:simtime],   simresults.networkHistory.V_soma');
legend(num2str(Plist))


p = plot([1:simtime],   simresults.networkHistory.V_soma');
R = replayResults(simresults,[],0)


% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'),legend(num2str(Plist)),title('V vs Calcium\_r')


 pos = subplus(simresults.networkHistory.I_cx36');
 neg = -subplus(-simresults.networkHistory.I_cx36');

area([pos ; neg])

	receiving = area(pos,'edgecolor','none'), hold on
	donating  = area(neg,'edgecolor','none')


