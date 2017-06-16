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


% [=================================================================]
%  global parameters
% [=================================================================]
clear
rng(0,'twister')
gpu = 1;

debugging = 0;
% [=================================================================]
%  simulation parameters
% [=================================================================]
delta = .02;

% cell_function = 'devel';
cell_function = 'vanilla'; % 'devel'
nconn = 3;

steadystate_time = 1000;
simtime  = 1000;

% currentstep = 5; %uA/cm^2 -> x .1 nA for a cell with 10000um^2

% [=================================================================]
%  state variables to report
% [=================================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36' 'curr_noise_pert'};
vs = {'V_soma' ,'V_dend','V_axon' };
to_report = currents;


% [=================================================================]
%  parameter grid
% [=================================================================]
% gaps = [];
gaps = [0];

pspace_type = 'pgrid';
% pspace_type = 'randomized';%
 % pspace_type = '2p_sweep';
switch pspace_type


	case 'randomized'

	noneurons = 100;
	netsize = [noneurons 1 1];
	def_neurons = createDefaultNeurons(noneurons,'celltypes','randomized2');
	Plist = def_neurons.Plist;
	Pnames = def_neurons.Pnames;

	case '2p_sweep'
	
	% def_neurons = createDefaultNeurons(1,'celltypes','param_sweep','Pnames',{'g_CaL' ; 'g_int'});
	def_neurons = createDefaultNeurons(1,'celltypes','psweep_gh_gcal');
	
	Plist = def_neurons.Plist;
	Pnames = def_neurons.Pnames;
	noneurons = length(Plist);
	netsize = [noneurons 1 1];
	

	case 'pgrid'

	% 9 Dimensional GRID: parameter ranges
	p1 = [1 2 3]; 		% CalciumL - conductance range
	p2 = [0];      	    % g_h_s
	p3 = [.15]; 	% g_int
	p4 = [0 1];      	% g_h
	p5 = [.15]; % ratio soma dendrite
	p6 = [35 55 85];	% Ca act Potassium: not voltage dependent 
	p7 = [4.5 6.5]; % Ca High threshold
	p8 = [0.013 0.016];    % leak
	p9 = [0.013 0.016];; % leak
	p10 = [0];; % arbitrary

	[p{1} p{2} p{3} p{4} p{5} p{6} p{7} p{8} p{9} p{10}] = ndgrid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);

	Pnames = {'g_CaL' 'g_h_s' 'g_int' 'g_h' 'p1' 'g_K_Ca' 'g_CaH' 'g_ld' 'g_ls' 'gbar_ampa_soma'}';
	Plist = [p{1}(:) p{2}(:) p{3}(:) p{4}(:) p{5}(:) p{6}(:) p{7}(:) p{8}(:) p{9}(:) p{10}(:)]; 

	noneurons = length(p{1}(:));
	netsize = [1 noneurons 1];noneurons = prod(netsize);

	def_neurons = createDefaultNeurons(noneurons);
	def_neurons.g_CaL    = p{1}(:);
	def_neurons.g_h_s    = p{2}(:);
	def_neurons.g_int 	 = p{3}(:);
	def_neurons.g_h 	 = p{4}(:);
	def_neurons.p1  	 = p{5}(:);
	def_neurons.g_K_Ca   = p{6}(:);       
	def_neurons.g_CaH    = p{7}(:);     % High-threshold calcium
	% def_neurons.arbitrary= p{8}(:);
	def_neurons.g_ld = p{8}(:);
	def_neurons.g_ls = p{9}(:);
	def_neurons.gbar_ampa_soma = p{10}(:);


end


	W = createW('all to all', [1 noneurons 1], [], 1, 0, 0, 0, nconn, []);


% [=================================================================]
%  perturbation
% [=================================================================]

ounoise = [0 0 0 0];
% ounoise = [1/20  .2 0 0];
% ounoise = [1/20 .2 0 0];



% apply some current to check the behavior of the cells
I_app = [];
% I_app = ones(noneurons, simtime*1/delta)*2;
% I_app(:,(100*(1/delta):110*(1/delta))) = currentstep; % nAmpere 20/dt [nA/s.cm^2] 
% I_app(:,(500*(1/delta):510*(1/delta))) = -currentstep;  % nAmpere 20/dt [nA/s.cm^2] 


% mask = create_input_mask(...)

% perturbation_onsets{2} = 100;
% perturbation_mask{2} 	= glu_mask; % selection of stimulated neurons
% perturbation_type{2} = 'ampa';
% perturbation_pulse_duration{2} = exc_pulse_dur;

% perturbation_onsets{3} = gaba_onsets_dend;
% perturbation_mask{3} 	= gaba_mask_dend; % selection of stimulated neurons
% perturbation_type{3} = 'gaba_dend';
% perturbation_pulse_duration{3} =  gaba_pulse_dur;

% perturbation_onsets{4} = gaba_onsets_soma;
% perturbation_mask{4} 	= gaba_mask_soma; % selection of stimulated neurons
% perturbation_type{4} = 'gaba_soma';
% perturbation_pulse_duration{4} =  gaba_pulse_dur;					

% pert.triggers{1} = .001; % proportion of neurons receiveing ampa per bin
% pert.mask{1} 	= create_input_mask(netsize,'all'); % selection of stimulated neurons
% pert.type{1} = 'ampa_noise';
% pert.duration{1} =  1;					


pert = [];



%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta,'ou_noise', [0 0 0 0],'debug',debugging);
	st_st.Plist = Plist;
	st_st.Pnames = Pnames;

	% replayResults_3(st_st, 'plotallfields',1)
	% resultstable = profile_sim(st_st);
	drawnow

end



% [===========================================================================================================]
s = 0;
for gap = gaps
	s = s+1;
   [simresults{s}] = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, ...
   	'networksize', netsize,'time',simtime ,'W', W.W*gap ,'ou_noise', ounoise , ...
   	'to_report', to_report ,'gpu', gpu , 'appCurrent', I_app, ...
   	'cell_function', cell_function ,'delta',delta,'sametoall', 0 , 'displaytext', ['sim: ' num2str(s)]);

   simresults{s}.Plist = Plist;
   simresults{s}.gaps = gap;

   
	%  pos = subplus(simresults{s}.networkHistory.I_cx36');
	%  neg = -subplus(-simresults{s}.networkHistory.I_cx36');

	% figure
	% area([pos ; neg])

	% receiving = area(pos,'edgecolor','none'), hold on
	% donating  = area(neg,'edgecolor','none')


end
% [===========================================================================================================]
R = profile_sim(simresults{1},'tslice', [1:1000], 'plotme',1);
sel_fields = {'g_CaL', 'g_int', 'g_h', 'g_ld' ,'g_ls', 'ampl', 'freq_each', 'meanVm'};
NDscatter(R.allneurons(:,sel_fields),1)


prunecells = 0;
if prunecells
	

	% criteria for cell selection
	ampl   = table2array(R.allneurons(:,'ampl'))  <30;
	freq   = table2array(R.allneurons(:,'freq_each'))  >2 | table2array(R.allneurons(:,'freq_each'))<10 ;
	maxV   = table2array(R.allneurons(:,'maxV'))  <-30;
	meanVm = table2array(R.allneurons(:,'meanVm'))<-55;


	sel_cel_idx = ampl & freq & maxV & meanVm;

	waterfall(simresults{1}.networkHistory.V_soma(sel_cel_idx,:)')
	plot(simresults{1}.networkHistory.V_soma(sel_cel_idx,:)')

	sel_fields = {'g_CaL', 'g_K_Ca', 'g_int', 'p1',  'ampl', 'freq_each', 'maxV', 'meanVm'}
	sel_table = R.allneurons(sel_cel_idx,sel_fields);
	% sel_table = R.allneurons(:,sel_fields);
	NDscatter(sel_table, 1)


	save cellset_devel_3 sel_cel_idx simresults resultstable
else
		sel_fields = {'g_CaL', 'g_int', 'p1',  'ampl', 'freq_each', 'meanVm'}
		sel_table = R.allneurons(:,sel_fields);
		% sel_table = R.allneurons(:,sel_fields);
		NDscatter(sel_table, 1)
end

