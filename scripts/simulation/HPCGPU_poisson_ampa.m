% Input correlations and center peak sharpnesss?

% pspace tests for distributions of ampa input over a wee time


% [================================================]
% 		 simulation parameters
% [================================================]
% clear
randomseed = 1;

cell_function = 'vanilla'; % 'devel'
steadystate_time = 1000; %ms
simtime  = 20000; %ms
delta = .025;
gpu = 1;


activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};
gapcur= {'V_soma' 'I_cx36' 'Calcium_r'};

% variables to store
to_report = gapcur;


% [================================================]
% 		 define parameter space to create grid
% [================================================]
% 8 Dimensional GRID: parameter ranges
p1 = [.5:.1:1]; 		% CalciumL - conductance range
p2 = [240];     	    % g_Na_a
p3 = [.12 .15]; 		% g_int
p4 = [.12:.12:.48];     % g_h
p5 = [-38];       		% V_h
p6 = [45 55];			% Ca act Potassium: not voltage dependent 
p7 = [4.5];
p8 = [0.012 0.013];      % leak
[p{1} p{2} p{3} p{4} p{5} p{6} p{7} p{8}] = ndgrid(p1,p2,p3,p4,p5,p6,p7,p8);

Plist = [p{1}(:) p{2}(:) p{3}(:) p{4}(:) p{5}(:) p{6}(:) p{7}(:) p{8}(:)]; 

noneurons = length(p{1}(:));
netsize = [1 noneurons 1];
 

% [================================================]
% 		 create default neurons
% [================================================]
def_neurons = createDefaultNeurons(noneurons);
twins = createDefaultNeurons(noneurons);
def_neurons.g_CaL    = p{1}(:);
def_neurons.g_Na_a   = p{2}(:);
def_neurons.g_int 	 = p{3}(:);
def_neurons.g_h 	 = p{4}(:);
def_neurons.V_h 	 = p{5}(:);
def_neurons.g_K_Ca   = p{6}(:);       
def_neurons.g_CaH    = p{7}(:);     % High-threshold calcium
def_neurons.g_ld = p{8}(:);
def_neurons.g_ls = p{8}(:);
def_neurons.g_la = p{8}(:);

def_neurons.gbar_ampa_soma = ones(noneurons,1)*.15; % linspace(0,.3,noneurons); %
def_neurons.gbar_gaba_soma = ones(noneurons,1)*1;

% compensation for gap junctions
gap_neurons = def_neurons;
gap_neurons.g_CaL = def_neurons.g_CaL + .1;
gap_neurons.g_int = def_neurons.g_int -0.04;
gap_neurons.g_ld = def_neurons.g_ld - 0.003;
gap_neurons.g_K_Ca = def_neurons.g_K_Ca + 10;


% [================================================]
% 		gap connections
% [================================================]
 W = zeros(noneurons);
 % W = createW(noneurons);

gaps = [0 .025 .05];
W = createW('all to all', [1 noneurons 1], [], 1, 0, 0, 0, 10, []);

% [================================================]
% 		 input
% [================================================]

% gnoise = [3 5 0 0];
gnoise = [0 0 0 0];
sametoall = 0;

% [================================================]
%  Distribute Ampa Perturbation over time and masks
% [================================================]

% create overlapping ampa masks

numberofmasks = [10];
masktype = 'all';
synapseprobability = .1;
stim_interval = [25 50 100 250];
onset_of_stim = 400;
n_of_pulses = 100; % make it a function of the number of pulses == energy delivery
last_stim_time = 100; % make it a function of the total stimlus time ?

%%================================================]
% 		 compute transients/steadystate
%=================================================]
if ~exist('st_st','var')
	disp('calculating transients')
	st_st = IOnet('cell_function', cell_function ,'networksize', netsize, 'cell_parameters', def_neurons, 'time', steadystate_time ,'gpu', gpu,'to_report', to_report ,'delta',delta);
	st_st.Plist = Plist;
end


% [===========================================================================================================]
 ns = 0;
 for g = gaps
 	
 	rng(randomseed,'twister') % random seed

 	for si = stim_interval

	 	ns = ns +1;
	 	thissimtime = simtime;
 		
 		for nm = 1:numberofmasks

			pert.mask  	  {nm} = create_input_mask(netsize, masktype, 'synapseprobability', synapseprobability);
			pert.amplitude{nm} = 1;
			pert.triggers {nm} = onset_of_stim + round(si*rand) + cumsum(poissrnd(si,n_of_pulses,1)) ;
			pert.duration {nm} = 1;
			pert.type	  {nm} = 'ampa';
			
		end 

		thissimtime = max(max(cell2mat(pert.triggers))) + 1000
		pert.all_pulses = cell2mat(pert.triggers);
		pert.inp_per_cell = sum(cell2mat(pert.mask));

		simtime = max(pert.all_pulses(:));

		if g > 0
			neurons = gap_neurons ;
		else 
			neurons = def_neurons ;
		end	



	   [simresults{ns}] = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons, ...
	   	'networksize', netsize,'time',thissimtime ,'W', W.W*g ,'ou_noise', gnoise , ...
	   	'to_report', to_report ,'gpu', gpu , 'perturbation', pert, ...
	   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall);

	   simresults{ns}.Plist = Plist;
	end
end
% [===========================================================================================================]

save simresults simresults

% spks{1} = spikedetect(simresults{1}, 0,0);
% spks{2} = spikedetect(simresults{2}, 0,0);


% figure
% imagesc(st_st.networkHistory.V_soma,[-80 -20]), colorbar
% set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
% legend(num2str(Plist))


% figure
% ca = axis;
% set(0,'defaultaxescolororder', linspecer(length(Plist)))
% pl = plot([1:steadystate_time],   st_st.networkHistory.V_soma');
% legend(num2str(Plist))



% figure
% imagesc(simresults{1}.networkHistory.V_soma,[-80 -20]), colorbar
% set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
% title([num2str(spks{1}.popfrequency) ' Hz'])


% figure
% imagesc(simresults{2}.networkHistory.V_soma,[-80 -20]), colorbar
% set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
% title([num2str(spks{2}.popfrequency) ' Hz'])


% figure
% ca = axis;
% set(0,'defaultaxescolororder', linspecer(length(Plist)));
% pl = plot([1:simtime],   simresults{2}.networkHistory.V_soma');
% legend(num2str(Plist))



% figure
% ca = axis;
% set(0,'defaultaxescolororder', linspecer(length(Plist)));
% pl = plot([1:simtime],   simresults{2}.networkHistory.V_soma');
% legend(num2str(Plist))


% pl = plot([1:simtime],   simresults{1}.networkHistory.V_soma');
% R = replayResults(simresults{1},[],0)


% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'),legend(num2str(Plist)),title('V vs Calcium\_r')


 pos = subplus(simresults{2}.networkHistory.I_cx36');
 neg = -subplus(-simresults{2}.networkHistory.I_cx36');

area([pos ; neg])

	receiving = area(pos,'edgecolor','none'); hold on
	donating  = area(neg,'edgecolor','none');






% [=================================================================]
%  
% [=================================================================]
% si = 200;

% for nm = 1:numberofmasks

% 	pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', .25);
% 	pert.amplitude{nm} = 1;
% 	pert.triggers {nm} = onset_of_stim + round(si*rand) + cumsum(poissrnd(si,n_of_pulses,1)) ;
% 	pert.duration {nm} = 1;
% 	pert.type	  {nm} = 'ampa';
	
	
% end 

% P = cell2mat(pert.triggers);



