% createDefaultNeurons.m
function cell_parameters = createDefaultNeurons(varargin)
% createDefaultNeurons(noneurons, celltypes, gapcompensation)

%TODO: Check for correlations in the random draws - if there are, shuffle/redraw

	ip = inputParser;
	ip.addOptional('noneurons',0)
	
	ip.addParameter('celltypes', 'none') 
	ip.addParameter('Pnames', {'g_CaL' ;'g_int'; 'g_h'; 'g_K_Ca'; 'g_ld' ;'p1'} )
	ip.addParameter('gapcompensation', 0) 
	ip.addParameter('nogapcompensation', 0);
	ip.addParameter('shuffle', 0) 
	ip.addParameter('addrand',0)
	ip.addParameter('rng', rng(0))
	
	ip.parse(varargin{:});

	noneurons = ip.Results.noneurons;
	celltypes = ip.Results.celltypes;
	gapcompensation = ip.Results.gapcompensation;
	shuffle = ip.Results.shuffle;
	nogapcompensation = ip.Results.nogapcompensation;
	addrand = ip.Results.addrand;
	Pnames = ip.Results.Pnames;

	rng(ip.Results.rng)

	
cell_parameters = defneurons(noneurons);

switch celltypes
	case 'none'
		% do nothing

	case 'randomized'

		cell_parameters = defneurons(noneurons);
		
		cell_parameters.g_CaL    			= cell_parameters.g_CaL   +  rand(noneurons,1)*(-.6);
		cell_parameters.g_int 	 			= cell_parameters.g_int   +  rand(noneurons,1)*(-.02);
		cell_parameters.g_h 	 			= cell_parameters.g_h 	  +  rand(noneurons,1)*(1);
		cell_parameters.g_K_Ca   			= cell_parameters.g_K_Ca  +  rand(noneurons,1)*10;       
		cell_parameters.g_ld     			= cell_parameters.g_ld    +  rand(noneurons,1)*(-0.003);
		cell_parameters.g_la     			= cell_parameters.g_la    +  rand(noneurons,1)*(-0.003);
		cell_parameters.g_ls     			= cell_parameters.g_ls    +  rand(noneurons,1)*(-0.003);
		% cell_parameters.gbar_ampa_soma      = .5   -  rand(noneurons,1)*(.15);


	case 'randomized2'

		cell_parameters = defneurons(noneurons);
		
		cell_parameters.g_CaL    = cell_parameters.g_CaL   - .3  + rand(noneurons,1)*.6;
		cell_parameters.g_int 	 = cell_parameters.g_int   + .07 + rand(noneurons,1)*.4;
		cell_parameters.g_h 	 = cell_parameters.g_h 	   +  rand(noneurons,1)*(1);
		cell_parameters.g_K_Ca   = cell_parameters.g_K_Ca  +  rand(noneurons,1)*10;     
		cell_parameters.g_ld     = cell_parameters.g_ld    +  rand(noneurons,1)*(-0.003);
		cell_parameters.g_la     = cell_parameters.g_la    +  rand(noneurons,1)*(-0.003);
		cell_parameters.g_ls     = cell_parameters.g_ls    +  rand(noneurons,1)*(-0.003);
		cell_parameters.p1       = cell_parameters.p1      - .15 +  rand(noneurons,1)*(0.2);
		cell_parameters.g_Kdr_s  = cell_parameters.g_Kdr_s  -3 + rand(noneurons,1)*6;



	case 'randomized3'

		cell_parameters = defneurons(noneurons);
		
		cell_parameters.g_CaL    = cell_parameters.g_CaL   - .3  + rand(noneurons,1)*1.5;
		cell_parameters.g_int 	 = cell_parameters.g_int   + .07 + rand(noneurons,1)*.4;
		cell_parameters.g_h 	 = cell_parameters.g_h 	   +  rand(noneurons,1)*(1);
		cell_parameters.g_K_Ca   = cell_parameters.g_K_Ca  +  rand(noneurons,1)*10;       
		cell_parameters.g_ld     = cell_parameters.g_ld    +  rand(noneurons,1)*(-0.001);
		cell_parameters.g_la     = cell_parameters.g_la    +  rand(noneurons,1)*(-0.001);
		cell_parameters.g_ls     = cell_parameters.g_ls    +  rand(noneurons,1)*(-0.001);
		cell_parameters.p1       = cell_parameters.p1      - .15 +  rand(noneurons,1)*(0.1);
		cell_parameters.g_Kdr_s  = cell_parameters.g_Kdr_s  -3 + rand(noneurons,1)*6;



	case 'randomized4'

		cell_parameters = defneurons(noneurons);
		
		cell_parameters.g_CaL    = cell_parameters.g_CaL  -.6 + rand(noneurons,1)*1.5;
		cell_parameters.g_int 	 = cell_parameters.g_int   + .07 + rand(noneurons,1)*.4;
		cell_parameters.g_h 	 = cell_parameters.g_h 	   +  rand(noneurons,1)*(1);
		cell_parameters.g_K_Ca   = cell_parameters.g_K_Ca  +  rand(noneurons,1)*10;       
		cell_parameters.g_ld     = cell_parameters.g_ld    +  rand(noneurons,1)*(-0.002);
		% cell_parameters.g_la     = cell_parameters.g_la    +  rand(noneurons,1)*(-0.001);
		cell_parameters.g_ls     = cell_parameters.g_ls    +  rand(noneurons,1)*(-0.002);
		cell_parameters.p1       = cell_parameters.p1      - .15 +  rand(noneurons,1)*(0.1);
		cell_parameters.g_Kdr_s  = cell_parameters.g_Kdr_s  -3 + rand(noneurons,1)*6;


	case 'permuted'

		p1 = [.5:.1:1.1]; 		% CalciumL - conductance range
		p2 = [.0];      	    % g_h_s
		p3 = [.11 3]; 		% g_int
		p4 = [.12:.12:.48];      	% g_h
		p5 = [-38];       	% V_h
		p6 = [45 55];		% Ca act Potassium: not voltage dependent 
		p7 = [4.5];
		p8 = [.013];    % leak
		[p{1} p{2} p{3} p{4} p{5} p{6} p{7} p{8}] = ndgrid(p1,p2,p3,p4,p5,p6,p7,p8);

		Plist = [p{1}(:) p{2}(:) p{3}(:) p{4}(:) p{5}(:) p{6}(:) p{7}(:) p{8}(:)]; 

		psweepnoneurons = length(p{1}(:));

		cell_parameters = defneurons(psweepnoneurons);
		
		cell_parameters.g_CaL    = p{1}(randi(noneurons,noneurons,1));
		cell_parameters.g_h_s    = p{2}(randi(noneurons,noneurons,1));
		cell_parameters.g_int 	 = p{3}(randi(noneurons,noneurons,1));
		cell_parameters.g_h 	 = p{4}(randi(noneurons,noneurons,1));
		cell_parameters.V_h 	 = p{5}(randi(noneurons,noneurons,1));
		cell_parameters.g_K_Ca   = p{6}(randi(noneurons,noneurons,1));       
		cell_parameters.g_CaH    = p{7}(randi(noneurons,noneurons,1));     % High-threshold calcium
		cell_parameters.g_ld     = p{8}(randi(noneurons,noneurons,1));
		


	case 'devel'

	% CAL    IH     g_int    g_l     p1     g_K_Ca    arb    freq_each     ampl     meanVm     spks    supth
 %    ___    ___    _____    ____    ___    ______    ___    _________    ______    _______    ____    _____
 %    8      1.2    0.3      0.01    0.1    50        0.5    0.006676     72.685    -63.776    5       0.236
 % 	  4      1.2    0.05     0.01    0.1    40        0.5    0.0094888    14.896    -55.374    0       0.465

		cell_parameters = defneurons(noneurons);
		cell_parameters.g_CaL    			= ones(noneurons,1)* 4;
		cell_parameters.g_int 	 			= ones(noneurons,1)* .05;
		cell_parameters.g_h 	 			= ones(noneurons,1)* 1.2;
		cell_parameters.g_K_Ca   			= ones(noneurons,1)* 40;
		cell_parameters.g_ld     			= ones(noneurons,1)* .01;
		cell_parameters.g_la     			= ones(noneurons,1)* .01;
		cell_parameters.g_ls     			= ones(noneurons,1)* .01;
		cell_parameters.p1    	 			= ones(noneurons,1)* .1;
		cell_parameters.arbitrary			= ones(noneurons,1)* .5;
		
		
	case 'cellset_devel'
		 if exist('dev_cells_pspace.mat')
		 	load dev_cells_pspace
		 else
		 	disp('couldnt find the cell set')
		 end

		rtan = resultstable.allneurons;
		sel_cel_idx =find( table2array(rtan(:,9))>=5 & table2array(rtan(:,9))<20 & table2array(rtan(:,8))>0.004 & table2array(rtan(:,10))<-55 );

		sel_cel_idx = sel_cel_idx(1:noneurons);
		 
		for ff = fields(st_st.cellParameters)'
			str = ['cell_parameters.' ff{1} ' = simresults{1}.cellParameters.' ff{1} '(sel_cel_idx);'];
			eval(str)
		end


	case 'cellset_vanilla'
		 if exist('cellset_vanilla_3.mat')
		 	load cellset_vanilla_3
		 else
		 	disp('couldnt find the cell set')
		 	keyboard
		 	return
		 end

		if noneurons > length(sel_cel_idx)
			disp('Insufficient number of neurons with criterion. Padding with random picks.')
			disp('max num of cells')
			length(sel_cel_idx)

			sel_cel_idx = sel_cel_idx(randi(length(sel_cel_idx)	,1,noneurons));

		else
			sel_cel_idx = sel_cel_idx(randperm(noneurons));
		end

		for ff = fields(simresults{1}.cellParameters)'
			str = ['cell_parameters.' ff{1} ' = simresults{1}.cellParameters.' ff{1} '(sel_cel_idx);'];
			eval(str)
		end


	case 'param_sweep'
		

		p1 = [.5:.1:1.5]; 		% CalciumL - conductance range
		p2 = [.1:.1:3]; 		% g_int
		
		[p{1} p{2}  ] = ndgrid(p1,p2);

		Plist = [p{1}(:) p{2}(:) ]; 

		psweepnoneurons = length(p{1}(:));

		cell_parameters = defneurons(psweepnoneurons);
		
		cell_parameters.g_CaL    = p{1}(:);
		cell_parameters.g_int 	 = p{2}(:);

	case 'psweep_gh_gcal'

		% parameter ranges
		p1 = [.5:.1:1.1]; 		% CalciumL - conductance range
		p2 = [.12:.12:1.2];      	% g_h
		p3 = [45 55];		% Ca act Potassium: not voltage dependent 
		p4 = [4.5 5.5];

		[p{1} p{2} p{3} p{4}] = ndgrid(p1,p2,p3,p4);

		Plist = [p{1}(:) p{2}(:) p{3}(:) p{4}(:) ]; 

		psweepnoneurons = length(p{1}(:));

		cell_parameters = defneurons(psweepnoneurons);
		
		cell_parameters.g_CaL    = p{1}(:);
		cell_parameters.g_h 	 = p{2}(:);
		cell_parameters.g_K_Ca   = p{3}(:);       
		cell_parameters.g_CaH    = p{4}(:);     % High-threshold calcium
	


	otherwise

		disp('Legacy parameter or celltypes case not found.')


end

% [=================================================================]
%  parameter list
% [=================================================================]

string = [];

for param = Pnames'

	string = [string 'cell_parameters.' param{1} '(:) '];

end



% [=================================================================]
%  randomizer
% [=================================================================]

if addrand
	cell_parameters = jitter_cell_parameters(cell_parameters,.05);
end
		

if gapcompensation
		% compensation for gap junctions
		
		% cell_parameters.g_CaL = cell_parameters.g_CaL + .1*gapcompensation;

		cell_parameters.g_ld  = cell_parameters.g_ld  - 0.001*gapcompensation;
		cell_parameters.g_la  = cell_parameters.g_la  - 0.001*gapcompensation;
		cell_parameters.g_ls  = cell_parameters.g_ls  - 0.001*gapcompensation;
end

if nogapcompensation~=0

		% cell_parameters.g_CaL = cell_parameters.g_CaL - .1;
		% cell_parameters.g_int = cell_parameters.g_int + 0.03;
		% cell_parameters.g_ld  = cell_parameters.g_ld  + 0.001*nogapcompensation
		cell_parameters.g_ls  = cell_parameters.g_ls  + 0.001*nogapcompensation
		% cell_parameters.g_K_Ca= cell_parameters.g_K_Ca - 10;
end




cell_parameters.Plist = eval([ '[' string  ']' ] );
cell_parameters.Pnames = Pnames;




% [=================================================================]
%  create default neurons
% [=================================================================]

function cell_parameters = defneurons(noneurons)

O = ones(noneurons,1);

%% Cell properties

% Capacitance
cell_parameters.C_m    =   1 .* O; % mF/cm^2 
 

% [================================================]
% 		 conductances
% [================================================]
% Defaults

% Somatic conductances (mS/cm2)
cell_parameters.g_CaL    =  1.1   .*O; 
cell_parameters.g_Na_s   =  150   .*O;      % Sodium
cell_parameters.g_Kdr_s  =  9.0   .*O;      % Potassium
cell_parameters.g_K_s    =  5     .*O;      % Potassium
cell_parameters.g_ls     =  0.016 .*O;      % Leaks
    
% Dendritic conductances (mS/cm2)
cell_parameters.g_K_Ca   =  35      .*O;       % Potassium: not voltage dependent 
cell_parameters.g_CaH    =  4.5     .*O;     % High-threshold calcium
cell_parameters.g_ld     =  0.016  .*O;   % Leak
cell_parameters.g_h      =  .12    .*O;    % H current .12
cell_parameters.g_h_s    =  .12    .*O;    % H current, somatic


% Axon hillock conductances (mS/cm2)
cell_parameters.g_Na_a   =  240		.*O;      % Sodium
cell_parameters.g_K_a    =  240		.*O;      % Potassium
cell_parameters.g_la     =  0.016	.*O;      % Leak
    
% Cell morphology
cell_parameters.p1     = 0.25		.*O;        % Cell surface ratio soma/dendrite
cell_parameters.p2     = 0.15 		.*O;        % Cell surface ratio axon(hillock)/soma

cell_parameters.g_int  = 0.13		.*O;        % Cell internal conductance 


% synaptic conductances

cell_parameters.gbar_gaba_dend  = O*.25;
cell_parameters.gbar_gaba_soma  = O*.5;
cell_parameters.gbar_ampa_soma 	= O*.1;
cell_parameters.gbar_ampa_dend 	= O*.1;


%% Reversal potentials
cell_parameters.V_Na =  55 .* O;       % Sodium
cell_parameters.V_K  = -75 .* O;       % Potassium
cell_parameters.V_Ca = 120 .* O;       % Calcium
cell_parameters.V_h  = -43 .* O;       % H current
cell_parameters.V_l  =  10 .* O;       % Leak

cell_parameters.V_gaba_dend = -70 .*O; % from Devor and Yarom, 2002
cell_parameters.V_gaba_soma = -63 .*O; % from Devor and Yarom, 2002
cell_parameters.V_ampa_soma = 0 	  .*O; % from Cian McDonnel et al 2012
cell_parameters.V_ampa_dend = 0 	  .*O; % from Cian McDonnel et al 2012


cell_parameters.arbitrary = 1 .* O;


