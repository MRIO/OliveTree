 %    _   __________________  _____ ____  ___   ____________
%    / | / / ____/_  __/ __ \/ ___// __ \/   | / ____/ ____/
%   /  |/ / __/   / / / /_/ /\__ \/ /_/ / /| |/ /   / __/
%  / /|  / /___  / / / ____/___/ / ____/ ___ / /___/ /___
% /_/ |_/_____/ /_/ /_/    /____/_/   /_/  |_\____/_____/


% [=================================================================]
%  script parameters
% [=================================================================]
compute_simulations        	  = 1;
make_summary_tables 		  = 1;

% [================================================]
%  common parameters
% [================================================]
delta = .02; 
Fs = 1000; % NOTE: sampling rate for output != delta

if ~exist('simtime') ;simtime = 10;end
if ~exist('nameprefix');nameprefix =[];end
if ~exist('seed');seed = 0; end
gpu = 1;
rng(seed,'twister');
cell_function = 'vanilla'; % 'devel'

savehist = 1;
	plotfigs = 0;



% [=================================================================]
%  parameter space grid
% [=================================================================]

if ~exist('p13','var')
	conntype = 'iso'
	% 8 Dimensional GRID: parameter ranges
	p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
	p2  = [0];		% sametoall (noise correlation) of ornstein uhlenbeck process
	p3  = [0];			% noiseamp of ornstein uhlenbeck process
	p4  = [0.04 eps]; % average gap leak per cell
	p5  = [2];			% single cell connection radius 
	p6  = [20];	% clustersize (iff conntype='cluster')
	p7  = [.9];			% intraclusterP (iff conntype='cluster')
	p8  = [.1];			% extraclusterP (iff conntype='cluster')
	p9  = [8];			% mean number of connections
	p10 = [6];			% N, size of N x N network
	p11 = [2];			% depth in Z
	p12 = [0];		    % whether using gap compensation in cells
	p13 = [.5];			% variability of noise injected

end

	[p{1} p{2} p{3} p{4} p{5} p{6} p{7} p{8} p{9} p{10} p{11} p{12} p{13}] = ...
						ndgrid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13);

	Plist = ...
		[p{1}(:) p{2}(:) p{3}(:) p{4}(:) p{5}(:) p{6}(:) p{7}(:) p{8}(:) p{9}(:) p{10}(:) p{11}(:) p{12}(:) p{13}(:)]; 

	fieldnames = {...
			'tau';
			'sametoall';
			'noiseamp';
			'gap';
			'radius';
			'clustersize';
			'intraclusterP';
			'extraclusterP';
			'meannoconn';
			'numneurons'
			'depth';
			'gapcomp';
			'noisesig'}';



% for connectivity
plotconn  = 0;
normleak  = 1;
randomize = 1; %!check
scaling   = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;

I_app = [];

gapcompensation_f = @(x)(40*x)-1.6;




% [=================================================================]
%  readouts
% [=================================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};

to_report = vsoma;



%===========================================================%
%     ____  __  ___   __   _____ ______  ________
%    / __ \/ / / / | / /  / ___//  _/  |/  / ___/
%   / /_/ / / / /  |/ /   \__ \ / // /|_/ /\__ \ 
%  / _, _/ /_/ / /|  /   ___/ // // /  / /___/ / 
% /_/ |_|\____/_/ |_/   /____/___/_/  /_//____/  
% 
%===========================================================%            


% [===========================================================================================================]

if compute_simulations


  for simcount = 1:size(Plist,1)

  		if simcount == 1
  			saveappliednoise = 1;
  		else
  			saveappliednoise = 0;
  		end

		tau 			= Plist(simcount,1);
		sametoall 		= Plist(simcount,2);
		noiseamp  		= Plist(simcount,3);
		gap 			= Plist(simcount,4);
		radius 			= Plist(simcount,5);
		clustersize 	= Plist(simcount,6);
		intraclusterP 	= Plist(simcount,7);
		extraclusterP	= Plist(simcount,8);
		meannoconn 		= Plist(simcount,9);
		numneurons 		= Plist(simcount,10);		
		depth   		= Plist(simcount,11);
		gapcomp 		= Plist(simcount,12);
		noisesig 		= Plist(simcount,13);

		netsize = [depth numneurons numneurons];
		noneurons = prod(netsize);
		clusterparameters = [1 clustersize intraclusterP extraclusterP];

   	 	displaytext = { [fieldnames{:}] ; [ num2str(Plist(simcount,:))] } ;


   	 	rng(seed,'twister')
		switch conntype
			case 'iso'
				W  = createW('3d_chebychev', netsize, radius, scaling, randomize, plotconn, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
			case 'cluster'
				W  = createW('3d_chebychev', netsize, radius, scaling, randomize, plotconn, maxiter, meannoconn, somatapositions, symmetrize, clusterparameters, normleak);
			case 'onecluster'
				W  = createW('one_cluster', netsize, radius, scaling, randomize, plotconn, maxiter, meannoconn, somatapositions, symmetrize, clusterparameters, normleak);
		end

		% [=================================================================]
		%  Neurons
		% [=================================================================]


		rng(seed,'twister')
		if gap>eps
			neurons = createDefaultNeurons(noneurons,'celltype', 'randomized');
		else
			neurons = createDefaultNeurons(noneurons,'celltype', 'randomized','nogapcompensation',gapcomp);
		end

	   	gnoise = [1/tau noisesig noiseamp seed];


   	 	% [================================================]
   	 	%  COMPUTE!
   	 	% [================================================]
   	 	

		transients{simcount} = IOnet('cell_parameters', neurons , ...
		   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W.W*gap ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu ,  ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', sametoall,'saveappliednoise',saveappliednoise,...
		   	'displaytext', displaytext);


		% [================================================]
		%  finish up
		% [================================================]
		
		   transients{simcount}.scriptParameters.Plist = Plist(simcount,:);
		   transients{simcount}.scriptParameters.fieldnames = fieldnames;
		   transients{simcount}.scriptParameters.cellPlist = neurons.Plist;
		   transients{simcount}.W = W;

		   transients{simcount}.spikes = spikedetect(transients{simcount});

		   if not(savehist)
		   		transients{simcount}.networkHistory = [];
		   end

			% eval(['save noisecorr_dataset_' date])

			fname = [nameprefix '_netpspace' num2str(simcount) '_' conntype '_' num2str(simtime) '_' date ];

			eval(['save ' fname ' -v7.3'])
			eval(['!rm ' nameprefix '_netpspace' num2str(simcount-1) '_' conntype '_' num2str(simtime) '_' date '.mat' ])

	end

end
   

% [===========================================================================================================]

% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'), legend(num2str(Plist)),title('V vs Calcium\_r')

%============================= Make Table with Parameters and Results ==============================%


if make_summary_tables
	ci = [0.05 .5 .95];
	tslice = [300:1300]; % only used for frequency estimations

	sims			= [1:size(Plist,1)]';
	tau 			= Plist(:,1);
	sametoall 		= Plist(:,2);
	noiseamp  		= Plist(:,3);
	gap 			= Plist(:,4);
	radius 			= Plist(:,5);
	clustersize 	= Plist(:,6);
	intraclusterP 	= Plist(:,7);
	extraclusterP	= Plist(:,8);
	meannoconn 		= Plist(:,9);
	numneurons 		= Plist(:,10);
	depth	  		= Plist(:,11);
	gapcomp  		= Plist(:,12);
	noisesig  		= Plist(:,13);

	PTable = table(sims, tau, sametoall, noiseamp, noisesig, gap, radius, clustersize, intraclusterP, extraclusterP, meannoconn, depth, numneurons, gapcomp );

	 for s = 1:size(Plist,1)
				

		T = transients{s};
		V = T.networkHistory.V_soma;
		CAL = T.cellParameters.g_CaL;
		IH  = T.cellParameters.g_h;

		ampl(s,:) 		= quantile(max(V, [], 2) - min(V, [], 2), ci);
		meanVm(s,:) 	= quantile(mean(V, 2), ci);
		caL(s,:) 		= quantile(CAL, ci);

		K 		= measureGlobalSync(transients{s},'plotme', 0,'duration',tslice); %, 'duration', tslice
		spikes  = spikedetect(transients{s});

		FO_sync (s, :)  = K.stats.firstordersync;
		SO_sync (s, :)  = K.stats.secondordersync;
		All_sync(s, :)  = K.stats.overallsync;

		try
			freq(s, :) = quantile(K.frequency, ci);
		catch
			freq(s, :) = 0;
		end

		pop_r(s, 1)  = spikes.popfrequency;
		prop_f(s, 1) = spikes.propspkneurons;

		W = transients{s}.networkParameters.connectivityMatrix;
			W(W==0) = [];
			gaps(s,:) = quantile(W, ci) ;

		W = transients{s}.networkParameters.connectivityMatrix;
			nconn(s,:) = quantile(sum(W>0), ci) ;

	end

	RTable = table(gaps, pop_r, prop_f,  ampl, meanVm, caL, freq, nconn, FO_sync, SO_sync, All_sync);

			fname = [nameprefix '_netpspace' num2str(simcount) '_' conntype '_' num2str(simtime) '_' date 'results'];

			eval(['save ' fname ' -v7.3'])
			

end



% [=================================================================]
%  plot traces
% [=================================================================]

if 0

		tslice = [1:simtime];

	figure
		imagesc(transients{simcount}.networkHistory.V_soma,[-80 -20]), colorbar
		set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
		title([num2str(transients{simcount}.spikes.popfrequency) ' Hz'])

	figure
		ca = axis;
		set(0,'defaultaxescolororder', linspecer(length(Plist)));
		p = plot([1:simtime], transients{simcount}.networkHistory.V_soma');
		legend(num2str(neurons.Plist))
	

end

