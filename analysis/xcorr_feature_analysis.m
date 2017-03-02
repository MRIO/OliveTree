% xcorr_feature_analysis.m #jochen

numneurons = 200;

path_to_sims = ['/Users/M/Synced/Projects/Experiment/Olive/model/simresults/periodic_ampa/'];
computeselectedxcorr = 1;

sims2p = [1 2 3];
nwins = 1;
plotxcorrs = 1;
	plot_selected_neurons  = 1;
	nselneurons = 20;
	if not(exist('selected_neurons'))
		selected_neurons = 'neighbors';
		% selected_neurons = 'nextneighbors';
		% selected_neurons = 'stimulated';
		% selected_neurons = 'acrosscluster';
		% selected_neurons = 'highgapconn';
		% selected_neurons = 'leastgapconn';
		% selected_neurons = 'all';
	end

if ~exist('loaded'); loaded = 0; end


if 1
	if ~exist('Joinedsim')


		simfiles = {...
		'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017.mat';
		% 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_gallop_50000_4_17-Jan-2017.mat';
		'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_spont_50000_4_17-Jan-2017.mat'
		};

		
			load( simfiles{1});
			Joinedsim{1}  = joinsim(simresults,[1:3]); 
			Joinedsim{2}  = joinsim(simresults,[5:7]); 

			load( simfiles{2});
			Joinedsim{3}  = joinsim(simresults,[1:3]); 	
			Joinedsim{4}  = joinsim(simresults,[5:7]); 			

		simlegend = ...
		{'1Hz iso WT'; '1Hz iso MT'; 'spont iso MT'; 'spont iso WT';}

	end

	loaded = 1;
	simresults = Joinedsim;


end


% periodic iso + more oscillation
if 0
	if ~exist('Joinedsim')
		% clustered
		simfiles = {...
			'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_1Hz_50000_2_20-Jun-2016.mat';
			'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_gallop_50000_2_20-Jun-2016.mat';
			'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_spont_50000_2_20-Jun-2016.mat'
			};

		nruns = 2;
		for i = 1:3
			load([path_to_sims simfiles{i}]);

			Joinedsim{i}  = joinsim(simresults,[1:2]); 
			% Joinedsim{(i-1)*2+2}  = joinsim(simresults,[5:8]);
		end

		simlegend = ...
		{'1Hz iso WT'; 'spont iso WT'; 'gallop iso MT' }

	end

	loaded = 1;
	simresults = Joinedsim;

end




% periodic clustered
if 0
	if ~exist('Joinedsim')
		% clustered
		simfiles = {...
			'periodic_ampa_2_iso_0.04_1Hz_50000_2_20-Jun-2016.mat';
			'periodic_ampa_2_iso_0.04_spont_50000_2_20-Jun-2016.mat';
			'periodic_ampa_2_iso_0.04_gallop_50000_2_20-Jun-2016.mat';
			};

		nruns = 2;
		for i = 1:3
			load([path_to_sims simfiles{i}]);

			Joinedsim{i}  = joinsim(simresults,[1:2]); 
			% Joinedsim{(i-1)*2+2}  = joinsim(simresults,[5:8]);
		end

		simlegend = ...
		{'1Hz iso WT'; 'spont iso WT'; 'gallop iso MT' }

	end

	loaded = 1;
	simresults = Joinedsim;

end



% periodic clustered
if 0 
	if ~exist('Joinedsim')
		% clustered
		simfiles = {...
			'periodic_ampa_2_cluster_1Hz_50000_2_07-Jun-2016.mat';
			'periodic_ampa_2_cluster_1Hz_50000_2_08-Jun-2016.mat';
			'periodic_ampa_2_cluster_spont_50000_2_07-Jun-2016.mat';
			'periodic_ampa_2_cluster_spont_50000_2_08-Jun-2016.mat';
			};

		nruns = 2;
		for i = 1:4
			load([path_to_sims simfiles{i}]);

			Joinedsim{i}  = joinsim(simresults,[1:2]); 
			% Joinedsim{(i-1)*2+2}  = joinsim(simresults,[5:8]);
		end

		
		simlegend = ...
		{'1Hz clustered WT'; '1Hz clustered MT'; 'spont cluster WT'; 'spont cluster MT' }

		
	end

	loaded = 1;
	simresults = Joinedsim;

end

% periodic iso
if 0
	if ~exist('Joinedsim')
		% clustered
		simfiles = {...
				'periodic_ampa_8_cluster_1Hz_50000_4_03-Jun-2016.mat';
				'periodic_ampa_8_cluster_spont_50000_4_04-Jun-2016.mat';
				'periodic_ampa_8_iso_1Hz_50000_4_02-Jun-2016.mat';
				'periodic_ampa_8_iso_spont_50000_4_01-Jun-2016.mat';
				};

		for i = 1:4
			load([path_to_sims simfiles{i}]);

			Joinedsim{(i-1)*2+1}  = joinsim(simresults,[1:4]);  % with gap
			Joinedsim{(i-1)*2+2}  = joinsim(simresults,[5:8]);  % no gap

		end

		simresults = Joinedsim;

		simlegend = ...
		{'1Hz clustered WT'; '1Hz clustered MT'; 'spontaneous cluster WT'; 'spontaneous clustered MT' ; ...
		'1Hz isotropic WT'; '1Hz isotropic MT'; 'spontaneous isotropic WT'; 'spontaneous isotropic MT'}

		
	end

	loaded = 1;

end

if 0

	if ~loaded

		simfiles = ...
			{ 'periodic_ampa_16_iso_1Hz_50000_8_25-May-2016.mat';
			  'periodic_ampa_16_iso_spont_50000_8_23-May-2016.mat';
			  'periodic_ampa_16_iso_spont_50000_8_22-May-2016.mat';
			  'periodic_ampa_16_cluster_1Hz_50000_8_21-May-2016.mat';
			  'periodic_ampa_16_iso_1Hz_50000_8_20-May-2016.mat';
			  'periodic_ampa_16_iso_1Hz_50000_8_19-May-2016.mat';
			  'periodic_ampa_16_cluster_gallop_50000_8_18-May-2016.mat';}

		for i = 1:7
			load([path_to_sims simfiles{i}]);

			Joinedsim{(i-1)*2+1}  = joinsim(simresults,[1:8]); % withgap 
			Joinedsim{(i-1)*2+2}  = joinsim(simresults,[9:16]); %nogap

		end
		simlegend = simfiles;
		simresults = Joinedsim;
		loaded = 1;
	end

		
end


neighbors = retrieveNeuronsByClass(simresults{1}, 'neighbors');
stimulated= retrieveNeuronsByClass(simresults{1}, 'stimulated');

selectedneurons = neighbors;

	% centerneuron_index = 85;

	% CM = simresults{sims2p(1)}.networkParameters.connectivityMatrix;
	% 	CM(find(eye(size(CM))))=0;
	% 	[i j v] = find(CM);
	% 	MN = full([i j v]);
	% 	S = sortrows(MN,3);
	% 	NPmost = S(end-5:end, [1 2]);
	% 	NPleast = S(end-3:end, [1 2]);
	% 	NPall = [NPmost ; NPleast];

			
	% 	NPl = reshape(NPmost,1,[]);

	% 	centerneuron_b = zeros (numneurons,1); centerneuron_b(centerneuron_index) = 1;
	% 	neighbors = find(CM*centerneuron_b); 
	% 	neighbors_b = zeros(numneurons,1);
	% 	neighbors_b(neighbors) = 1;
	% 	nextneighbors = find(CM*neighbors_b); 
		
	% 	S = setdiff(nextneighbors, neighbors);

	% 	stimulated = find(simresults{sims2p(1)}.perturbation.mask{1});

	% 	% remove neurons that don't fire enough spikes
	% 	lowfr = find(sum(simresults{sims2p(1)}.spikes.spikespercell')<20);

		
	% 	switch selected_neurons

	% 		case 'neighbors' % johny and his neighbors
	% 			selectedneurons = setdiff(neighbors, stimulated);
	% 			selectedneurons = setdiff(selectedneurons, lowfr);
	% 			R = randperm(length(selectedneurons));
	% 			nselneurons  = min(nselneurons, length(selectedneurons))
	% 			selectedneurons = [centerneuron_index ; selectedneurons(R(1:nselneurons))];

	% 		case 'nextneighbors'
	% 			selectedneurons = setdiff(nextneighbors, stimulated);
	% 			selectedneurons = setdiff(selectedneurons, lowfr);
	% 			R = randperm(length(selectedneurons));
	% 			nselneurons  = min(nselneurons, length(selectedneurons))
	% 			selectedneurons = selectedneurons(R(1:nselneurons))

	% 		case 'stimulated'
	% 			R = randperm(length(selectedneurons));
	% 			nselneurons  = min(nselneurons, length(selectedneurons))
	% 			selectedneurons = setdiff(stimulated, lowfr);
	% 			selectedneurons = stimulated(R(1:nselneurons));
				
	% 		case 'acrosscluster'
	% 			sel1 = setdiff(find(simresults{1}.W.stats.clusters==1), stimulated)  ;
	% 			sel1 = setdiff(sel1, lowfr);
	% 			sel2 = setdiff(find(simresults{1}.W.stats.clusters==2), stimulated);
	% 			R1 = randperm(length(sel1));
	% 			R2 = randperm(length(sel2));
	% 			nselneurons  = min(nselneurons, length(R1))

	% 			selectedneurons =  [sel1(R1(1)) ; sel2(R1(1:nselneurons/2))];
	% 			selectedneurons = setdiff(selectedneurons, lowfr);

	% 		case 'highgapconn'
	% 			selectedneurons = setdiff(NPmost, stimulated);
	% 			selectedneurons = setdiff(selectedneurons, lowfr);
	% 			R = randperm(length(selectedneurons));
	% 			nselneurons  = min(nselneurons, length(selectedneurons))
	% 			selectedneurons = selectedneurons(R(1:nselneurons));

	% 		case  'leastgapconn'
	% 			selectedneurons = setdiff(NPleast, stimulated);
	% 			selectedneurons = setdiff(selectedneurons, lowfr);
	% 			R = randperm(length(selectedneurons));
	% 			nselneurons  = min(nselneurons, length(selectedneurons))
	% 			selectedneurons = selectedneurons(R(1:nselneurons));

	% 		case 'all'
	% 			selectedneurons = 1:length(CM);
	% 			selectedneurons = setdiff(selectedneurons, lowfr);
	% 			R = randperm(length(selectedneurons));
	% 			selectedneurons = selectedneurons(R(1:nselneurons));


	% 	end
		


%=============================flush==============================%
clear collectX	
clear collectAmpl 
clear collectDelay
clear collectAsym 
%===========================================================%



% [=================================================================]
%  compute correlations from all pairs
% [=================================================================]
if computeselectedxcorr

	

	c = 0;
	for simidx = sims2p
		c = c+1;
		
		thisW  = simresults{simidx}.networkParameters.connectivityMatrix;
		
		
		X{c} = xcorr_summa(simresults{simidx},'nwins',nwins,'plotme',plotxcorrs ,'selectedneurons', selectedneurons );
		% X{c}.winsize = winsize;
		% X{c}.asymmetry = asym;
		% X{c}.xcorr_pairs = XCs;
		% X{c}.spikespercell = sum(VSB');
		% X{c}.pairs = pairs;
		% X{c}.delay = delay;
		% X{c}.amplitude = ampl;
		% X{c}.XcorrNoAc = XnoAC;
		% X{c}.XC = XCs;
		% X{c}.selectedneurons = selectedneurons;


		% find significant pairs (currently we use 1:numneurons as baseline. okay as an estimator.)
		XX = X{c}.xcorr_pairs{1}(:,1:numneurons);
		meanAmp = mean(XX,2);
		stdAmp  =  std(XX,[],2);

		sigs{c} = X{c}.amplitude' > meanAmp+2*stdAmp;
		sig_pairs{c} = find(X{c}.amplitude' > meanAmp+2*stdAmp);
		prop_sig_pairs{c} = length(sig_pairs{c})/size(XX,1);

		if plot_selected_neurons & c == 1 
	        depth   = [1:netsize(1)];
	        breadth = [1:netsize(2)];
	        height  = [1:netsize(3)];

	        % # compute adjacency matrix
	        [XX YY ZZ] = meshgrid(depth,breadth,height);
	        XX = XX(:); YY = YY(:); ZZ = ZZ(:);

	        idx = zeros(length(XX),1);
	        idx(X{c}.selectedneurons) = 1;

	        plotnetstruct(thisW,XX,YY,ZZ,idx,[ 1 1 1])
	        thiscmap = [ .7 .7 .7; .2 .2 1];
			colormap(thiscmap)

	        
	        if W.stats.clusters & 0%simresults{simidx}.W.stats.clusters

	        	C = W.stats.clusters;
				plotnetstruct(thisW,XX,YY,ZZ,C ,[ 1 1 1])
			end

			figure
				plotSelectedNeurons(simresults{simidx},[3000:4000],selectedneurons)
				title(num2str(simidx))
				axis tight

		end

	end

	collectX	 = [];
	collectAmpl  = [];
	collectDelay = [];
	collectAsym  = [];
	G = [];
	sig = [];
	
	for c = 1:length(sims2p)
		collectX(:,c)	 = X{c}.XcorrNoAc; % (find(not(isnan(sum(X{c}.XcorrNoAc)))),:) ;

		collectAmpl  = [collectAmpl ; X{c}.amplitude' ];
		collectDelay = [collectDelay; X{c}.delay' 	 ];
		collectAsym  = [collectAsym ; X{c}.asymmetry' ];
		G = [G ; ones(size(X{c}.amplitude))'*c];
		sig = [sig ; sigs{c}];


	end

	if 1
		figure
			
			axes
			set(gca,'colororder',jet(4))
			plot([-400:400], collectX)
			legend(num2str(sims2p'))
			axis tight
			xlabel('ms')
			ylabel('correlation (coeff)')
			xlabel('aggregate correlation')

			legend(simlegend{sims2p})

			leg = ['nonsig' ; simlegend];


		figure; boxplot([abs(collectAsym )  ], G.*sig, 'notch', 'on'); title('asym'); 		set(gca,'xticklabel', leg)
		figure; boxplot([abs(collectDelay)  ], G.*sig, 'notch', 'on'); title('delay'); 		set(gca,'xticklabel', leg)
		figure; boxplot([	 collectAmpl   	], G.*sig, 'notch', 'on'); title('amplitude'); 	set(gca,'xticklabel', leg)
	end

end





