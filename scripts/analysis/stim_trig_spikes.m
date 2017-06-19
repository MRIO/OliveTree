% stim_triggered_spikes.m


% periodic_ampa_8_cluster_gallop_1000_4_03-Jun-2016.mat
% periodic_ampa_8_iso_gallop_50000_4_31-May-2016.mat
% periodic_ampa_8_iso_spont_50000_4_01-Jun-2016
% periodic_ampa_8_cluster_1Hz_50000_4_03-Jun-2016.mat
% periodic_ampa_8_cluster_gallop_50000_4_03-Jun-2016.mat
% periodic_ampa_8_cluster_spont_50000_4_04-Jun-2016.mat
% F = 'periodic_ampa_8_iso_1Hz_50000_4_02-Jun-2016.mat';

% F1 = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_12-Jun-2016.mat';
% F2 = 'periodic_ampa_2_iso_0.04_spont_50000_2_12-Jun-2016.mat';
% F3 = 'periodic_ampa_2_iso_0.04_gallop_50000_2_12-Jun-2016.mat';
% F = 'periodic_ampa_2_iso_0.04_1Hz_50000_2_20-Jun-2016.mat';
% F = 'periodic_ampa_2_iso_0.04_gallop_50000_2_20-Jun-2016.mat';
% F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_1Hz_50000_2_28-Jun-2016.mat'
% F = 'periodic_ampa_moreoscillations_nocorr_2_iso_0.04_gallop_50000_2_29-Jun-2016.mat'
	
% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_gallop_50000_4_25-Sep-2016.mat';
% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0.04_1Hz_50000_4_25-Sep-2016.mat'; % seemingly wrong gap junction neighborhood
% F1 = 'periodic_ampa_replay_06_12_16_4_iso_0_1Hz_50000_4_25-Sep-2016.mat';

addpath('/Users/M/Synced/Titan/Bench2/periodic_ampa/')
addpath('/Users/M/Synced/Titan/Bench2/')
addpath('/Users/M/Synced/Titan/Bench/')

addpath('/Users/M/Projects/Experiments/Olive/model/simresults/periodic_ampa')
addpath('/Users/M/Synced/Projects/Experiments/Olive/model/simresults/periodic_ampa')

F1 = 'periodic_ampa_replay_06_12_16_with_spont_gaptest8_iso_1Hz_50000_4_17-Jan-2017.mat'; runs = [5:8];





% [=================================================================]
%  analysis to run
% [=================================================================]


spontaneous = 0;
% npulses = 3;
profilesim = 0;
plotstruct = 0;
plotstuff = 0;
plot_selected_neurons = 0;
computerasters = 1;
	partialcorrelations_and_responsescatters = 1;
calculate_xcorrs = 0;
stimtrigwaves = 1;


trigger = 1; % perturbation (according to 'pert')
cellselection = [105];
cellselection = [];
% cellselection = [105 115] ;  %[7 35 55 115]
% cellselection = [1:200];

% [=================================================================]
%  defaults
% [=================================================================]

fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');

[az el] = view(-30,30);

%=============================gather data==============================%

if not(exist('Joinedsim'))
	load (F1)
	Joinedsim{1}  = joinsim(simresults,runs); 
end

sim = Joinedsim{1};
simtime = sim.duration;
sims = 1;
allspikes = sim.spikes.spikes;
triggers = sim.perturbation.triggers{trigger}+35; % Subtract offset to shift triger to first ampa pulse
triggers(end-3:end) = [];
ntrigs = length(triggers);


spks = sim.spikes;

spkfreq = (spks.spikespercell+eps)/sim.duration*1e3;


Q = quantile(spks.spikespercell,[.25 , .5, .75])

hicount = find(spks.spikespercell>Q(3),1,'last')

triggers = sim.perturbation.triggers{trigger};

W = sim.networkParameters.connectivityMatrix;
coords = sim.W.coords;

neurons = sim.cellParameters;
CAL = sim.cellParameters.g_CaL;


gapneighborhood = full(sum(triu(W)))+eps;	
noneurons = prod(netsize);
%===========================================================%
if profilesim
	profiled{1} = profile_sim(Joinedsim{1},'tslice',[1000:5000]);
	NDscatter(profiled{1}.allneurons(:,{'g_CaL'; 'spks'}))
end




if plotstruct 

	        depth   = [1:netsize(1)];
	        breadth = [1:netsize(2)];
	        height  = [1:netsize(3)];

	        % # compute adjacency matrix
	        [XX YY ZZ] = meshgrid(depth,breadth,height);
	        XX = XX(:); YY = YY(:); ZZ = ZZ(:);

	        idx = CAL;

	        plotnetstruct(W,XX,YY,ZZ,idx,[ 1 1 1])
	        create_input_mask(netsize, 'dist_to_point', 'radius', 2,'cell_coordinates', [XX YY ZZ],'projection_center', netsize/2,'synapseprobability',1 ,'plotme',1);

	        thiscmap = [ .7 .7 .7; .2 .2 1];
			colormap(thiscmap)
end


if computerasters

	for i = 1:1:prod(netsize) % for all neurons



		if spks.spikespercell(i) >= 1 & ~spontaneous

			if plot_selected_neurons & ismember(i, cellselection) 
				figure(1000)
				clf
			end

			% try
				% trigrast = ETR(triggers, allspikes{i} , 'waves', sim.networkHistory.V_soma(i,:),'bin', 10, 'span', 3000,'plotQ',ismember(i, cellselection) ,'markertype', 'none');
			% catch
				trigrast = ETR(triggers, allspikes{i} ,'bin', 1, 'span', 1000,'plotQ',ismember(i, cellselection) );
			% end
			spksininterval = @(tr) length(find(trigrast.eventTriggeredRaster{tr}>0 & trigrast.eventTriggeredRaster{tr}<30 ));


			% resp(i,1) = sum(arrayfun(spksininterval, [1:length(triggers)])>=1 )/ length(triggers); % if we're checking the probability of any spike in window
			resp(i,1) = sum(arrayfun(spksininterval, [1:length(triggers)])>=1 )/ length(triggers); % if we're counting their number

			collectedHistogram(i,:) = trigrast.histogram{2};

			if plotstuff  & ismember(i, cellselection) 

				figure(2000)
				clf
				trigrast2 = ETR(triggers(1:ntrigs), allspikes{i} , 'waves', sim.networkHistory.V_soma(i,:),'bin', 10, 'span', 3000,'plotQ',1);
				
				TWVS = trigrast2.eventTriggeredWaveforms;
				TWVS(:,find(isnan(sum(TWVS))))=[];


				figure(2001)
				clf
				imagesc(TWVS',[-70 -10]), axis xy
				[SPKi SPKj v_ ] = find(TWVS>-10);
				line(SPKi', SPKj','marker', '.','linestyle', 'none','color', 'k','markersize',20)
				line(SPKi', SPKj','marker', '+','linestyle', 'none','color', 'w','markersize',5)
				colormap bone(64)




				txt = {['CaT: ' num2str(sim.cellParameters.g_CaL(i)) ]; ['cell: ' num2str(i)] ; ['nspikes: ' num2str(spks.spikespercell(i))] }
				tt = text(1800,0,txt,'backgroundcolor',[1 1 1]);
				tt.VerticalAlignment = 'baseline';
				% axis off;

			
			end			

		end


		if plot_selected_neurons  & ismember(i, cellselection) 
				
			figure(20000)
				[az el] = view(-30,30);
				subplot(1,2,1)
				cla
				scatter3( sim.cellParameters.g_CaL, gapneighborhood, spkfreq , 25, gapneighborhood,'filled'); view([az el])
				l1 = line( sim.cellParameters.g_CaL(i), gapneighborhood(i), spkfreq(i), 'marker', '+','color', 'r', 'markersize',50,'linewidth',4);
				
				XXX = [sim.cellParameters.g_CaL' ; sim.cellParameters.g_CaL'];
				YYY = [gapneighborhood ;gapneighborhood];
				ZZZ = [spkfreq*0 ; spkfreq]; 
				line(XXX, YYY, ZZZ,'color', 'k');


				xlabel('Ttype conductance'), ylabel( 'gapneighborhood'),zlabel( 'spike frequency (Hz)')
				view(az,el)
				axis tight


				subplot(1,2,2)
				cla
				scatter3 (coords(:,1), coords(:,2), coords(:,3), 100,spks.spikespercell'+eps,'filled'); view([az el]); axis equal
				
				l2 = line(coords(i,1), coords(i,2), coords(i,3), 'marker', '.','color', 'r', 'markersize',70);
		


				axis equal
				view(az,el)
		



				drawnow

				pause
		end
	
	% saveallfigs('prefix', ['cell_' num2str(i) '_'])
	


	end
end

% [=================================================================]
%  compute partial correlations of spiking responses
% [=================================================================]
if partialcorrelations_and_responsescatters
	 	% resp     = sum(collectedHistogram(:,1001:1030),2)./30*1e3 ; % response frequency in post stimulus window
		CAL = sim.cellParameters.g_CaL;
		IH  = sim.cellParameters.g_h;
		IINT = sim.cellParameters.g_int;
		GAMP = sim.cellParameters.gbar_ampa_soma;
		GN 	 = gapneighborhood';
		MASK = sim.perturbation.mask{1};

		T = [CAL IH IINT MASK GAMP GN spkfreq' resp];		

		[CORR PVAL] = partialcorr(T);

		figure
		subplot(121)
		% CORR(find(triu(ones(size(CORR))))) = -.01;
		imagesc(CORR,[-1 1]), colorbar
		set(gca,'xticklabel', {'g CaL' 'g h' 'g s-d' 'mask' 'g_ampa' 'gap leak' 'spk freq' 'resp freq'})
		set(gca,'yticklabel',{'g CaL' 'g h' 'g s-d' 'mask' 'g_ampa' 'gap leak' 'spk freq' 'resp freq'})
		title('partial correlation')

		subplot(122)
		imagesc(PVAL<0.05,[0 1]), colorbar
		set(gca,'xticklabel', {'g CaL' 'g h' 'g s-d' 'mask' 'g_ampa' 'gap leak' 'spk freq' 'resp freq'})
		set(gca,'yticklabel', {'g CaL' 'g h' 'g s-d' 'mask' 'g_ampa' 'gap leak' 'spk freq' 'resp freq'})
			title('pval<0.05 (0-no, 1-yes)')

		figure
		scatter(CAL, resp , 100, MASK,'filled')
		xlabel('CaT conductance (ms/cm^2)')
		ylabel('response frequency (Hz)')

		figure
		scatter(CAL, GAMP , resp*500, MASK,'filled')
		xlabel('CaT conductance (ms/cm^2)')
		ylabel('AMPA conductance')

		figure
		scatter3( sim.cellParameters.g_CaL, gapneighborhood, spkfreq , 30, resp*10000,'filled'); view([az el])
		title('response probability')
		xlabel('CaT conductance (mS/cm^2)')
		ylabel('gap neighborhood (mS/cm^2)')
		zlabel('spike frequency (Hz)')
				
		figure
		scatter3( sim.cellParameters.g_CaL, gapneighborhood, resp , spkfreq*10, MASK,'filled'); view([az el])
		title('response probability')
		xlabel('CaT conductance (mS/cm^2)')
		ylabel('gap neighborhood (mS/cm^2)')
		zlabel('response frequency (Hz)')
		
		figure
		scatter(CAL, resp, spkfreq*100 , MASK,'filled')
		xlabel('CaT conductance (ms/cm^2)')
		ylabel('spike frequency (Hz)')

		figure
		scatter(CAL, spkfreq , resp*50, MASK,'filled')
		xlabel('CaT conductance (ms/cm^2)')
		ylabel('spike frequency (Hz)')



	end


if calculate_xcorrs

	% calculate the cross corr between direct neighbors
	[i_ j_ v_] = find(triu(W));
	pairs = [i_ j_ v_]; 
	if isempty(pairs)
		[i_ j_ v_] = find(triu(ones(size(W))));
		pairs = [i_ j_ v_];
		pairs = pairs(1:100,:)
	end

	s = 0;
	for p = 1:length(pairs)
		% p

	% calculate the cross corr between distant neighbors	

		if spks.spikespercell(pairs(p,1)) > 20 & spks.spikespercell(pairs(p,2)) > 20
			s = s+1;
			% [spks.spikespercell(pairs(p,1))  0 spks.spikespercell(pairs(p,2)) ];
			
			resu{p} = Xcorr_dirac(spks.spikes{pairs(p,1)} , spks.spikes{pairs(p,2)});
			m = length(resu{p}.sup)/2;
			X(s,:) = resu{p}.xcorr(m -1000:m+1000 );

		% figure(1000)
			% plot(resu{p}.sup, resu{p}.xcorr)
			% xlim([-1000 1000])
			% title({num2str(v_(p)); [num2str(spks.spikespercell(pairs(p,1))) ' ' num2str(spks.spikespercell(pairs(p,2))) ]} );
			% drawnow
		% pause(.1)
			
		

		end
	end
end


if 0

	if ~exist('HT','var')
	 HT = measureGlobalSync(sim3D, [1:120000],1);
	 U = HT.hilbert.hilbert;
	end
	 

	 triggers = sim3D.perturbation.triggers{2};

	[H x] = hist(HT.hilbert.hilbert,100);

	Z = zeros(100,200);
	 for t = triggers(1:100)
	 	Z = Z + H(:, t-100:t+99);
	 end

	imagesc(Z)
	title('triggered phase concentration')

	K = zeros(1, 301); n =0; 
	for t = triggers(1:30)
		n = n+1;
		K(n,:) = HT.stats.order_parameter_all(t-150:t+150);
	end


	mask = sim3D.perturbation.mask{2};
	W = sim3D.networkParameters.connectivityMatrix>0;
	S =  sim3D.networkHistory.V_soma;
	fo_neighbors = setdiff(find(W * mask), find(mask));
	neither_nor  = setdiff([1:prod(netsize)], union(fo_neighbors, find(mask)));


	% triggered kuramoto parameter distribution 
	% #todo, make it dependent on phase
	
	n =0; 
	for t = triggers(1:100)
		n = n+1;
		interv = t-150:t+150;
		
		Km = U(find(mask), interv );
		KKm(n,:) = mean(exp(i*bsxfun(@minus, Km, circ_mean(Km))));


		K1 = U(fo_neighbors, interv );
		KK1(n,:) = mean(exp(i*bsxfun(@minus, K1, circ_mean(K1))));

		Ka = U(neither_nor, interv);
		KKa(n,:) = mean(exp(i*bsxfun(@minus, Ka, circ_mean(Ka))));
	end

	n =0; 
	for t = triggers(1:100)
		n = n+1;
		interv = t-150:t+150;
		V1 = S(fo_neighbors, interv);
		VV1(n,:) = mean(V1);

		Va = S(neither_nor, interv);
		VVa(n,:) = mean(Va);	

		Vm = S(find(mask), interv);
		VVm(n,:) = mean(Vm);

		Vnm = S(find(not(mask)), interv);
		VVnm(n,:) = mean(Vnm);	

	end

	figure, plot(VV1'), figure, plot(VVa'), figure, plot(VVm')
	figure, plot([mean(VV1); mean(VVa); mean(VVm)]')



end


if stimtrigwaves

	for n = cellselection
		w = [-1000 1000]; % window around trigger
		rw = 150; %response window

		spks= allspikes{n};
		V = sim.networkHistory.V_soma(n,:);

		spks(spks>=length(V)+w(1)) = [];
		spks(spks<=w(2)) = [];


		trgspks = [];
		for t = triggers'
			trgspks = union(trgspks, spks(find(spks>=t & spks <= t+rw,1,'first' )));
		end
		nontrigspks = setdiff(spks, trgspks);

		trgspks(trgspks>=length(V)+w(1)) = [];
		trgspks(trgspks<=w(2)) = [];

		trwvsnp = [];
		for t = trgspks'
			trwvsnp = [trwvsnp ; V(t+w(1):t+w(2)) ]  ;
		end

		nontrwvsnp = [];
		for t = nontrigspks'
			nontrwvsnp = [nontrwvsnp ; V(t+w(1):t+w(2)) ]  ;
		end

		figure(31415+n)
		clf
		 plot_mean_and_std([w(1):w(2)],trwvsnp, 'color', [1 0 0]);
		 hold on
		 plot_mean_and_std([w(1):w(2)],nontrwvsnp, 'color', [0 0 1]);
		alpha(.7)
		drawnow
		% saveallfigs('prefix', ['meantrigwvs_' num2str(n)])


		

	end


end

gallop = 0;
if gallop

	sim = simresults{1};
	plotstuff = 1;
	triggers = sim.perturbation.triggers{1};

	triggers1 = triggers(diff([0 ; triggers])==388); % T = 400; 4 pulses with 4ms difference
	triggers2 = triggers(diff([0 ; triggers])==238);
	triggers1(end) = []; % prune last stimulation trial
	triggers2(end) = []; % prune last stimulation trial

	allspikes = sim.spikes.spikes;
	binaryspikes = sim.spikes.binaryspikes;

	close all
	for ccc = [50:60];

		unit =ccc;
		ETB = @(t,unit) binaryspikes(unit, t-500:t+500);
		bla1 = arrayfun(ETB, triggers1, unit*ones(size(triggers1)),'uniformoutput', 0)
		bla2 = arrayfun(ETB, triggers2, unit*ones(size(triggers2)),'uniformoutput', 0)
		br1 = cell2mat(bla1');
		br2 = cell2mat(bla2');
		figure
			bar( [-500:500]', [sum(br1) ; sum(br2)]')
		% plot([-500:500], sum(br1)), hold on, plot(sum(br2),'r')
		legend({'short' 'long'})
		title(num2str(unit))


	end



end

