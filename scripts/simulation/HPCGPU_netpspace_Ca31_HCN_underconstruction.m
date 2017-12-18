% noise_corr_gap_comparisons.m


% conntype = 'cluster'; HPCGPU_netpspace_noise_corr_gap;
% conntype = 'iso'    ; HPCGPU_netpspace_noise_corr_gap;

% [=================================================================]
%  script parameters
% [=================================================================]

if ~exist('conntype');conntype = 'iso';end

compsim = 1; gpu = 1;
	savehist = 1;
summarize = 0;
	calcfreq = 0;
plotsummaries = 0;
plottraces = 0;
plotfigs = 0;
parametertable = 0;
computeselectedxcorr = 0;
	sims2p = [16 24];

rndseed = 0;
rng(rndseed,'twister')

fname = ['pspace_noise_gap_' num2str(simcount) '_' conntype '_' num2str(simtime) '_' date ];


% [=================================================================]
%  simulation parameters
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
steadystate_time = 300;
simtime  = 10000;

delta = .02; 
Fs = 1000; % sampling rate =! delta

% [================================================]
%  net parameters
% [================================================]

netsize = [4 10 10];
	noneurons = prod(netsize);

rd = 3;
meannoconn = 10;

plotconn  = 1;
normleak  = 1;
randomize = 1; %!check
scaling   = 1;
plotthis  = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;


switch conntype
	case 'iso'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
	case 'cluster'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 30 .5 .05], normleak);
end

% [=================================================================]
%  Neurons
% [=================================================================]

rng(rndseed,'twister')
def_neurons = createDefaultNeurons(noneurons,'celltype', 'randomized','gapcompensation',1);
rng(rndseed,'twister')
gapless_neurons = createDefaultNeurons(noneurons,'celltype', 'randomized','nogapcompensation', 1);
Plist = def_neurons.Plist;



% [=================================================================]
%  parameter space grid
% [=================================================================]

gaps = [0.05 0.025 0.001];
noisecorr = [0 .1 .2 .3] ; % proportion of shared noise
noiseamps = [.5 1 ]; %uA/cm^2 -> x .1 nA for a cell with 10000um^2


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

tau = 10;
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
   for gs = gaps
   	for ncorr = noisecorr
   	 for namp = noiseamps


   	 	displaytext = [ conntype '_' num2str(gs) '_' num2str(ncorr)];

   	 	rng(0,'twister')
   	 	gnoise = [theta namp 0 5];
   	 	simcount = simcount+1;

   	 	if gs > 0.01
		   neurons = def_neurons;
		else
		   neurons = gapless_neurons;
		end

		transients{simcount} = IOnet('tempState', st_st.lastState ,'cell_parameters', neurons , ...
		   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W.W*gs ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu ,  ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', ncorr,'saveappliednoise',0,...
		   	'displaytext', displaytext);

		   transients{simcount}.scriptParameters.noiseAmplitude   = namp;
		   transients{simcount}.scriptParameters.noiseCorrelation = ncorr;
		   transients{simcount}.scriptParameters.gapAmplitudes    = gs;
		   transients{simcount}.spikes = spikedetect(transients{simcount});
		   transients{simcount}.Plist = Plist;

		   if not(savehist)
		   		transients{simcount}.networkHistory = [];
		   end

		    if plotfigs
				figure
				imagesc(transients{simcount}.networkHistory.V_soma,[-80 -20]), colorbar
				set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
				title([num2str(transients{simcount}.spikes.popfrequency) ' Hz'])

				figure
				ca = axis;
				set(0,'defaultaxescolororder', linspecer(length(Plist)));
				p = plot([1:simtime], transients{simcount}.networkHistory.V_soma');
				legend(num2str(Plist))
			end

			% eval(['save noisecorr_dataset_' date])

			eval(['save ' fname ' -v7.3'])
				
   	  end
	 end
	end


end

disp(['computed ' num2str(simcount) ' simulations.'])
   

% [===========================================================================================================]

% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'), legend(num2str(Plist)),title('V vs Calcium\_r')

%=============================Make Parameter Table ==============================%


if parametertable

	% pspace
	simcount= 0;
	   for gs = gaps
	   	for ncorr = noisecorr
	   	 for namp = noiseamps
	   	 	simcount = simcount +1;
				   

				sc(simcount,1) = simcount;
				gp(simcount,1) = gs;
				nc(simcount,1) = ncorr;
				na(simcount,1) = namp;

				transients{simcount}.spikes = spikedetect(transients{simcount});

				fr(simcount,1) = transients{simcount}.spikes.popfrequency;

			end
		end
	end
end

Ptable = [sc gp nc na fr ];


% [=================================================================]
%  plot trace matrix
% [=================================================================]

if plottraces & ~isemtpy(transients{ind}.networkHistory)
	ind = 0;
	tslice = [1:simtime];
	pind = 1;
	set(0,'defaultaxescolororder', linspecer(length(transients{1}.Plist)));
	% set(0,'defaultaxescolormap', linspecer(length(transients{1}.Plist)));
	set(0,'defaultfigurecolor', [1 1 1]);
	for gs = gaps
		
	   	for ncorr = noisecorr
	   	
	   		for namp = noiseamps


		   		ind = ind+1;

		   		if transients{ind}.scriptParameters.noiseAmplitude(1) == 1

					figure(1)
					ax1(pind) = subplot(length(gaps), length(noisecorr), pind);
					pl = plot(transients{ind}.networkHistory.V_soma');
					title(['gaps:' num2str(transients{ind}.scriptParameters.gapAmplitudes) 'ncorr:' num2str(transients{ind}.scriptParameters.noiseCorrelation)]) 
					axis tight
					xlim([tslice(1) tslice(end)])	

					figure(2)
					ax2(pind) = subplot(length(gaps), length(noisecorr), pind);
					imagesc(transients{ind}.networkHistory.V_soma,[-80 -20])
					% set(gca,'ytick', [1:56],'yticklabel', num2str(transients{ind}.Plist),'fontsize',8)
					title([num2str(transients{ind}.spikes.popfrequency) ' Hz'])

					pind = pind+1;

				end


			end
		end
	end
	linkaxes(ax1,'x')
	linkaxes(ax2,'x')

end

% [=================================================================]
%  plot correlations
% [=================================================================]

if plotcorrelations
	
	ind = 0; 
	minspkevents = 20;
	spkthreshold = -10;
		
	lags = 300;


	
	if not(exist('X'))
		for gs = gaps
		   	for ncorr = noisecorr
		   		for namp = noiseamps
			   		ind = ind+1;
			   		
			   			t{ind}= transients{ind}.networkHistory.V_soma';
						T{ind} = t{ind}(:,find(sum(t{ind}>spkthreshold,1))>minspkevents);
						ss = T{ind}>spkthreshold;
						% sss = padarray(diff(ss,2), [0 1], ;

						if isempty(T{ind})
							continue
						end
						X{ind} = xcorr(T{ind},lags,'coeff');
						% X{ind} = xcorr(T{ind},lags,'unbiased');
					
				end
			end
		end
	end




	

	% if not(exist('X'))
	% 	for gs = gaps
	% 	   	for ncorr = noisecorr
	% 	   		for namp = noiseamps
	% 		   		ind = ind+1;
			   		
	% 		   			t{ind}= transients{ind}.networkHistory.V_soma';
	% 					T{ind} = t{ind}(:,find(sum(t{ind}>spkthreshold,1))>minspkevents);
	% 					ss = T{ind}>spkthreshold;
	% 					% sss = padarray(diff(ss,2), [0 1], ;

	% 					if isempty(T{ind})
	% 						continue
	% 					end
	% 					X{ind} = xcorr(T{ind},lags,'coeff');
	% 					% X{ind} = xcorr(T{ind},lags,'unbiased');
					
	% 			end
	% 		end
	% 	end
	% end


	pind = 0; ind = 0;
	for gs = gaps
	   	for ncorr = noisecorr
	   		for namp = noiseamps
	   			ind= ind+1;
	
		   		if transients{ind}.scriptParameters.noiseAmplitude(1) == 1
		   			pind = pind +1;
					figure(1)
					ax1(pind) = subplot(length(gaps), length(noisecorr), pind);
					pl = plot([-lags:lags], mean(X{ind}'));
					axis tight
					title(['gaps:' num2str(transients{ind}.scriptParameters.gapAmplitudes) 'ncorr:' num2str(transients{ind}.scriptParameters.noiseCorrelation)]) 
					
					meanxcorr(pind,:) = mean(X{ind}');
					legends{pind} =  ['gaps:' num2str(transients{ind}.scriptParameters.gapAmplitudes) 'ncorr:' num2str(transients{ind}.scriptParameters.noiseCorrelation)];

				end

			end
		end
	end
	figure(2)
	set(0,'defaultaxescolororder',linspecer(6))
	plot([-lags:lags], meanxcorr);
	legend(legends)

end

% [=================================================================]
%  measure sync
% [=================================================================]

window_for_sync_measure = [1001:2000];

if summarize
	simcount= 0;
	ii = 0; jj = 0; kk = 0;
	for gs = gaps
		ii = ii + 1;
	   	for ncorr = noisecorr
	   	jj = jj + 1;
	   	 for namp = noiseamps
	   	 	kk = kk + 1;
	   	 	simcount = simcount+1;

			 Zpropfiring(ii,jj,kk) = transients{simcount}.spikes.propspkneurons;
			 Zpopfreq   (ii,jj,kk) = transients{simcount}.spikes.popfrequency;
			 warning off
			 transients{simcount}.perturbation.mask{1} = ones(1,noneurons);
			 transients{simcount}.perturbation.onsets{1} = 10; 
			 transients{simcount}.perturbation.triggers{1} = 10; 

			 if calcfreq
				 K = measureGlobalSync(transients{simcount},'duration', window_for_sync_measure);
				 
			 	 InstFreq = diff(unwrap(K.hilbert.hilbert)') * Fs/(2*pi);
			 	 [HistInstFreq x] = hist(InstFreq,201); % per cell
			 	 [HistInstFreq_overall x] = hist(InstFreq(:),x); %InstFreq
			 	 
			 	 % surf([1:3000],x, HistInstFreq);
			 	 % title({'peak frequency'  ; ['gap:' num2str(gs) '\gamma:' num2str(ncorr)]} )
			 	 % ylim([-150 150])
			 	 % xlabel('time (ms)')
			 	 % shading interp
			 	 view(-90, 95)
			 	 peakIF = median(InstFreq(:));
			 	 kurtIF = kurtosis(InstFreq(:));


			 	 warning on
				 Zfreq      	(ii,jj,kk) = peakIF;
				 Zfreqkurtosis  (ii,jj,kk) = kurtIF;
				 Zkuraparam		(ii,jj,kk) = mean(abs(K.hilbert.order_parameter));
				 Zstdhilb		(ii,jj,kk) = mean(circ_std(K.hilbert.hilbert(:,:)));

			 end
			
	   	  end
	   	  kk = 0;
		 end
		 jj = 0;
	end
end





if plotsummaries

	figure
	subplot(1,2,1)
	surf(Zfreq(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq (Hz)')
	title('1pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zfreq(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq (Hz)')
	title('1.2pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zfreqkurtosis(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq KURTOSIS (Hz)')
	title('1pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zfreqkurtosis(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq KURTOSIS (Hz)')
	title('1.2pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zpropfiring(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('prop firing neurons')
	title('1pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zpropfiring(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('prop firing neurons')
	title('1.2')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zpopfreq(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('average spike frequency (Hz)')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	title('1pA')
	subplot(1,2,2)
	surf(Zpopfreq(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('average spike frequency (Hz)')
	title('1.2pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})


	figure
	subplot(1,2,1)
	surf(Zkuraparam(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('synchrony (kuramoto parameter)')
	title('1pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zkuraparam(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('synchrony (kuramoto parameter)')
	title('1.2pA')
	set(gca,'xtick', [1:length(noisecorr)], 'ytick', [1:length(gaps)], 'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})




end





if computeselectedxcorr
% [12 16 19 24]
	
	XnoAC  = [];
	lag = 300;
	centerwin = 25;
	noneur = 30;
	xcorrwin = 1:50000;
	nwins = 5; windur = 10000;



	c = 0;
	for simcount = sims2p
		c = c+1;
		xcorr_summa(transients{simcount})

	end


end


	





% legend({'.25 .1' ; '.25 .4' ; '0 .1' ; '0 .4'})
% legend(num2str([12 16 19 24]'))



	% Ptable =

	%     1.0000    0.0500         0    1.0000    0.0155
	%     2.0000    0.0500         0    1.2000    0.1165
	%     3.0000    0.0500    0.1000    1.0000    0.0114
	%     4.0000    0.0500    0.1000    1.2000    0.0659
	%     5.0000    0.0500    0.2000    1.0000    0.0234
	%     6.0000    0.0500    0.2000    1.2000    0.0885
	%     7.0000    0.0500    0.4000    1.0000    0.1510
	%     8.0000    0.0500    0.4000    1.2000    0.3759
	%     9.0000    0.0250         0    1.0000    0.4264
	%    10.0000    0.0250         0    1.2000    1.1095
	%    11.0000    0.0250    0.1000    1.0000    0.2401
	%    12.0000    0.0250    0.1000    1.2000    0.7203
	%    13.0000    0.0250    0.2000    1.0000    0.1792
	%%%    14.0000    0.0250    0.2000    1.2000    0.5219
	%    15.0000    0.0250    0.4000    1.0000    0.2953
	%    16.0000    0.0250    0.4000    1.2000    0.6309
	%    17.0000    0.0010         0    1.0000    2.3214
	%    18.0000    0.0010         0    1.2000    3.0339
	%    19.0000    0.0010    0.1000    1.0000    1.9656
	%    20.0000    0.0010    0.1000    1.2000    2.6470
	%    21.0000    0.0010    0.2000    1.0000    1.6602
	%%%    22.0000    0.0010    0.2000    1.2000    2.2909
	%    23.0000    0.0010    0.4000    1.0000    1.2583
	%    24.0000    0.0010    0.4000    1.2000    1.8228
