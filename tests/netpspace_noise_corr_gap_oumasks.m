% noise_corr_gap_comparisons.m

% script parameters
compsim = 0;
summarize = 0;
plotsummaries = 0;
plottraces = 0;
plotfigs = 0;
plotcorrelations = 1;



% clear
rndseed = 0;
rng(rndseed,'twister')

% [=================================================================]
%  parameter space grid
% [=================================================================]

gaps = [0.05 0];
noisecorr = [0 .25 .5] ; % proportion of shared noise
noiseamps = [.5 1]; %uA/cm^2 -> x .1 nA for a cell with 10000um^2

disp(['updated!'])
disp(['using:' num2str(gaps)])
disp(['using:' num2str(noisecorr)])
disp(['using:' num2str(noiseamps)])

gpu = 1;

% [=================================================================]
%  readouts
% [=================================================================]

activations =  {'V_soma','V_dend','V_axon','Calcium_l', 'Calcium_r', 'Ca2Plus', 'Potassium_s', 'Hcurrent_q', 'Hcurrent_q','Sodium_m_a', 'Sodium_h_a','Potassium_x_a'};
currents = {'V_soma','V_dend','V_axon', 'I_CaL', 'I_ds', 'I_as', 'I_Na_s', 'I_ls', 'I_Kdr_s', 'I_K_s', 'I_CaH', 'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_h_s', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'};
vsoma = {'V_soma'};

to_report = vsoma;

% [=================================================================]
%  simulation parameters
% [=================================================================]

cell_function = 'vanilla'; % 'devel'
steadystate_time = 100;
simtime  = 1000;

delta = .025;
Fs = 1000; % sampling rate =! delta

% [=================================================================]
%  OU noise parameters
% [===========================================================q======]

tau = 5;
theta = 1/tau;

%% cell parameters

netsize = [4  3  10];
noneurons = prod(netsize); 
% netsize = [1 noneurons 1]; 
def_neurons = createDefaultNeurons(noneurons,'celltype', 'randomized','gapcompensation', 0);
Plist = def_neurons.Plist;


% [================================================]
% 		gap connections
% [================================================]
 W = zeros(noneurons);
 % W = createW(noneurons);

noconnections = 8;
randomize = 1;
symmetrize = 1;
radius = 3;
iter = 3;
scaling = 1;
plotconns = 0;
normneighbors = 1;
% W = createW('all to all', netsize, [],  1,  randomize,         0,        0,       noconnections, [],symmetrize);
W = createW('3d_euclidean_rndwalk', [4 3 10], radius, scaling, randomize, plotconns, iter , noconnections,[],symmetrize,[0 0 0 0],normneighbors)

% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)


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
   	 	rng(0,'twister')
   	 	gnoise = [theta namp 0 5];
   	 	simcount = simcount+1;




   	 	% [================================================]
		%  Ornstein Uhlenbeck Perturbation over masks
		% [================================================]
	


		% create overlapping ou-masks
		pert = [];

		numberofmasks = 3;

		onset_of_stim = 100;
		stim_dur      = 800;
		offset_stim   = 50;
		synapseprobability = .5;

		th =	 1/10 ; % decay time parameter
		sig = 	 .2 ; % pA
		mu = 	 0  ; % pA
		mix =    .8 ; % amount of noise shared between neurons in the mask 

			for nm = 1:numberofmasks

				pert.mask  	  {nm} = rand(noneurons,1)>synapseprobability;
				% pert.mask  	  {nm}(nm) = 1;

				% pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', synapseprobability);
				pert.amplitude{nm} = 1;
				pert.triggers {nm} = onset_of_stim + (nm-1)*offset_stim;
				pert.duration {nm} = stim_dur;
				pert.type	  {nm} = 'ou_noise';

				pert.param{nm}(1)  = th  ;
				pert.param{nm}(2)  = mu  ;
				pert.param{nm}(3)  = sig ;
				pert.param{nm}(4)  = mix;


			end 



		   transients{simcount} = IOnet('tempState', st_st.lastState ,'cell_parameters', def_neurons, ...
		   	'networksize', netsize,'appCurrent',I_app,'time',simtime ,'W', W.W*gs ,'ou_noise', gnoise , ...
		   	'to_report', to_report ,'gpu', gpu ,  ...
		   	'cell_function', cell_function ,'delta',delta,'sametoall', ncorr);



		   transients{simcount}.scriptParameters.noiseAmplitude   = namp;
		   transients{simcount}.scriptParameters.noiseCorrelation = ncorr;
		   transients{simcount}.scriptParameters.gapAmplitudes    = gs;
		   transients{simcount}.spikes = spikedetect(transients{simcount}, 0,0);
		   transients{simcount}.Plist = Plist;

		    if plotfigs
				figure
				imagesc(transients{simcount}.networkHistory.V_soma,[-80 -20]), colorbar
				set(gca,'ytick', [1:noneurons],'yticklabel', num2str(Plist),'fontsize',8)
				title([num2str(transients{simcount}.spikes.popfrequency) ' Hz'])

				figure
				ca = axis;
				set(0,'defaultaxescolororder', linspecer(length(Plist)));
				p = plot([1:simtime],   transients{simcount}.networkHistory.V_soma');
				legend(num2str(Plist))
			end


   	  end
	 end
	end

	save noisecorr_gap_pspace -v7.3

end
   
% [===========================================================================================================]




% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Hcurrent_q'),legend(num2str(Plist)), title('V vs q (Hcurrent)')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Sodium_h_a'),legend(num2str(Plist)), title('V vs Sodium\_h axon')
% figure, plot(transients.networkHistory.V_soma',transients.networkHistory.Calcium_r'),legend(num2str(Plist)),title('V vs Calcium\_r')




% [=================================================================]
%  plot trace matrix
% [=================================================================]

if plottraces
	ind = 0;
	tslice = [1:simtime];
	ind = 0; pind = 1;
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
						% X{ind} = xcorr(T{ind},lags,'coeff');
						X{ind} = xcorr(T{ind},lags,'unbiased');
					
				end
			end
		end
	end


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

calcfreq = 0;
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
				 K = measureGlobalSync(transients{simcount});
				 
			 	 InstFreq = diff(unwrap(K.hilbert.hilbert)') * Fs/(2*pi);
			 	 [HistInstFreq x] = hist(InstFreq,201); % per cell
			 	 [HistInstFreq_overall x] = hist(InstFreq(:),x); %InstFreq
			 	 kuraparam(ii,jj,kk) = mean(abs(K.hilbert.order_parameter));
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
			end

			
	   	  end
	   	  kk = 0;
		 end
		 jj = 0
		end
end



if plotsummaries

	figure
	subplot(1,2,1)
	surf(Zfreq(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq (Hz)')
	title('5pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zfreq(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq (Hz)')
	title('10pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zfreqkurtosis(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq KURTOSIS (Hz)')
	title('5pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zfreqkurtosis(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('STO freq KURTOSIS (Hz)')
	title('10pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zpropfiring(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('prop firing neurons')
	title('5pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(Zpropfiring(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('prop firing neurons')
	title('10pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})

	figure
	subplot(1,2,1)
	surf(Zpopfreq(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('average spike frequency (Hz)')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	title('5pA')
	subplot(1,2,2)
	surf(Zpopfreq(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('average spike frequency (Hz)')
	title('10pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})


	figure
	subplot(1,2,1)
	surf(kuraparam(:,:,1))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('synchrony (kuramoto parameter)')
	title('5pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})
	subplot(1,2,2)
	surf(kuraparam(:,:,2))
	ylabel('gap conductance (mS/cm^2)')
	xlabel('noise correlations')
	zlabel('synchrony (kuramoto parameter)')
	title('10pA')
	set(gca,'yticklabel',{num2str(gaps')}, 'xticklabel', {num2str(noisecorr')})




end



parametertable = 0;
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


			end
		end
	end
end