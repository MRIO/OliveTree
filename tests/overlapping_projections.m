% overlapping_projections.m

%script parameters
prep     = false;
simulate = false;
analyze  = true;


%sim parameters 
dt = 0.025;
simtime = 1200;


% create network
netsz = [3 20 20];
rd = 2;
meannoconn = 10;
noneurons = prod(netsz);
gaps = [0 0.05];
plotconn = 0;

if prep
	W = createW('3d_euclidean_rndwalk', netsz, rd, 1, 1, plotconn,  3,meannoconn);

	% create neurons
	def_neurons = createDefaultNeurons(noneurons);
		g_CaL = linspace(.5,1.1,noneurons);
		g_CaL = g_CaL(randperm(noneurons))';
		def_neurons.g_CaL = g_CaL;

	% noise levels
	sametoall = 0.1;
	noise_level = [3 10 0 0];

	% create perturbations
	onsets1 = [ 120 ];
	onsets2 = [-50:10:50 ];
	% onsets2 = [-50:20:50];

	connrad  = 5;
	M1 = create_input_mask(netsz, 'dist_to_point', connrad, W.coords,[1.5 8 10]);
	M2 = create_input_mask(netsz, 'dist_to_point', connrad, W.coords,[1.5 13 10]);

		% plot stimulation pattern, uncomment
		% 
		% scatter3(W.coords(:,1), W.coords(:,2), W.coords(:,3), 100, M1+M2,'filled'), axis equal
		% cm = [133 131 115; 150 147 130 ; 250 247 221; 255 255 255]/255;
		% colormap(cm)

	pert.mask{1}  	   = M1;
	pert.amplitude{1} = 2;
	pert.type{1}	  = 'gaba_soma';
	pert.duration{1}  = 5;

	pert.mask{2}  	   = M1;
	pert.amplitude{2} = 2;
	pert.type{2}	  = 'gaba_dend';
	pert.duration{2}  = 5;


	pert.mask{3}  	   = M2;
	pert.amplitude{3} = 2;
	pert.type{3}	  = 'ampa';
	pert.duration{3}  = 1;
end

if simulate
	transienttime = 1000;
	if ~exist('transients')
		noise_level_transients = [3 10 0 0];
		dt = 0.025;
		[transients] = IOnet( 'networksize', netsz ,'time',transienttime,'delta',dt,'cell_parameters', def_neurons ,'W',W.W*gaps(1),'ou_noise', noise_level_transients, 'sametoall',sametoall);
	end

	s = 0;
	for g = gaps
		for on1 = onsets1
			for on2 = onsets2
				s = s+1;

				pert.triggers{1}  = on1 + [1:3];
				pert.triggers{2}  = on1 + [1:3];
				pert.triggers{3}  = on1 + on2 + [1:3];

				sim3D{s} = IOnet('networksize', netsz,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'tempState',transients.lastState,'W',W.W*g ,'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall);
			end
		end
	end

end
		
% for ss = [1:s]
% 	animate_volume(sim3D{s},[],1,1)
% 	eval(['!mv volume.mp4 ' num2str(ss) '.mp4'])
% end


% analysis
if analyze
	files  =1:100
	for f = files
		eval(['load /Users/M/Public/BitSync/Bench/overlapping_projs/sim_overlap_proj' num2str(f) '.mat']);
		allneu = 1:noneurons;
		VS = sim3D.networkHistory.V_soma;
		H = hilbert_of_membranepotential(VS); 
		U = H.hilbert;

		pert_mask1 = sim3D.perturbation.mask{1};
		pert_mask2 = sim3D.perturbation.mask{2};

		pert_time1 = sim3D.perturbation.triggers{1}(1)
		pert_time2 = sim3D.perturbation.triggers{2}(1)

		% stimulation masks
		group1 = (pert_mask1==1);
		group2 = (pert_mask2==1);

		group12 = pert_mask1 & pert_mask2;
		group1 = group1 & ~group12;
		group2 = group2 & ~group12;
		other = ~group12;

		% group by stimulation mask
		G1 = U(group1,:);
		G2 = U(group2,:);
		G12 = U(group12,:);
		GO =  U(other,:);

		order_parameter_G1 = mean( exp(i*(G1+pi)));
		order_parameter_G2 = mean( exp(i*(G2+pi)));
		order_parameter_G12 = mean( exp(i*(G12+pi)));
		order_parameter_GO = mean( exp(i*(GO+pi)));

		figure(1)
		clf
		plot(abs(order_parameter_G1),'r')
		hold on
		plot(abs(order_parameter_G2),'b')
		plot(abs(order_parameter_G12),'color', [1 1 0])
		plot(abs(order_parameter_GO),'k')
		ylim([0.5 1])
		xlim([0 500])

		export_fig(['OverlapplingProjections_kparam' num2str(pert_time1) '-' num2str(pert_time2) '_' num2str(f) '.png'])

		figure(2)
		clf
		plot(imag(order_parameter_G1),'r')
		hold on
		plot(imag(order_parameter_G2),'b')
		plot(imag(order_parameter_G12),'color', [1 0 1])
		plot(imag(order_parameter_GO),'k')
		
		xlim([0 500])

		export_fig(['OverlapplingProjections_phase' num2str(pert_time1) '-' num2str(pert_time2) '_' num2str(f) '.png'])

		% disp('paused')
		
		
		result{f}.order_parameter_G1 = order_parameter_G1;
		result{f}.order_parameter_G2 = order_parameter_G2;
		result{f}.order_parameter_G12 = order_parameter_G12;
		result{f}.order_parameter_GO = order_parameter_GO;

		[mean(abs(order_parameter_G1)) mean(abs(order_parameter_G2)) mean(abs(order_parameter_G12)) mean(abs(order_parameter_GO))]
		R(f,:) = [mean(abs(order_parameter_G1)) mean(abs(order_parameter_G2)) mean(abs(order_parameter_G12)) mean(abs(order_parameter_GO))];



	end





end

