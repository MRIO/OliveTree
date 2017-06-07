% overlapping_projections.m

%sim parameters 
dt = 0.025;
simtime = 5000;


% create network
netsz = [3 20 20];
rd = 2.5;
meannoconn = 10;
noneurons = prod(netsz);
gaps = [0.05 0.01 0];
plotconn = 0;

W = createW('3d_euclidean_rndwalk', netsz, rd, 1, 1, plotconn,  3,meannoconn)

% create neurons
def_neurons = createDefaultNeurons(noneurons);
	g_CaL = linspace(.5,1.1,noneurons);
	g_CaL = g_CaL(randperm(noneurons));
	def_neurons.g_CaL = g_CaL;

% noise levels
sametoall = 0.1;
noise_level = [3 3 0 0];

% create perturbations
onsets1 = [ -50:20:50 ];
onsets2 = [-30:10:30];
% onsets2 = [-50:20:50];

connrad  = 5;

	% scatter3(W.coords(:,1), W.coords(:,2), W.coords(:,3), 100, M1+M2,'filled'), axis equal
	% cm = [133 131 115; 150 147 130 ; 250 247 221; 255 255 255]/255;
	% colormap(cm)

pert.mask{1}  	   = M1;
pert.amplitude{1} = 2;
pert.type{1}	  = 'ampa';
pert.duration{1}  = 1;

pert.mask{2}  	   = M2;
pert.amplitude{2} = 2;
pert.type{2}	  = 'gaba_soma';
pert.duration{2}  = 5;


transienttime = 2000;
if ~exist('transients')
	noise_level_transients = [3 5 0 0];
	[transients] = IOnet( 'networksize', netsz ,'time',transienttime,'delta',dt,'cell_parameters', def_neurons ,'W',W.W*gaps(1),'ou_noise', noise_level_transients, 'sametoall',sametoall);
	[continuation] = IOnet( 'networksize', netsz ,'time',500,'delta',dt,'cell_parameters', def_neurons ,'W',W.W*gaps(1),'ou_noise', noise_level_transients, 'sametoall',sametoall', 'tempState',transients.lastState);

	[v pks] = findpeaks(mean(continuation.networkHistory.V_soma),'minpeakdistance',80);
	pko = pks(2);
	save sim_overlap_proj_transients transients continuation
end

s = 0;
for g = gaps
	for on1 = onsets1
		for on2 = onsets2

			s = s+1;
			M1 = create_input_mask(netsz, 'dist_to_point', connrad, W.coords,[1.5 8 10]);
			M2 = create_input_mask(netsz, 'dist_to_point', connrad, W.coords,[1.5 13 10]);


			pert.triggers{1}  = pko + on1 + [1:5];
			pert.triggers{2}  = pko + on1 + on2;

			sim3D = IOnet('networksize', netsz,'time',simtime,'delta',dt,'cell_parameters',def_neurons,'tempState',transients.lastState,'W',W.W*g ,'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall);
			
			eval(['save sim_overlap_proj' num2str(s) ' sim3D'])
		end
	end
end


	
% for ss = [1:s]
% 	animate_volume(sim3D{s},[],1,1)
% 	eval(['!mv volume.mp4 ' num2str(ss) '.mp4'])
% end

