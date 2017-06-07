% simulate_3d_layer_perturbation.m

compute_transients = 1;
display_results    = 1;

savemovies = 1;

rng(1, 'twister')

% defaults
dt =0.01;
fs = 1/dt;
depth = 5; breadth = 5;  height = 5;
netsize= [depth breadth height];
noneurons = prod(netsize);



def_neurons = createDefaultNeurons(noneurons);
	alpha  = [7]; %edges for alpha parameter
	beta   = [4];       %edges for beta parameter
	sup    = [.5 1.2];

	GCAL = BetaDistributions('alpha', alpha, 'beta', beta, 'no_draws',noneurons, 'support', sup); 
		g_CaL = GCAL.sampleDraws{1}';
		distribution_parameters = GCAL.parameters{1};
		def_neurons.g_CaL = g_CaL;

rad = 3;
W_3d = createW('3d_euclidean',netsize, rad, 1,1,1,1);
		W = W_3d.W;
		
ou_noise = [7 3 0 0];
	sametoall = [1 :-.2: 0];

gaps = [0.01 .001 0];

desync_time = 2000;


if ~exist('noiselesstransient')
	[noiselesstransient] = IOnet('networksize', netsize,'time',1000,'delta',0.05,'cell_parameters',def_neurons ,'W',W*1);
end
	
for sta  = sametoall;
	rng(1, 'twister')

	for g = gaps

		[desync_transients] = IOnet('networksize', netsize, 'time',desync_time,'delta',dt,...
			'cell_parameters', def_neurons, 'W',W*g, 'tempState',noiselesstransient.lastState, 'ou_noise', ou_noise);


		flush = 1;
		if flush
			save(['desync_transients_' num2str(sta) '_' num2str(g)],  'desync_transients')
			clear transients
		end


	end


end
	




% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 




