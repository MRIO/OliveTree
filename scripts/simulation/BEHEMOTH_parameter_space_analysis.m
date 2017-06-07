% parameter_space_analysis.m
prep_conditions    = 1;
compute_transients = 1;
volumetric_activity = 0;

savemovies = 0; flush = 1; %flush saves the results

% global 
dt =0.025; %ms
fs = 1/dt; %kHz

% reset noise seed
rng(1, 'twister')


% fieldstosave =  {'V_soma', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r', ...
 % 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a' , 'g_GABA_soma','g_GABA_dend', 'g_AMPA', 'Ca2_soma','curr_noise'};

% parameters
fieldstosave = {'V_soma'};%, 'V_dend'};
cell_function = 'vanilla';
steadystatetime = 5; transienttime = 10000;
gapset    = [0 0.005 0.01 0.05]; % gapset = [0 0.001 0.01 0.1];
connset   = [3];   
radiuses  = [2.5 10];
noneurons = 1000;
alphap    = [7];    % alpha parameter
betp      = [3];    % beta parameter
ou_theta  = [1 3];  % ou noise generator parameter (theta)
ou_sigma  = [2 5 10];   % ou noise generator parameter (sigma)
ou_mu	  = [0];    % ou noise generator parameter (mu)
sametoall = [0 .5 1];

betasup  = [.5 1.2]% support for the beta distribution generating low threshold calcium conductance values.

gh_dist = .12 + rand(noneurons,1)*  .01 ; % adds variability to the ih distribution
gint_dist = .11 + rand(noneurons,1)*.02 ;  % adds variability to the internal conductance
gl_dist = .013 + rand(noneurons,1)* .003;  % adds variability to the leak conductances


% g_l  	  = 0.012 + rand(noneurons,1)*.008; % apply larger leak to restitute excitability when gaps are present
% GCAL = BetaDistributions('alpha', [5 6 7], 'beta', [1 2 3], 'no_draws', 1000, 'support', [0 1],'plot_distributions',1);

%                                     __          __                        _            __      
%   _________  ____ ___  ____  __  __/ /____     / /__________ _____  _____(_)__  ____  / /______
%  / ___/ __ \/ __ `__ \/ __ \/ / / / __/ _ \   / __/ ___/ __ `/ __ \/ ___/ / _ \/ __ \/ __/ ___/
% / /__/ /_/ / / / / / / /_/ / /_/ / /_/  __/  / /_/ /  / /_/ / / / (__  ) /  __/ / / / /_(__  ) 
% \___/\____/_/ /_/ /_/ .___/\__,_/\__/\___/   \__/_/   \__,_/_/ /_/____/_/\___/_/ /_/\__/____/  
%                    /_/                                                                         



fff = figure% ,set(fff,'visible', 'off');
set(fff,'visible', 'off')
sim = 0;



for randomseed = [1]
	rng(randomseed, 'twister')

for sta = sametoall
for bet = betp
	def_neurons = createDefaultNeurons(noneurons);
    GCAL = BetaDistributions('alpha', alphap, 'beta', bet, 'no_draws', ...
       noneurons, 'support', betasup,'plot_distributions',0);
       g_CaL = GCAL.sampleDraws{1}';
       distribution_parameters = GCAL.parameters{1};
       def_neurons.g_CaL = g_CaL;
       def_neurons.g_h = gh_dist;
       def_neurons.g_int = gint_dist;

		for gaps = gapset
			for conntype = connset
				for radius = radiuses
                switch conntype
                    case 1
                        netsize = [1 25 40];
                        W_3d_trans1 = createW('3d',netsize, radius ,gaps,1,0);
                        W_3d_trans = W_3d_trans1;

                    case 2
                        netsize = [5 10 20]; 
                        W_3d_trans2 = createW('3d_euclidean_rndwalk',netsize, radius, gaps,1,0);
                        W_3d_trans = W_3d_trans2; 

                    case 3
                        netsize = [5 10 20];
                        W_3d_trans3 = createW('3d_euclidean_rndwalk',netsize, radius, gaps,1,0);
                        W_3d_trans = W_3d_trans3;

                    case 4
                        netsize = [5 10 20];
                        W_3d_trans4 = createW('3d',netsize, radius   ,gaps,1,0);
                        W_3d_trans = W_3d_trans4;
	                
	                case 5
                        netsize = [5 4 5];
                        W_3d_trans4 = createW('3d',netsize, radius   ,gaps,1,0);
                        W_3d_trans = W_3d_trans4;

                end
                rng(randomseed, 'twister')
                [st_st] = IOnet('cell_function', cell_function,'networksize', netsize,'time',steadystatetime,'delta',0.05,'cell_parameters',def_neurons,'W',W_3d_trans.W,'to_report', fieldstosave);

                for theta = ou_theta
				% connectivity

					for sigma = ou_sigma
					
			
						sim = sim+1;

                        rng(randomseed, 'twister')
					    [transients] = IOnet('cell_function', cell_function, 'networksize', netsize,'time',transienttime,'delta',dt,'cell_parameters',def_neurons ,'W',W_3d_trans.W,'ou_noise', [theta sigma ou_mu 0],'to_report', fieldstosave,'tempState', st_st.lastState,'sametoall',sta);%,'Potassium_n', 'Potassium_x_s', 'Ca2_soma', 'Ca2Plus', 'V_dend','curr_noise'});


					    results{sim} = replayResults(transients, [transienttime-1:transienttime], 0,0,fff);
					    transients.results = results{sim};
					    
					    if 1
						    figure(fff)
						    export_fig(num2str(sim),fff)
						    clf
						else
							figure(fff)
							clf
						end
					    
					    if not(transients.failed)
						    orderparameter = measureGlobalSync(transients,[1:transienttime],0);
						else
							orderparameter =0;
						end
				

					    RESULTS(sim,1 ) = gaps;
					    RESULTS(sim,2 ) = sigma;
					    RESULTS(sim,3 ) = theta;
					    RESULTS(sim,4 ) = bet;
					    RESULTS(sim,5 ) = radius;
					    RESULTS(sim,6 ) = results{sim}.propspkneurons;
					    RESULTS(sim,7 ) = mean(sum(results{sim}.spikespercell)/(transienttime*prod(netsize)))*1e3;
					    RESULTS(sim,8 ) = mean(results{sim}.medfreq(results{sim}.medfreq>0));
					    RESULTS(sim,9 ) = mean(mean(transients.networkHistory.V_soma,2));
					    RESULTS(sim,10) = sum(W_3d_trans.W(W_3d_trans.W~=0))/noneurons;
						RESULTS(sim,11) = median(W_3d_trans.stats.connections);
						RESULTS(sim,12) = min(W_3d_trans.stats.connections);
						RESULTS(sim,13) = max(W_3d_trans.stats.connections);
						RESULTS(sim,14) = orderparameter.stats.firstordersync(end);
						RESULTS(sim,15) = orderparameter.stats.secondordersync(end);
						RESULTS(sim,16) = orderparameter.stats.overallsync(end);
						RESULTS(sim,17) = mean(W_3d_trans.stats.clustercoeff.bu);
					    RESULTS(sim,18) = transients.failed;
					    RESULTS(sim,19 ) = randomseed;

					    overallspkspercell{sim} = results{sim}.spikespercell;
					   
					   if flush
							save(['transients_' num2str(sim)],  'transients')
							clear transients
							gpuDevice(1)
							save RESULTS RESULTS
						end
					end
				end
			end
		end
	end
end
end
end


summarize = 0;
if summarize

	%%
	figure
	scatter3(RESULTS(:,1), RESULTS(:,3), RESULTS(:,2),100*RESULTS(:,5)+eps, RESULTS(:,6)*100+eps,'filled')

	clear Z
	dim1 = RESULTS(:,1); 
	dim2 = RESULTS(:,2);
	dim3 = RESULTS(:,3);

	vals1 = unique(RESULTS(:,1));
	% vals2 = [5.3 7.3 10.3 12.3 14.3] % unique(RESULTS(:,2));
	vals2 = [5.6 7.6 10.6 12.6 14.6];
	vals3 = unique(RESULTS(:,3));

	res1  = RESULTS(:,6);

	clamp = 4;

	for i = 1:length(vals1)
		for j = 1:length(vals2)
			
				Z(i,j) = res1(find(dim1==vals1(i) & dim2==vals2(j) & dim3 == clamp ));
			
		end
	end


	% idn = find(RESULTS(:,2) == 7.3 );
	% mesh(RESULTS(idn,1), RESULTS(idn,3), RESULTS(idn,2),100*RESULTS(idn,6)+eps)

end

