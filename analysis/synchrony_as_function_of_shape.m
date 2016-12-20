% synchrony_as_function_of_shape.m
rng(0,'twister')

% connectivity_comparisons.m
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions)
load('/Users/M/Public/Dropbox/WholeDamnOlive/DensityAnalysis/JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;


% examples of connectivity
 rd = 2.5;
 ni = 10;
 gap = .01;
 plotfig =0;

W1= createW('3d_euclidean_rndwalk', [5 10 20], rd, gap, 1, plotfig, 1, ni); set(gcf, 'color', [1 1 1]);
W2= createW('3d_euclidean',  		[25 20 1], rd, gap, 1, 0,1,ni); set(gcf, 'color', [1 1 1]);
W3= createW('3d_euclidean',  		[5 10 20], rd, gap, 1, 0,1,ni); set(gcf, 'color', [1 1 1]);
W4= createW('3d_euclidean', 		[10 10 10], rd, gap, 1, 0,1,ni); set(gcf, 'color', [1 1 1]);
W5= createW('3d_reconstruction', 	[], rd*15, gap, 1, 0, [], ni, somatapositions); set(gcf, 'color', [1 1 1]);

netsize1 = [5 10 20];
netsize2 = [25 20 1];
netsize3 = [5 10 20];
netsize4 = [10 10 10];
netsize5 = [1 length(somatapositions) 1];

transienttime = 3000; dt = .05; noise_parameters = [3 1 0 0]; sametoall = .0;
if 0
	for o = [1:5]

		W_ = eval(['W' num2str(o)]);
		W = W_.W;
		netsize = eval(['netsize' num2str(o)]);

		noneu = prod(netsize);

		defneu = createDefaultNeurons(noneu,1);
		defneu.g_CaL = .7  + rand(noneu,1)*.35;
		defneu.g_int = .12 + rand(noneu,1)*.1;

		[steady_state{o}] = IOnet( 'networksize', netsize ,'time',transienttime,'delta',dt,'cell_parameters', defneu ,'W', W,'ou_noise', noise_parameters, 'sametoall',sametoall);

		save sync_func_shape
	end
end



transienttime = 3000; dt = .025; noise_parameters = [3 1 0 0]; sametoall = .1;
rng(0,'twister')
if 1
	for o = [1:5]

		W_ = eval(['W' num2str(o)]);
		W = W_.W;
		netsize = eval(['netsize' num2str(o)]);

		noneu = prod(netsize);

		defneu = createDefaultNeurons(noneu,1);
		defneu.g_CaL = .7  + rand(noneu,1)*.35;
		defneu.g_int = .12 + rand(noneu,1)*.1;

		[steady_state{o}] = IOnet('tempState', steady_state{o}.lastState, 'networksize', netsize ,'time',transienttime,'delta',dt,'cell_parameters', defneu ,'W', W,'ou_noise', noise_parameters, 'sametoall',sametoall);

		save sync_func_shape_noise
	end
end



if 0

	for o = [1:5]

		W_ = eval(['W' num2str(o)]);
		W = W_.W;
		netsize = eval(['netsize' num2str(o)]);

		noneu = prod(netsize);

		% replayResults(steady_state{o},[],0)
		ss_sync{o} = measureGlobalSync(steady_state{o},[2000:2800],1);
		title(num2str(o))

	end
end


if 1 
	for o = [1:5]
		
		msync_all(o) = ss_sync{o}.stats.overallsync(1);
		msync_fo(o) = ss_sync{o}.stats.firstordersync(1);
	end
end






