% connectivity_comparisons.m
% (Note: for the real data we scale connectivity radius by median(mindist)*
% * to include only the closest neighbors

% [=================================================================]
%  models to try
% [=================================================================]
% models = {'reconstruction' '3dbrick'  '3drandom' 'random' '2d_euclidean' 'reconstruction_clusterized'};
% models = {'reconstruction' '3dbrick' '3drandom' '2d_chebychev' '2d_euclidean' 'random'};
% models = {'reconstruction'  '2d_euclidean' '2d_chebychev' '3dbrick' '3dcube' 'random'};
models = {'reconstruction' };



% [=================================================================]
%  parameters
% [=================================================================]
n_connections = [20];
radiuses = [5];

% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions)



% examples of connectivity

%  rd = 2.5;
%  ni = 10;
% createW('3d_euclidean_rndwalk', [5 10 20], rd, 1, 1, plotfig, 1, ni); set(gcf, 'color', [1 1 1])
% createW('3d_euclidean', [25 20 1], rd, 1, 1, 1,1,ni); set(gcf, 'color', [1 1 1])
% createW('3d_euclidean', [5 10 20], rd, 1, 1, 1,1,ni); set(gcf, 'color', [1 1 1])
% createW('3d_euclidean', [10 10 10], rd, 1, 1, 1,1,ni); set(gcf, 'color', [1 1 1])
% createW('3d_reconstruction', [], rd*40, 1, 1, 1, [], ni, somatapositions); set(gcf, 'color', [1 1 1])

load('JM394_horizontal_coordinates-MAO.mat')
somatapositions = JM394_horizontal_coordinates;
somatapositions(1,:) = [];
out = createW('3d_reconstruction', [], 5*40, 1, .01, blaplot, [1 20 1 0], 20, somatapositions,1,[1 20 1 0]);

 plotfig = 0;



 % [=================================================================]
 %  load reconstruction data
 % [=================================================================]

addpath('/Users/M/Projects/Experiments/Olive/Experiments/Anatomy/Somata/')
% addpath('/Users/M/Synced/Projects/Experiments/Olive/Experiments/Anatomy/Somata/')

% MAO slice
reconstruction_data = 'MAO';
% reconstruction_data = 'all';

switch reconstruction_data
	case 'MAO'
		load('JM394_horizontal_coordinates-MAO.mat')
		somatapositions = JM394_horizontal_coordinates;
		somatapositions(1,:) = [];
	case 'all'
		% whole half olive
		load('PDX-AI9-AllCells_merged_5umMinDist.mat')
		somatapositions = R(:,[3 4 8]);
		[nan_i nan_j] = find(isnan(somatapositions));
		somatapositions(nan_i,:) = [];
end

% [=================================================================]
%  cleanup and organize somata positions
% [=================================================================]

span = max(somatapositions)-min(somatapositions);
numcells = length(somatapositions);

% prune somata closer than X um

mindist = 7;
% empirical_minimum_distance
Dmat = squareform( pdist(somatapositions, 'euclidean')  );
Dmat(find(isnan(Dmat)|Dmat==0)) = Inf;
[i j] = find(Dmat<mindist);
somatapositions(i,:) = [];
Dmat = squareform( pdist(somatapositions, 'euclidean')  );

Dmat(find(isnan(Dmat)|Dmat==0)) = Inf;

% [=================================================================]
%  normalize somata positions from reconstruction data
% [=================================================================]

NZTU = triu(ones(size(Dmat))) - eye(size(Dmat));
edges = [1:2:300];
[H] = hist(Dmat(find(NZTU)),edges);
[peakH idx] = max(H);

% empirical_minimum_distance = min(Dmat(find(NZTU)))

D2 = Dmat.*NZTU;
D2(find(tril(ones(size(Dmat)))))= nan;
mD2 = min(D2,[],2); mD2(isnan(mD2))=[];

empirical_minimum_distance = median(mD2);

% somatapositions = somatapositions/empirical_minimum_distance;
span = max(somatapositions)-min(somatapositions);



% [=================================================================]
%  compute connectivity models
% [=================================================================]

 % n_connections = [5:5:40];
 % radiuses = [1.1:5];


if sum(ismember(models, '3drandom'))
	% 3d_euclidean_rndwalk
	i = 0; j = 0; k = 0; n =0; iter = 5;
	for rd = radiuses
		
		for ni = n_connections
			
			for rep = [1]
				
				n = n+1;
				out = createW('all_to_all', [10 10 10], rd, 1, 1, plotfig, iter, ni);
				
				R_3d_rnd(n,1) = rd;
				R_3d_rnd(n,2) = ni;
				R_3d_rnd(n,3) = rep;
				R_3d_rnd(n,4) = mean(out.stats.clustercoeff.bu);
				R_3d_rnd(n,5) = mean(out.stats.clustercoeff.wd);
				R_3d_rnd(n,6) = mean(out.stats.connections);
			end

		end
	end
end


if sum(ismember(models, '3dbrick'))
	% 3d euclidean sheet
	i = 0; j = 0; k = 0; n =0;
	for rd = radiuses
		
		for ni = n_connections;
			
			for rep = [1]
				
				n = n+1;
				out = createW('3d_euclidean', [5 10 20], rd, 1, 1, 0,0,ni);
				
				R_3d_brick(n,1) = rd;
				R_3d_brick(n,2) = ni;
				R_3d_brick(n,3) = rep;
				R_3d_brick(n,4) = mean(out.stats.clustercoeff.bu);
				R_3d_brick(n,5) = mean(out.stats.clustercoeff.wd);
				R_3d_brick(n,6) = mean(out.stats.connections);
			end

		end
	end
end


if sum(ismember(models, '3dcube'))
	% 3d euclidean vol
	i = 0; j = 0; k = 0; n =0;
	for rd = radiuses
		
		for ni = n_connections;
			
			for rep = [1]
				
				n = n+1;
				out = createW('3d_chebychev', [10 10 10], rd, 1, 1, 0,0,ni);
				
				R_3d_cube(n,1) = rd;
				R_3d_cube(n,2) = ni;
				R_3d_cube(n,3) = rep;
				R_3d_cube(n,4) = mean(out.stats.clustercoeff.bu);
				R_3d_cube(n,5) = mean(out.stats.clustercoeff.wd);
				R_3d_cube(n,6) = mean(out.stats.connections);
			end

		end
	end
end


% 2d sheet chebychev
if sum(ismember(models, '2d_chebychev'))
	i = 0; j = 0; k = 0; n =0;
	for rd = radiuses
		
		for ni = n_connections;
			
			for rep = [1]
				
				n = n+1;
				out = createW('3d_chebychev', [25 40 1], rd, 1, 1, 0,0,ni);
				
				R_2d_8(n,1) = rd;
				R_2d_8(n,2) = ni;
				R_2d_8(n,3) = rep;
				R_2d_8(n,4) = mean(out.stats.clustercoeff.bu);
				R_2d_8(n,5) = mean(out.stats.clustercoeff.wd);
				R_2d_8(n,6) = mean(out.stats.connections);
			end

		end
	end
end




if sum(ismember(models, 'reconstruction'))
	% 3d reconstruction
	n = 0;
	for rd = radiuses*empirical_minimum_distance
		
		for ni = n_connections;
			
			for rep = [1]
				n = n+1;
				out = createW('3d_reconstruction', [], rd, 1, 1, 0, [], ni, somatapositions);


				R_reconst(n,1) = rd;
				R_reconst(n,2) = ni;
				R_reconst(n,3) = rep;
				R_reconst(n,4) = mean(out.stats.clustercoeff.bu);
				R_reconst(n,5) = mean(out.stats.clustercoeff.wd);
				R_reconst(n,6) = mean(out.stats.connections);
			end

		end
	end

end






if sum(ismember(models, 'random'))
	
	n = 0;
	for rd = radiuses
		
		for ni = n_connections;
			
			for rep = [1]
				n = n+1;
				out = createW('random', [1000 1 1], rd, 1, 1, 0, [], ni, []);


				R_random(n,1) = rd;
				R_random(n,2) = ni;
				R_random(n,3) = rep;
				R_random(n,4) = mean(out.stats.clustercoeff.bu);
				R_random(n,5) = mean(out.stats.clustercoeff.wd);
				R_random(n,6) = mean(out.stats.connections);
			end

		end
	end

end




if sum(ismember(models, '2d_euclidean'))
	% 3d reconstruction
	n = 0;
	for rd = radiuses
		
		for ni = n_connections;
			
			for rep = [1]
				n = n+1;
				out = createW('3d_euclidean', [25 40 1], rd, 1, 1, 0, [], ni, []);


				R_2d_euclid(n,1) = rd;
				R_2d_euclid(n,2) = ni;
				R_2d_euclid(n,3) = rep;
				R_2d_euclid(n,4) = mean(out.stats.clustercoeff.bu);
				R_2d_euclid(n,5) = mean(out.stats.clustercoeff.wd);
				R_2d_euclid(n,6) = mean(out.stats.connections);
			end

		end
	end

end


blaplot = 1;

if sum(ismember(models, 'reconstruction_clusterized'))
	% 3d reconstruction
	n = 0;
	for rd = radiuses*empirical_minimum_distance
		
		for ni = n_connections;
			
			for rep = [1]
				n = n+1;
				out = createW('3d_reconstruction', [], rd, 1, 1, blaplot, [], ni, somatapositions,1,[1 ni 1 0]);


				R_reconst_clusterized(n,1) = rd;
				R_reconst_clusterized(n,2) = ni;
				R_reconst_clusterized(n,3) = rep;
				R_reconst_clusterized(n,4) = mean(out.stats.clustercoeff.bu);
				R_reconst_clusterized(n,5) = mean(out.stats.clustercoeff.wd);
				R_reconst_clusterized(n,6) = mean(out.stats.connections);
			end

		end
	end

end





	% Gather clustering results _ bu
	i = 0; j = 0; k = 0; n =0;
	for rd = radiuses
		i  = i+1;
		
		for ni = n_connections;
			j = j+1;
			
			for rep = [1]
				n = n+1;
				if sum(ismember(models, 'reconstruction'))
					clusterSURF1_bu(i,j)   = R_reconst (n,4);
					clusterSURF1_wd(i,j)   = R_reconst (n,5);
					connectivitySURF1(i,j) = R_reconst (n,6);
				end
				if sum(ismember(models, '3dbrick'))
					clusterSURF2_bu(i,j)   = R_3d_brick(n,4);
					clusterSURF2_wd(i,j)   = R_3d_brick(n,5);
					connectivitySURF2(i,j) = R_3d_brick(n,6);
				end
				if sum(ismember(models, '3drandom'))
					clusterSURF3_bu(i,j) = R_3d_rnd  (n,4);
					clusterSURF3_wd(i,j) = R_3d_rnd  (n,5);
					connectivitySURF3(i,j) = R_3d_rnd  (n,6);
				end
				if sum(ismember(models, '2d_chebychev'))
					clusterSURF4_bu(i,j) = R_2d_8    (n,4);
					clusterSURF4_wd(i,j) = R_2d_8    (n,5);
					connectivitySURF4(i,j) = R_2d_8    (n,6);
				end
				if sum(ismember(models, '3dcube'))
					clusterSURF5_bu(i,j) = R_3d_cube  (n,4);
					clusterSURF5_wd(i,j) = R_3d_cube  (n,5);
					connectivitySURF5(i,j) = R_3d_cube(n,6);
				end
				if sum(ismember(models, 'random'))
					clusterSURF6_bu(i,j) = R_random  (n,4);
					clusterSURF6_wd(i,j) = R_random  (n,5);
					connectivitySURF6(i,j) = R_random  (n,6);
				end
				if sum(ismember(models, '2d_euclidean'))
					clusterSURF7_bu(i,j)   = R_2d_euclid  (n,4);
					clusterSURF7_wd(i,j)   = R_2d_euclid  (n,5);
					connectivitySURF7(i,j) = R_2d_euclid  (n,6);
				end
				if sum(ismember(models, 'reconstruction_clusterized'))
					clusterSURF8_bu(i,j)   = R_reconst_clusterized (n,4);
					clusterSURF8_wd(i,j)   = R_reconst_clusterized (n,5);
					connectivitySURF8(i,j) = R_reconst_clusterized (n,6);
				end



			end

		end
		j = 0;
	end









figure
	if sum(ismember(models, 'reconstruction'))
	M1_wd = surf( n_connections, radiuses, clusterSURF1_wd), hold on; set(M1_wd, 'edgecolor','g','facecolor','g') ; 
	end
	if sum(ismember(models, '3dbrick'))
	M2_wd = surf( n_connections, radiuses, clusterSURF2_wd), hold on; set(M2_wd, 'edgecolor','b','facecolor','b');  
	end
	if sum(ismember(models, '3drandom'))
	M3_wd = surf( n_connections, radiuses, clusterSURF3_wd), hold on; set(M3_wd, 'edgecolor','r','facecolor','r');  
	end
	if sum(ismember(models,'2d_chebychev'))
	M4_wd = surf( n_connections, radiuses, clusterSURF4_wd), hold on; set(M4_wd, 'edgecolor','c','facecolor','c');  
	end
	if sum(ismember(models,'3dcube'))
	M5_wd = surf( n_connections, radiuses, clusterSURF5_wd), hold on; set(M5_wd, 'edgecolor','b','facecolor','r');  
	end
	if sum(ismember(models,'random'))
	M6_wd = surf( n_connections, radiuses, clusterSURF6_wd), hold on; set(M6_wd, 'edgecolor','k','facecolor','k');  
	end
	if sum(ismember(models,'2d_euclidean'))
	M7_wd = surf( n_connections, radiuses, clusterSURF7_wd), hold on; set(M7_wd, 'edgecolor','y','facecolor','y');  
	end
	if sum(ismember(models,'reconstruction_clusterized'))
	M8_wd = surf( n_connections, radiuses, clusterSURF8_wd), hold on; set(M8_wd, 'edgecolor','y','facecolor','w');  
	end

	camlight left; %
	alpha(.3)
	title('cluster coeff weighted directed')
	ylabel('connection radius')
	xlabel('number of connections')
	zlabel('custer coefficient')
	legend(models)


figure
	if sum(ismember(models, 'reconstruction'))
	M1_bu = surf( n_connections, radiuses, clusterSURF1_bu), hold on; set(M1_bu, 'edgecolor','g','facecolor','g') ; 
	end
	if sum(ismember(models, '3dbrick'))
	M2_bu = surf( n_connections, radiuses, clusterSURF2_bu), hold on; set(M2_bu, 'edgecolor','b','facecolor','b');  
	end
	if sum(ismember(models, '3drandom'))
	M3_bu = surf( n_connections, radiuses, clusterSURF3_bu), hold on; set(M3_bu, 'edgecolor','r','facecolor','r');  
	end
	if sum(ismember(models,'2d_chebychev'))
	M4_bu = surf( n_connections, radiuses, clusterSURF4_bu), hold on; set(M4_bu, 'edgecolor','c','facecolor','c');  
	end
	if sum(ismember(models,'3dcube'))
	M5_bu = surf( n_connections, radiuses, clusterSURF5_bu), hold on; set(M5_bu, 'edgecolor','b','facecolor','r');  
	end
	if sum(ismember(models,'random'))
	M6_bu = surf( n_connections, radiuses, clusterSURF6_bu), hold on; set(M6_bu, 'edgecolor','k','facecolor','k');  
	end
	if sum(ismember(models,'2d_euclidean'))
	M7_bu = surf( n_connections, radiuses, clusterSURF7_bu), hold on; set(M7_bu, 'edgecolor','y','facecolor','y');  
	end
	if sum(ismember(models,'reconstruction_clusterized'))
	M8_bu = surf( n_connections, radiuses, clusterSURF8_bu), hold on; set(M8_bu, 'edgecolor','y','facecolor','y');  
	end
	camlight left; %
	alpha(.3)
	title('cluster coeff binary undirected')
	ylabel('connection radius')
	xlabel('number of connections')
	zlabel('custer coefficient')
	legend(models)


figure
	if sum(ismember(models, 'reconstruction'))
		M1 = surf( n_connections ,radiuses, connectivitySURF1), hold on;  set(M1, 'edgecolor','g','facecolor','g') ; 
	end
	if sum(ismember(models, '3dbrick'))
		 M2 = surf( n_connections ,radiuses, connectivitySURF2), hold on; set(M2, 'edgecolor','b','facecolor','b');  
	end

	if sum(ismember(models, '3drandom'))
		 M3 = surf( n_connections ,radiuses, connectivitySURF3), hold on; set(M3, 'edgecolor','r','facecolor','r');  
	end

	if sum(ismember(models,'2d_chebychev'))
		 M4 = surf( n_connections ,radiuses, connectivitySURF4), hold on; set(M4, 'edgecolor','c','facecolor','c');  
	end

	if sum(ismember(models,'3dcube'))
		 M5 = surf( n_connections ,radiuses, connectivitySURF5), hold on; set(M5, 'edgecolor','b','facecolor','r');  
	end

	if sum(ismember(models,'random'))
		 M6 = surf( n_connections ,radiuses, connectivitySURF6), hold on; set(M6, 'edgecolor','k','facecolor','k');  
	end
	if sum(ismember(models,'2d_euclidean'))
		 M7 = surf( n_connections, radiuses, connectivitySURF7), hold on; set(M7, 'edgecolor','y','facecolor','y');  
	end
	if sum(ismember(models,'reconstruction_clusterized'))
		M8 = surf( n_connections, radiuses, connectivitySURF8), hold on; set(M8, 'edgecolor','y','facecolor','y');  
	end

	camlight left; %_wd
	alpha(.3)
	title('number of connections')
	xlabel('number of connections  (parameter)')
	ylabel('connection radius')
	zlabel('number of actual connections')
	legend(models)











% [================================================]
% 		 only generate with radiuses
% [================================================]

dothis  = 1;
if dothis

	fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


	distfunc = 'euclidean';
	% distfunc = 'chebychev';


	% 3d_euclidean_rndwalk
	i = 0; j = 0; k = 0; n =0;

	for rd = radiuses
			
			n = n+1;
			W1 = createW('3d_euclidean', 		 [5 10 20] , rd, 0, 0, plotfig, []  ,    Inf);
			W2 = createW('3d_chebychev', 		 [5 10 20] , rd, 0, 0, plotfig, []  ,    Inf);
			W3 = createW('3d_euclidean', 		 [25 40 1] , rd, 0, 0, plotfig, []  ,    Inf);
			W4 = createW('3d_chebychev', 		 [25 40 1] , rd, 0, 0, plotfig, []  ,    Inf);
			W5 = createW('3d_reconstruction', 	 [] , rd*empirical_minimum_distance, 0, 0, plotfig, []  ,    Inf, somatapositions);
			connav1(n,:) = sum(squareform(pdist(W1.coords, distfunc))<=rd);
			connav2(n,:) = sum(squareform(pdist(W2.coords, distfunc))<=rd);
			connav3(n,:) = sum(squareform(pdist(W3.coords, distfunc))<=rd);
			connav4(n,:) = sum(squareform(pdist(W4.coords, distfunc))<=rd);
			connav5(n,:) = sum(squareform(pdist(W5.coords, distfunc))<=rd*empirical_minimum_distance);
			[poss_conn1(n,:)] = quantile(connav1(n,:) , [.25 .5 .75]);
			[poss_conn2(n,:)] = quantile(connav2(n,:) , [.25 .5 .75]);
			[poss_conn3(n,:)] = quantile(connav3(n,:) , [.25 .5 .75]);
			[poss_conn4(n,:)] = quantile(connav4(n,:) , [.25 .5 .75]);
			[poss_conn5(n,:)] = quantile(connav5(n,:) , [.25 .5 .75]);
	end



	figure
	
	fill_between_lines(radiuses, poss_conn1(:,1)' , poss_conn1(:,3)','r'), alpha(.2), hold on;
	fill_between_lines(radiuses, poss_conn2(:,1)' , poss_conn2(:,3)','b'), alpha(.2);
	fill_between_lines(radiuses, poss_conn3(:,1)' , poss_conn3(:,3)','g'), alpha(.2);
	fill_between_lines(radiuses, poss_conn4(:,1)' , poss_conn4(:,3)','y'), alpha(.2);
	fill_between_lines(radiuses, poss_conn5(:,1)' , poss_conn5(:,3)','k'), alpha(.2);
;
	plot(radiuses, poss_conn1(:,2),'r'), hold on;
	plot(radiuses, poss_conn2(:,2),'b');
	plot(radiuses, poss_conn3(:,2),'g');
	plot(radiuses, poss_conn4(:,2),'y');
	plot(radiuses, poss_conn5(:,2),'k');



	xlabel('radius'); ylabel('connections within radius')

	legend({'3d euclidean', '3d chebychev', '2d chebychev', '2d euclid', 'reconstr'})


end






	% fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


	% distfunc = 'euclidean';
	% % distfunc = 'chebychev';


	% i = 0; j = 0; k = 0; n =0; iter = 2;
	% for rd = radiuses
			
	% 		n = n+1;
	% 		W5 = createW('3d_reconstruction', 	 []        , rd*empirical_minimum_distance, 0, 0, 0, []  ,    Inf, somatapositions);
	% 		connav5(n,:) = sum(squareform(pdist(W5.coords, distfunc))<=rd*empirical_minimum_distance);
	% 		[poss_conn5(n,:)] = quantile(connav5(n,:) , [.25 .5 .75])
	% end



	% figure
	% fill_between_lines(radiuses, poss_conn5(:,1)' , poss_conn5(:,3)','k'), alpha(.2)
	% plot(radiuses, poss_conn5(:,2),'k')


