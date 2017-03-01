% retrieveNeuronsByClass.m
function selectedneurons = retrieveNeuronsByClass(simresults, connectivityclass)

connectivityclass = 'neighbors';

	centerneuron_index = 85;

	CM = simresults.networkParameters.connectivityMatrix;
		CM(find(eye(size(CM))))=0;
		[i j v] = find(CM);
		MN = full([i j v]);
		S = sortrows(MN,3);
		NPmost = S(end-5:end, [1 2]);
		NPleast = S(end-3:end, [1 2]);
		NPall = [NPmost ; NPleast];

			
		NPl = reshape(NPmost,1,[]);

		centerneuron_b = zeros (numneurons,1); centerneuron_b(centerneuron_index) = 1;
		neighbors = find(CM*centerneuron_b); 
		neighbors_b = zeros(numneurons,1);
		neighbors_b(neighbors) = 1;
		nextneighbors = find(CM*neighbors_b); 
		
		S = setdiff(nextneighbors, neighbors);

		stimulated = find(simresults{sims2p(1)}.perturbation.mask{1});

		% remove neurons that don't fire enough spikes
		lowfr = find(sum(simresults{sims2p(1)}.spikes.spikespercell')<20);

		
		switch selected_neurons

			case 'neighbors' % johny and his neighbors
				selectedneurons = setdiff(neighbors, stimulated);
				selectedneurons = setdiff(selectedneurons, lowfr);
				R = randperm(length(selectedneurons));
				nselneurons  = min(nselneurons, length(selectedneurons))
				selectedneurons = [centerneuron_index ; selectedneurons(R(1:nselneurons))];

			case 'nextneighbors'
				selectedneurons = setdiff(nextneighbors, stimulated);
				selectedneurons = setdiff(selectedneurons, lowfr);
				R = randperm(length(selectedneurons));
				nselneurons  = min(nselneurons, length(selectedneurons))
				selectedneurons = selectedneurons(R(1:nselneurons))

			case 'stimulated'
				R = randperm(length(selectedneurons));
				nselneurons  = min(nselneurons, length(selectedneurons))
				selectedneurons = setdiff(stimulated, lowfr);
				selectedneurons = stimulated(R(1:nselneurons));
				
			case 'acrosscluster'
				sel1 = setdiff(find(simresults{1}.W.stats.clusters==1), stimulated)  ;
				sel1 = setdiff(sel1, lowfr);
				sel2 = setdiff(find(simresults{1}.W.stats.clusters==2), stimulated);
				R1 = randperm(length(sel1));
				R2 = randperm(length(sel2));
				nselneurons  = min(nselneurons, length(R1))

				selectedneurons =  [sel1(R1(1)) ; sel2(R1(1:nselneurons/2))];
				selectedneurons = setdiff(selectedneurons, lowfr);

			case 'highgapconn'
				selectedneurons = setdiff(NPmost, stimulated);
				selectedneurons = setdiff(selectedneurons, lowfr);
				R = randperm(length(selectedneurons));
				nselneurons  = min(nselneurons, length(selectedneurons))
				selectedneurons = selectedneurons(R(1:nselneurons));

			case  'leastgapconn'
				selectedneurons = setdiff(NPleast, stimulated);
				selectedneurons = setdiff(selectedneurons, lowfr);
				R = randperm(length(selectedneurons));
				nselneurons  = min(nselneurons, length(selectedneurons))
				selectedneurons = selectedneurons(R(1:nselneurons));

			case 'all'
				selectedneurons = 1:length(CM);
				selectedneurons = setdiff(selectedneurons, lowfr);
				R = randperm(length(selectedneurons));
				selectedneurons = selectedneurons(R(1:nselneurons));


		end
