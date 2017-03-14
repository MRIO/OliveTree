function out = get_neighbors(W, neurons, order, varargin)
% W is a square matrix
% neurons is a vector of at least two elements
% 
W = triu(W);

Ind = zeros(length(W),1);
non_overlap = cell(length(neurons), order);

o = 1;
for n = neurons
	
	Ind(n) = 1;
	neigh{n, o} = find(W * Ind);
	non_overlap{n, o} = [n neigh{n, o}'];


	Ind = Ind*0;
	for o = 2:order
		Ind(non_overlap{n,o-1}) = 1;
		
		neigh{n, o} = find(W * Ind);
		non_overlap{n, o}= unique([non_overlap{n, o-1} neigh{n, o}']);

	end
	Ind = Ind*0;

end


all_neighbors = cell(1,o);
for o = 1:order
	
		for n = neurons
			all_neighbors{o} = [all_neighbors{o} non_overlap{n,o}];
		end
	
	unique_neigh{o} = unique(all_neighbors{o});
	neighbors_per_order.shared{o} = intersect(all_neighbors{o}, unique_neigh{o});
end

out.neighbors_per_order = unique_neigh;
out.neighbors_per_neuron_per_order = non_overlap;
