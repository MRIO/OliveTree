% test_ensure_conn.m


% depth   = [1:10];
% breadth = [1:3];
% height  = [1:1];

% % # compute adjacency matrix
% [X Y Z] = meshgrid(depth,breadth,height);
% X = X(:); Y = Y(:); Z = Z(:);

% D = squareform( pdist([X Y Z], 'chebychev'));

function W = test_ensure_conn(W,D)


    desired = 10;
    tolerance = 1;
    disttol = 2;
    maxiter = 10;

if mean(sum(W,2)) >= desired + tolerance | mean(sum(W,2)) <= desired - tolerance
	withintolerance = 0;
else
	withintolerance = 1;
end

W = triu(W)+triu(W)';

for iter = 1:maxiter
	iter

    ID(iter,:) = sum(W,2);

    while ~withintolerance & iter < maxiter

        fewest = find(ID(iter,:) == min(ID(iter,:)));
        most   = find(ID(iter,:) == max(ID(iter,:)));

        % find closest nodes that are not connected to fewest


        candidate = (D(fewest(1),:) == min(D(fewest(1),:))) & ~W(fewest(1),:) 
		candidate = find(candidate);
		

        W(fewest,candidate) = 1;
        W(candidate,fewest) = 1;

        candidate = (D(most(1),:) == min(D(most(1),:))) & W(most(1),:);
		candidate = find(candidate)

        most   = find(ID(iter,:) == max(ID(iter,:)));

        W(most(1),candidate) = 1;
        W(candidate,most(1)) = 1;






        if mean(sum(W,2)) >= desired + tolerance | mean(sum(W,2)) <= desired - tolerance
			withintolerance = 0;
		else
			withintolerance = 1;
		end


    end
end



