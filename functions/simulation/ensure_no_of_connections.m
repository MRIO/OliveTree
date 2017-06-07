% ensure_no_of_connections.m
function out =  ensure_no_of_connections(W, meanconn)
	% prune connectivity matrix such that number of connections is meanconn*length(W)
	% assumes absence of self connections!


   desiredconn =  meanconn * length(W);
   [i j v] = find(triu(W.*(ones(size(W))-eye(size(W)))));
   idmat = [i j];

   delta = desiredconn/2 - length(v);

   if delta>0 % we need more connections
    %to implement
    warning('could not ensure connectivity. Need more connections to start with. Try increasing meanconn.')
    return

   elseif delta < 0 % too many connections

    randidx = randperm(length(v));
    toprune = -delta;
warning off
    W(sub2ind(size(W), i(randidx(1:toprune)) , j(randidx(1:toprune) ))) =0;
    W(sub2ind(size(W), j(randidx(1:toprune)) , i(randidx(1:toprune) ))) =0;
warning on

    disp(['pruned ' num2str(toprune*2) ' connections of total ' num2str(length(v)*2)]);


   else
       
   end
    


out = W;