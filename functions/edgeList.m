function [ W_hat , W_bar, InvDeg ] = edgeList( W )
%this does that

W = double(W - diag(diag(W)));
degrees = sum(W>0 , 2);
InvDeg = diag(degrees.^-1);

W_bar = spalloc(length(W), sum(degrees), sum(degrees));
W_hat = spalloc(sum(degrees),length(W), sum(degrees)*2);
counter = 1;
for i = 1 : size(W,1)
    [~,H,~] = find(W(i,:));
    for j = 1 : degrees(i)
        W_hat(counter,H(j)) = W(i,H(j));
        W_hat(counter,i) = - W(i,H(j));
        counter = counter + 1;
    end
end

counter = 0;
for i = 1 : size(W,1)
    
    W_bar(i, 1 + counter : counter + degrees(i)) = ones(degrees(i),1);
    counter = counter + degrees(i);
end

W_hat = -full(W_hat);
W_bar = full(W_bar);


end

