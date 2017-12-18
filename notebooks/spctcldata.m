function [v, d] = spctcldata(z)
% [v, d] = spctcldata(z) computes eigenvectors and eigenvalues of a
% graph laplacian given an adjacency matrix
%
% Input:
%   z: Adjacency matrix
%
% Output:
%   v: Matrix whose columns are eigenvectors of the laplacian
%   d: Eigenvalues
%
% Written by Sungho Hong, CNS unit, OIST
% September 2017

cc = compute_cc(z);
ll = laplacian(cc);
[v, d] = eig(ll);
d = diag(d);
end

function cc = compute_cc(z)

cc = (pdist(z));
cc = exp(-squareform(cc/std(cc)/2).^2);
end

function ll = laplacian(cc)

dd2 = diag(sqrt(1./sum(cc)));
ll = dd2*cc*dd2;
ll = (ll+ll')/2;

end


