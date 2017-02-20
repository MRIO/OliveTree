function c = spect_clust(cc, ncl_limit)
% c = spect_clust(cc, ncl_limit) 
% Spectral clustering from the affinity matrix. See von Luxburg, 2007 for the detail.
%
% Input:
% 	cc: (N,N) affinity matrix whose elements range from 0 to 1. The diagonal elements will be ignored.
%   ncl_limit: Maximal number of clusters
% 
% Output:
%   c: result data structure
%   c.labs: (N,ncl_limit) matrix containing cluster labels
%   c.v: eigenvectors of the laplacian
%   c.d: eigenvalues of the laplacian
%   c.gc : eigenvalue gap...
%   c.nc : number of clusters corresponding c.gc
%
% Written by Sungho Hong, CNU, OIST
% 2014

cc0 = cc-diag(diag(cc));  % The diagonals are ignored.
%cc0 = (cc0+cc0')/2;
degree = sum(cc0,2);
%dmat = diag(sqrt(1./degree));
%lapl = eye(size(cc0,1))-dmat*cc0*dmat; % Laplacian
dmat = diag((1./degree));
lapl = eye(size(cc0,1))-dmat*cc0; % Laplacian
%lapl = (lapl+lapl')/2; % lapl should be symmetric
[v, d] = eig(lapl);
d = real(diag(d));
% d = d(end:-1:1);
% v = v(:,end:-1:1);
cl_labs = ones(size(cc,1), ncl_limit);

nc = zeros(ncl_limit,1);
gc = zeros(ncl_limit,1);

for i=1:ncl_limit
    [g1, n1] = max(diff(d.^i));
    nc(i) = n1+1;
    gc(i) = g1;
end

for i=2:ncl_limit
    yi = v(:,1:i);
    yi = bsxfun(@times, yi, 1./sqrt(sum(yi.^2,2)));
%     cl_labs(:,i) = cluster(linkage(yi,'ward'),'maxclust',i); 
    cl_labs(:,i) = kmeans(yi,i,'replicates',10); 
end

c = struct('label', cl_labs, 'v', v, 'd', d, 'nc', nc, 'gc', gc);


end
