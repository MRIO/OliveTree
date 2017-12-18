function idx = spctcl(v, n, rep)
% [v, d] = spctcl(v, n, rep) computes spectral clustering based on
% given eigenvectors of the laplacian and number of clusters.
%
% Input:
%   v: Matrix whose columns are eigenvectors of the laplacian
%   n: Number of clusters
%   rep: Number of replicates in the kmeans step
%
% Output:
%   idx: Cluster label
%
% Written by Sungho Hong, CNS unit, OIST
% September 2017

y = v(:, (end-n+1):end);
ny = sqrt(sum(y.^2,2));
y = bsxfun(@times, y, 1./ny);
idx = kmeans(y, n, 'Replicates', rep);

end

