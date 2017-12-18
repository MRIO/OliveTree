function [dsts, rdsts] = cij2(clidx, idxz)
% function [r, rr] = cij2(clidx, idxz) normalized (cosine) distance
% between each pair of clustering results.
%
% Input:
%   clidx: 1xn cells containing clustering labels
%   idxz:  1xn cells containing the indices used for random sampling
%
% Output:
%   dsts: distance results
%   rdsts: distances with randomly permuted labels, for control
%
% Written by Sungho Hong, CNS unit, OIST
% September 2017

n = length(clidx);
npair = n*(n-1)/2;
dsts = zeros(1, npair);
rdsts = zeros(1, npair);

k = 1;
for i=1:(n-1)
    idxzi = (idxz{i}>0);
    for j=(i+1):n
        idx_common = (idxzi.*(idxz{j}>0))>0; % use only the points present in both data sets.
        idi = idxz{i}(idx_common);
        idj = idxz{j}(idx_common);
        cli = clidx{i}(idi);
        clj = clidx{j}(idj);
        ciji = double(pdist(cli)==0);
        cijj = double(pdist(clj)==0);
        dsts(k) = (1-dot(ciji, cijj)/norm(ciji)/norm(cijj));

        cli = cli(randperm(length(cli)));
        clj = clj(randperm(length(clj)));
        ciji = double(pdist(cli)==0);
        cijj = double(pdist(clj)==0);
        rdsts(k) = (1-dot(ciji, cijj)/norm(ciji)/norm(cijj));

        k = k+1;
    end
end

end
