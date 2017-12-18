function [zn, idx] = make_rand_samples2(z, n, r)
% [zn, idx] = make_rand_samples2(z, n, r)
%
% Input:
%   z: NxN (partial) correlation matrix
%   n: Number of samples to generate
%   r: Fraction
%
% Output:
%   zn: 1xn cells containing resampled data
%   idx: 1xn cells containing index for zn
%
% Written by Sungho Hong, CNS unit, OIST
% September 2017

zn = cell(1, n);
idx = cell(1, n);

for i=1:n
    idx1 = (randn(1, size(z,1))<r);
    zn{i} = z(idx1, :);
    idx1 = int64(idx1);
    idx1(idx1>0) = 1:sum(idx1);
    idx{i} = idx1;
end

end
