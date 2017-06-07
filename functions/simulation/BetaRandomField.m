function [ F ] = BetaRandomField( rows,columns,correlation,FiltVal, alpha,beta )
% builds a Gaussian Markov Random Field with marginal Beta distribution
% code based on randomfield by Qiqi Wang and Paul G. Constantine

% alpha and beta are Beta distribution parameters

%   filter:     Set filter to a number between 0 and 1 to capture a
%               percentage of the energy of the field as determined by the
%               square roots of the eigenvalues of the covariance matrix.
%               (Default 1 means no filter.)

corr.name = 'exp';% exponential decay of the correlation with distance
% It's possible to specify a Gaussian ('gauss') decay of the correlation.
% The Exponential  decay has thicker tails than 'gauss'

corr.c0 = correlation; % isotropic (rotationally invariant) correlation a row vector like [0.3, 0.3];
% It's possible to specify anisotropic correlations, i.e. stretching
% preferrentially in one direction.

x = 1:rows;
y = 1:columns;
[X,Y] = ndgrid(x,y); mesh = [X(:) Y(:)]; % 2-D mesh

% set a spatially varying variance (must be positive!)
corr.sigma = 1;
mu = 0;
[F,~] = randomfield(corr,mesh,'filter',FiltVal);
F = normcdf(F,mu, corr.sigma);
F = betainv(F,alpha,beta);
F = reshape(F,rows,columns);
end