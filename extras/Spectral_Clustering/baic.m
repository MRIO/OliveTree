function [bic aicc aic] = baic(data, c)
% [BIC, AICC, AIC] = baic(data, c)
% Information criterions for the clustering result.
%
% Input:
% 	data: (N,L) original data.
%	c: clustering result
%
% Output:
%	BIC: Bayesian information criterion
%	AICC: Akaike information criterion with the data size correction
%   AIC: Akaike information criterion 
%
% Written by Sungho Hong, CNU, OIST
% 2014

cls = c.label;

n = size(cls,2);
bic = zeros(n,1);
aicc = bic;
aic = bic;

for i=1:n
    [bic(i), aicc(i), aic(i)] = baic1(data, cls(:,i)); % compute ICs for every partitionings in c
end

end

function [bic aicc aic] = baic1(data, cl)

k = numel(unique(cl));
n = size(data,1);
labs = unique(cl); % find all the labels

sig2 = 0;
for i=1:numel(labs)
    data1 = data(cl==labs(i),:);
    resid = bsxfun(@minus, data1, mean(data1,1));
    sig2 = sig2 + sum(sum(resid.^2)); % here we compute SSE.
end
sig2 = sig2/n;

aic = n*log(sig2) + 2*k;
bic = n*log(sig2) + k*log(n);
aicc = aic + 2*k*(k+1)/(n-k-1);

end
