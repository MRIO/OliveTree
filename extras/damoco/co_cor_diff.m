function [COR,DIFF] = co_cor_diff(Q1,Q2)
% DAMOCO Toolbox, function CO_COR_DIFF, version 28.02.14
%
% The function computes correlation and a normalized measure 
% of the difference between two coupling functions on a grid
% (for the definition of the measure see Ref 5.)
% Input: Q1 and Q2 are coupling functions on a grid of same size
%
Q1 = Q1 - mean(Q1(:)); % Eliminating mean value of functions    
Q2 = Q2 - mean(Q2(:));
Q1 = Q1(1:end-1,1:end-1);
Q2 = Q2(1:end-1,1:end-1);
autCor1 = sqrt(trapz(trapz(Q1.^2)));    % Autocorrelation
autCor2 = sqrt(trapz(trapz(Q2.^2)));    % Autocorrelation
COR = trapz(trapz(Q1.*Q2)) / (autCor1 * autCor2);  % Correlation
% Normalized difference measure
DIFF = sqrt(trapz(trapz((Q1-Q2).^2))) / (autCor1 + autCor2);  
end