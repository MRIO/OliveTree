function [COR]=co_gcfcor(q1,q2)
% DAMOCO Toolbox, function CO_GCFCOR, version 06.03.14
%
% This function computes the correlation between two coupling functions
% given on a grid. 
% 
% Form of call: 
%               [COR] = co_gcfcor(q1, q2)
% Input:
%       q1:             The coupling function of the first oscillator, it is used as
%                       the reference. 
%                      
%       q2:             The coupling function of the 2nd oscillator. This one is shifted
%                       in the test, so the shifts reported refer to this function. 
%                     
% Output:
%       COR:            Maximal correlation of the coupling functions
%
q1 = q1(1:end-1,1:end-1);    % Deleting last points in the matrix; 
q2 = q2(1:end-1,1:end-1);    % by convention these are identical to the first one.
q1 = q1 - mean(q1(:));       % Deleting constant = natural frequency
q2 = q2 - mean(q2(:));

autCor1 = trapz(trapz(q1.*q1)); % Autocorrelation
autCor2 = trapz(trapz(q2.*q2)); % Auocorrelation
COR = trapz(trapz(q1.*q2)) /(sqrt(autCor1*autCor2)); % Computing Correlation
             
end
