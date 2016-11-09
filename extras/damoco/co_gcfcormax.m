function [COR, Delta_self, Delta_ext, q2_shift]=co_gcfcormax(q1,q2)
% DAMOCO Toolbox, function CO_GCFCORMAX, version 06.03.14
%
% This function computes the correlation between two coupling functions
% given on a grid. This correlation is  defined as the maximum over possible 
% relative shifts of the coupling functions. 
% q1 serves as the reference function, while q2 is shifted. 
% 
% Form of call: 
%               [COR] = co_gcfcormax(q1, q2)
%               [COR, Delta_self] = co_gcfcormax(q1, q2) 
%               [COR, Delta_self, Delta_ext]=co_gcfcormax(q1, q2) 
%               [COR, Delta_self, Delta_ext, q2_shift]=co_gcfcormax(q1, q2) 
% Input:
%       q1:             The coupling function of the first oscillator, used as
%                       reference. 
%                      
%       q2:             The coupling function of the 2nd oscillator. This one is shifted
%                       in the test, so the shifts reported refer to this function. 
%                     
% Output:
%       COR:            Maximal correlation of the coupling functions
%       Delta_self:     Phase shift of the own phase of q2 which provides
%                       the maximal correlation of two coupling functions.
%       Delta_ext:      Phase shift of the external phase of q2 which provides
%                       the maximal correlation of the two coupling functions.
%       q2_shift:       q2 shifted according to Delta_self, Delta_ext, so that
%                       correlation of q1 and q2_shift is maximal without
%                       shift.
%
q1 = q1(1:end-1,1:end-1);    % Deleting last points in the matrix; 
q2 = q2(1:end-1,1:end-1);    % by convention these are identical to the first one.
q1 = q1 - mean(q1(:));       % Deleting constant = natural frequency
q2 = q2 - mean(q2(:));

ngrid = size(q1); ngrid = ngrid(1); % Computing size of the grid
Q2 =[q2 q2; q2 q2]; %Definig periodic function on a larger grid, to be shifted
autCor1 = trapz(trapz(q1.*q1)); % Autocorrelation
autCor2 = trapz(trapz(q2.*q2)); % Auocorrelation
m1=0; m2=0; COR=0;
for n = 0 : ngrid - 1;
    for m = 0:ngrid -1;
        A = Q2(ngrid+1-n:2*ngrid-n, ngrid+1-m:2*ngrid-m); % Shifting Q2
        a = trapz(trapz(q1.*A)) /(sqrt(autCor1*autCor2));   % Correlation for all different shifts
            if a > COR;                                     % Selecting maximum of correlations
                m1 = n;
                m2 = m;
                COR = a;
                q2_shift = A;
            end;
    end
end;
Delta_ext = mod(2*pi - 2*pi*(m2) / ngrid, 2*pi); % Computinng shift in terms of 2pi
Delta_self  = mod(2*pi - 2*pi*(m1) / ngrid, 2*pi);

q2_shift(:,end+1)=q2_shift(:,1); % Adding last points using the  periodicity.
q2_shift(end+1,:)=q2_shift(1,:); 
q2_shift(end,end)=q2_shift(1,1); 
end
