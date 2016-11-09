function [Nrmq, omega] = co_gnorm(q)
% DAMOCO Toolbox, function CO_GNORM, version 04.03.14
%
% Given the coupling function q(phi1,phi2) on a grid
% this functions returns the norm and average of q
%
% Form of call: [Nrmq, omega] = co_gnorm(q)
%
% Input:        
%       q:      model of phase dynamics on grid
% Output:
%       Nrmq:   Norm of the coupling function
%       omega:  Estimate of the natural frequency

S=size(q);  ngrid=S(1); ng1_2=(ngrid-1)*(ngrid-1);
omega = mean(q(:));
qm = q-omega;
Nrmq=sqrt(trapz(trapz(qm.*qm))/ng1_2);
end