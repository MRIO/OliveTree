function [COR]=co_fcfcor(Qcoef1, Qcoef2) 
% DAMOCO Toolbox, function CO_FCFCOR, version 06.03.14
%
% This function computes the correlation between two coupling functions
% given in terms of Fourier coefficients Qcoef1, Qcoef2. 
% 
% Form of call: 
%               COR = co_fcor(Qcoef1, Qcoef2)
% Input:
%               Qcoef1, Qcoef2: The Fourier coefficients of both coupling functions
% Output:
%               COR:            correlation of the coupling functions 
S = size(Qcoef1); S = S(1);
N = (S-1)/2;
% Autonomous frequencies are deleted if they are still present
Qcoef1(N+1, N+1)=0;
Qcoef2(N+1, N+1)=0;
%Computing Correlation
autCor1 = real(trapz(trapz(Qcoef1.*conj(Qcoef1)))); % Autocorrelation
autCor2 = real(trapz(trapz(Qcoef2.*conj(Qcoef2)))); 
COR = trapz(trapz(Qcoef1.*Qcoef2(2*N+1:-1:1,2*N+1:-1:1) )) / sqrt( autCor2 .* autCor1 );
COR = real(COR); %Eliminating numercal imaginary rest


