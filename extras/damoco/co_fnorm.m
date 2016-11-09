function [Nrmq, omega] = co_fnorm(Qcoef)
% DAMOCO Toolbox, function CO_FNORM, version 05.03.14
% Given the coefficients Qcoef of the Fourier expansion of the coupling
% function q, this functions returs the norm of q
% 
% Form of call: Nrmq = co_fbnorm(Qcoef)
% Input:        Qcoef are Fourier coefficients of the function
%
S=size(Qcoef);
N = (S(1)-1) / 2;
omega = real(Qcoef(N+1,N+1));            % defining estimate of natural frequency as Qcoef_0_0
Qcoef(N+1,N+1)=0;                        % Setting the constant term, omega, to zero  
Nrmq = sqrt(trapz(trapz (abs(Qcoef).^2)));   % Computing the norm of coupling function 
end
