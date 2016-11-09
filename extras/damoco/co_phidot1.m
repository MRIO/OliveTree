function [phi_dot,phi] = co_phidot1(phi,fsample)
% DAMOCO Toolbox, function CO_PHIDOT1, version 06.03.14
%
% CO_PHIDOT1 computes the derivative (instantaneous frequency)  
% of the phase using the Savitzky-Golay filter
% Parameter:   fsample is the sampling frequency
% Phase is truncated to avoid the boundary effect
%
norder=4;   % order of the fitting polynomial
sl=4 ;      % window semi-length
wl=2*sl+1;  % window length
[b,g] = sgolay(norder,wl);   % Calculate S-G coefficients
phi=unwrap(phi);
phi=phi(:); phi_dot=phi; 
for n=sl+1:length(phi)-sl
    phi_dot(n)=g(:,2)'*phi(n-sl:n+sl);
end
phi_dot=phi_dot(sl+1:end-sl)*fsample;
phi=phi(sl+1:end-sl);   % Truncating the phase in order to
                        % synchronize it with the derivative               
end

