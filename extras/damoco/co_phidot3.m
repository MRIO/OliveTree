function [phi1_dot,phi2_dot,phi3_dot,phi1,phi2,phi3] = co_phidot3(phi1,phi2,phi3,fsample)
% DAMOCO Toolbox, function CO_PHIDOT3, version 06.03.14
%
% CO_PHIDOT3 computes the derivatives (instantaneous frequencies)  
% of three phases using the Savitzky-Golay filter
% Parameter:   fsample is the sampling frequency
% All phases are truncated to avoid the boundary effect
%
norder=4;   % order of the fitting polynomial
sl=4 ;      % window semi-length
wl=2*sl+1;  % window length
[b,g] = sgolay(norder,wl);   % Calculate S-G coefficients
phi1=unwrap(phi1); phi2=unwrap(phi2); phi3=unwrap(phi3);
phi1=phi1(:); phi1_dot=phi1;  phi2_dot=phi1; phi3_dot=phi1; 
for n=sl+1:length(phi1)-sl
    phi1_dot(n)=g(:,2)'*phi1(n-sl:n+sl);
    phi2_dot(n)=g(:,2)'*phi2(n-sl:n+sl);
    phi3_dot(n)=g(:,2)'*phi3(n-sl:n+sl);
end
phi1_dot=phi1_dot(sl+1:end-sl)*fsample;
phi2_dot=phi2_dot(sl+1:end-sl)*fsample;
phi3_dot=phi3_dot(sl+1:end-sl)*fsample;
phi1=phi1(sl+1:end-sl); % Truncating both phases in order to 
phi2=phi2(sl+1:end-sl); % synchronize them with the derivative
phi3=phi3(sl+1:end-sl); % synchronize them with the derivative
end

