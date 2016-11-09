function [dgGABA_dt] = computefastGABA_state(stim,g)
% differential equations governing chnange in conductance of GABA_A.
% The same function is used for the spine and for the dendritic GABA_A
% receptors, to ensure that they are identical.

% Warning:
% TIME is seconds! To convert to ms, multiply outgoing derivatives by
% 1e-3.

alpha = 9.82337423115800e-07;
mu = 6.06975598757552e-05;
tau = 0.267272442329587;
gamma = 0.00192493295896767;
Beta_offset = 37e-3;
Beta_range = 32e-3;
% Ca_0 = 2e-3;


Opening= (3*1e-3 )^-1;%channel opening rate \alpha = 2msec^-1
% Closing = (Beta_offset + Beta_range * (1 / (1 + exp( -  alpha^-1 * (Ca - mu) ) ) )  )^-1; %channel closing rate \beta = 10msec^-1
Closing = (10e-3)^-1 ; %channel closing rate \beta = 10msec^-1
dgGABA_dt = stim * Opening*( 1 - g ) - Closing * g;

end