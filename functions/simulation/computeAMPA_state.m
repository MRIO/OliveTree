function [dAMPAS_dt] = computeAMPA_state( stim, g )
% differential equations governing change in conductance of AMPA.

% Warning:
% TIME was seconds! To convert to ms, multiplying outgoing derivatives by
% 1e-3.

%AMPA 5nS/1mum^3
% conductance value taken from O'Donnell v.Rossum 2011 J. Neurosci
Alpha=(.18*1e-3)^-1;%channel opening rate \alpha = 2msec^-1
Beta =(1.8*1e-3)^-1; %channel closing rate \beta = 10msec^-1
% POpen=(1/(1 + exp(- (stim))));
% dAMPAS_dt = Alpha*POpen*(1-g)-Beta*g;
dAMPAS_dt = stim*Alpha*(1-g)-Beta*g;

end