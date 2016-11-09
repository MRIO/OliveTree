function [Std_dPhi1, Std_dPhi1_synth, Std_Resid, dPhi1_synth, Resid] = ...
    co_resid_decomp(dPhi1, Phi1, Phi2, q)
% DAMOCO Toolbox, function CO_RESID_DECOMP, version 02.03.14
% The function co_resid_decomp computes the synthetic phase derivative,
% dPhi1_synth, from the fitted phase model. Form this the residual of the
% measured phase derivatve, dPhi1, is computet. In case of a perfect phase model, this residual should be zero, 
% while it grows when components of the measured phase are not captured in model. Hence, this function is a quality check of the phase model, q.  
% 
% Form of call: 
%          [Std_dPhi1, Std_dPhi1_synth, Std_Resid, dPhi1_synth, Resid] = co_resid_decomp(dPhi1, Phi1, Phi2, q)
%           
% Input:   dPhi1:   phase derivative 1st system
%          Phi1:    phase of the 1st system
%          Phi2:    phase of the second system
%          q:       phase model of the 1st system given on a grid
%          
% Output:  Std_dPhi1:       standarad deviation of the given phase derivative of the first system
%          Std_dPhi1_synth: standard deviation of the synthetic phase derivative of the first system, reconstructed from the phase model q.
%          Std_Resid:       standard deviation of the residual
%          dPhi1_synth:     the synthetic phase derivative, reconstructed from the phase model,q.
%          Resid:           the residual of the dPhi1

dPhi1_synth = dPhi(q, Phi1, Phi2); % Rekonstruction of dPhi from the coupling function q

Resid = dPhi1 - dPhi1_synth;

Std_dPhi1=std(dPhi1-mean(dPhi1));
Std_dPhi1_synth=std(dPhi1_synth-mean(dPhi1_synth));
Std_Resid = std(Resid-mean(Resid)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dPhi=dPhi(Q, phi1, phi2)
S=size(Q); S1= S(1)-1;
[X,Y]=meshgrid(2*pi*(0:S1)/S1, 2*pi*(0:S1)/S1);
dPhi= interp2(X,Y,Q,mod(phi2, 2*pi),mod(phi1, 2*pi));
end