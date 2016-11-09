function dirin=co_dirpar(Fcoef1, Fcoef2)
% DAMOCO Toolbox, function CO_DIRPAR, version 26.02.14
% Given the Fourier coefficients of the coupling functions, 
% CO_DIRPAR returns the directionality index dirin, computed 
% via the partial derivatives of the coupling function with 
% respect to the external protophase.
% The index is defined in a way that for symmetrical bidirectional 
% coupling dirin=0 holds, while purely 
% unidirectional coupling 1->2 yields dirin=1;
% unidirectional coupling 2->1 yields dirin=-1.  
% 
% Form of call:     dirin=co_dirpar(Fcoef1, Fcoef2)
%                  
% Input:
%       Fcoef1, Fcoef2:  Fourier coefficients for the model of phase
%                        dynamics for both systems
% Output:
%       dirin:           directionality index
%
S=size(Fcoef1); or=(S(1)-1)/2;
NP1=0;
NP2=0;
for n= -or : or;
    for m= -or : or;
        NP1 = NP1 + m^2*abs(Fcoef1(n+or+1,m+or+1))^2;
        NP2 = NP2 + m^2*abs(Fcoef2(n+or+1,m+or+1))^2;
    end;
end;
nrm1=sqrt(NP1);
nrm2=sqrt(NP2);
if nrm1+nrm2 < 0.02
    disp('Warning: the coupling is very weak or the systems are not coupled!');
    disp('Result on directionality index may be not reliable');
end
dirin= (nrm2 - nrm1) / (nrm1 + nrm2);
end
