function [Z,H,om_opt,min_err]=co_prciter(Q,ngrid,fignum)
% DAMOCO Toolbox, function CO_PRCITER, version 28.02.14
% This function decomposes the given coupling function Q(phi1,phi2) 
% into a product Q-omega=Z(phi1)*H(phi2), using iteration algorithm
%
% Input:  Q        is the coupling function computed on a grid,
%         ngrid    is the grid size,
%         fignum   if it is not zero, then error is plotted in
%                  figure fignum
% Output: Z(phi1)  is PRC on the same grid
%         H(phi2)  is forcing on the same grid
%
niter=10; 
% parameters c1, c2, npt_om shall be adjusted after checking the plot 
% of decomposition error vs omega
c1=1; c2=1; 
om1=c1*min(Q(:));  om2=c2*max(Q(:)); 
npt_om=200;                        % number of points between om1 and om2                   
omega=om1:(om2-om1)/(npt_om-1):om2;
decomp_err=zeros(length(omega),1);             % Error of decomposition
%
maxq=max(abs(Q(:)));         % to find the maximal value
[kmax,n]=find(abs(Q)==maxq); 
if length(kmax)>1            % if there are several maxima, choose first
   kmax=kmax(1);
end
for step=1:length(omega)
    F=Q-omega(step);
    [Z,H]=decomp(F,kmax,ngrid,niter);
    F=F-Z'*H;
    decomp_err(step)=std(F(:));
end
[min_err,ind]=min(decomp_err);
om_opt=omega(ind);
[Z,H]=decomp(Q-om_opt,kmax,ngrid,niter);
if fignum > 0 
    Qn=co_gnorm(Q);
    figure(fignum); plot(omega, decomp_err/Qn,'ro');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,H]=decomp(F,kmax,ngrid,niter)
H0=F(kmax,:);   % reference profile where the function Q is maximal
for i=1:niter
    for k=1:ngrid
        Z(k)=trapz(F(k,:).*H0);
    end
    Z=Z/ trapz(H0.*H0);
    for k=1:ngrid
        H(k)=trapz(F(:,k).*Z');
    end
    H=H/ trapz(Z.*Z);
    H0=H;
end
end