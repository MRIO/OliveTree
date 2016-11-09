function [phi,arg,sigma] = co_fbtrT(theta,ngrid)
% DAMOCO Toolbox, function CO_FBTRT, version 02.03.14
%
% Fourier series based transformation protophase theta --> phase phi 
% for one oscillator, with optimization according to 
% C. Tenreiro, J. Nonparametric Stat, v 23, 533, 2011
%
% Form of call:  co_fbtrT(theta); 
%                co_fbtrT(theta,ngrid); 
% Input parameters:  
%                theta is the protophase
%                ngrid is the grid size for computation of the 
%                      transformation function sigma
% Output:  phi = co_fbtrT(...) if only transformation is required.
% Output:  [phi,arg,sigma] = co_fbtrT(...) if also the transformation
%          function sigma is required; it can be plotted as
%          plot(arg,sigma); sigma is computed on the grid.
%          Default grid size is 50.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft=100;           % Maximal number of Fourier coefficients
Spl=zeros(nfft,1);  % Fourier coefficients 1,...,nfft
Hl=zeros(nfft,1);   % Tenreiro function to be minimized

IN = find(diff(mod(theta,2*pi))<0);  % correction for short time series:
npt=length(theta(IN(1) : IN(end)));  % only full periods are used

S=0; c=double(npt+1)/double(npt-1);
for k=1:nfft        % computing Fourier coefficients
    Spl(k)=sum(exp(-1i*k*theta(IN(1) : IN(end))))/npt;
    S=S+Spl(k)*conj(Spl(k))-1./double(npt);
    Hl(k)=k/npt-c*S;   % Tenreiro function
end
[~,indopt]=min(Hl);

phi=theta;     % Transformation theta --> phi
if nargout==3  % sigma is computed along with the transformation
  if nargin < 2, ngrid=50; end 
  arg=0:(ngrid-1); arg=arg*pi*2/(ngrid-1); arg=arg';
  sigma=ones(ngrid,1);
  for k=1:indopt
      sigma=sigma + 2*real(Spl(k)*exp(1i*k*arg));
      phi=phi+2*imag(Spl(k)*(exp(1i*k*theta)-1)/k);
  end
else           % only transformation is required
  for k=1:indopt
      phi=phi+2*imag(Spl(k)*(exp(1i*k*theta)-1)/k);
  end
end      
end
  