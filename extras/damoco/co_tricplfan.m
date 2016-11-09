function [COUP, NORM, OMEGA]=co_tricplfan(Qcoef1, Qcoef2, Qcoef3, meth, thresh)
% DAMOCO Toolbox, function CO_TRICPLFAN, version 12.03.14
%
% This function is a tool to analyze the coupling structure of
% a three oscillator network. Given the Fourier coefficients 
% Qcoef1, Qcoef2, Qcoef3 of the coupling functions, obtained by co_fcpltri,
% co_tricplfan returns the coupling structure in the form of matrix COUP.
%
% The n-th row contains measures of the effect ON the n-th oscillator.
% The m-th column contains measures of the effect CAUSED BY the m-th oscillator.
% For example, COUP(1,2) quantifies effect of node 2 on node 1, 2-->1
%                              0       COUP(1,2)       COUP(1,3)
%                          COUP(2,1)       0           COUP(2,3)
%                          COUP(3,1)   COUP(3,2)       0
% The diagonal is by definition empty, since it would present the effects of an 
% oscillator on itself; this place is however used as discussed below.
%
% COUP is computed by two different techniques, determined by the parameter 'meth' 
%
%   meth = 1 or meth = 2 : PARTIAL NORM
%
%       A non-diagonal element COUP(n,m) is partial norm, i.e. norm of the Fourier terms of the coupling
%       function Qcoefn, depending only on phases n,m (see Ref [4] for details).  
%       A diagonal value COUP(n,n) is the norm of the terms Qcoefn(n,n~0,k~0), i.e. of the terms, 
%       depending on both external phases. 
%       If diagonal terms are not small, the coupling cannot be explained by pairwise interaction 
%       only (see Ref [4]).
%
%
%   meth = 3 or meth = 4   : PARTIAL DERIVATIVE
%
%       The coefficients COUP are computed by means of the norm of the partial
%       derivative with respect to one of external phases. For example
%       the effect of the 2nd oscillator on the first is quantified as: 
%       COUP = sqrt(sum( m^2 Qcoef(n,m,k)^2) ). 
%       The drawback of this method is that noise is enhanced for large order n,m,k.
%       The diagonal terms here are zero.
%
% co_tricplfan also returns the norms of the all coupling functions in the 3-dimensional 
% array NORM and estimates of the autonomous frequencies of the oscillators in the 3-dimensional 
% array OMEGA. 
%
% To reduce numerical / noise effects, a threshold 'thresh' can be chosen in terms of percent 
% of the maximal absolute value of the coefficients Qcoef:
% all coefficients Qcoef with an absolute value smaller than thresh/100*max(abs(Qcoef)) 
% are set to zero. The default value is thresh=2%.
%
% From COUP, which quantifies the interaction of oscillators in absolute values, a relative 
% quantity can be defined by dividing COUP by the autonomous frequency OMEGA of the affected 
% oscillator: 
%           COUP_rel = COUP / OMEGA
% This measure quantify the strength of coupling in relation to autonomous frequency;
% hence, it quantify the tendency to synchronization.
% Output COUP_rel can be chosen by setting meth = 2 or meth = 4.
%
% Form of call: [COUP, NORM, OMEGA]=co_3osccoup(Qcoef1, Qcoef2, Qcoef3) [all parameters = default]
%               [COUP, NORM, OMEGA]=co_3osccoup(Qcoef1, Qcoef2, Qcoef3, meth): specifying method
%               [COUP, NORM, OMEGA]=co_3osccoup(Qcoef1, Qcoef2, Qcoef3, meth, thresh): specifying
%                                                                         both method and threshold
%               
% INPUT:
%       Qcoef1:         matrix of Fourier coefficients of the coupling function 
%                       of the 1st oscillator, produced by co_fcpltri
%       Qcoef2:         .. of the 2nd oscillator..
%       Qcoef3:         .. of the 3rd oscillator..
%       meth:           (1) (default)= partial norm, (2)= partial norm / omega, (3) = partial derivative, 
%                       (4) = partial derivative / omega 
%       thresh:         threshold for reduction of noise effects (also of computational noise) in terms 
%                       of percent of the maximal absolute value of the coefficients Qcoef: 
%                       thresh = 10 is good choice for dominant components only, the default
%                       value is thresh = 2%; typically this removes numerical artifacts, while keeping 
%                       all other terms. Feel free to vary this parameter [0,..100], to remove noise, etc..
%
%   OUTPUT:
%       COUP            matrix describing the coupling structure using the chosen method meth 
%       NORM:           norms of all coupling functions
%       OMEGA:          the autonomous frequencies of the 3 oscillators
%
if nargin == 4; 
    thresh = 2;                       % defining default value thresh=2%
end;
if nargin == 3;
    meth = 1;                         % default method =1: partial norm 
    thresh = 2;                       % defining default value thresh=2; 
end;
N = size(Qcoef1); N = (N(1)-1)/2;
COUP = zeros(3,3);
NORM = zeros(1,3);

OMEGA(1) = abs( Qcoef1(N+1,N+1,N+1)); Qcoef1(N+1,N+1,N+1)=0; % extracting autonomous frequencies
OMEGA(2) = abs( Qcoef2(N+1,N+1,N+1)); Qcoef2(N+1,N+1,N+1)=0;
OMEGA(3) = abs( Qcoef3(N+1,N+1,N+1)); Qcoef3(N+1,N+1,N+1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if meth<3; % choosing method
    [COUP(1,2), COUP(1,3), COUP(1,1), ~, ~, NORM(1)] = co_3to2(Qcoef1,N,thresh);
else
    [~, ~, ~, COUP(1,2), COUP(1,3), NORM(1)] = co_3to2(Qcoef1,N,thresh);
end
if meth==2 || meth==4;
  COUP(1,2) = COUP(1,2) / OMEGA(1);
  COUP(1,3) = COUP(1,3) / OMEGA(1);
  COUP(1,1) = COUP(1,1) / OMEGA(1); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if meth<3;
    [COUP(2,3), COUP(2,1), COUP(2,2), ~, ~, NORM(2)] = co_3to2(Qcoef2,N,thresh);
else
    [~, ~, ~, COUP(2,3), COUP(2,1), NORM(2)] = co_3to2(Qcoef2,N,thresh);
end
if meth==2 || meth==4;
  COUP(2,3) = COUP(2,3) / OMEGA(2);
  COUP(2,1) = COUP(2,1) / OMEGA(2);
  COUP(2,2) = COUP(2,2) / OMEGA(2); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if meth<3;
    [COUP(3,1), COUP(3,2), COUP(3,3), ~, ~, NORM(3)] = co_3to2(Qcoef3,N,thresh);
else
    [~, ~, ~, COUP(3,1), COUP(3,2), NORM(3)] = co_3to2(Qcoef3,N,thresh);
end
if meth==2 || meth==4;
  COUP(3,1) = COUP(3,1) / OMEGA(3);
  COUP(3,2) = COUP(3,2) / OMEGA(3);
  COUP(3,3) = COUP(3,3) / OMEGA(3); 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C12, C13, C123, PART12, PART13, TOT] = co_3to2(Q, N, thresh)  % this function sorts the coefficients and computes the norms
thresh=max(abs(Q(:)))*thresh/100;
ind= abs(Q)<thresh;   %all elements below threshold are set to zero
Q(ind)=0;
or1 = N+1; or21=2*N+1;
Q12 = zeros(or21,or21);
Q13 = zeros(or21,or21);
C123=0;
PART12=0;
PART13=0;
TOT=0;
for n = 1 : or21; % sorting coefficients
    for m = 1 : or21; 
        Q12(n, m) = Q(n,m,or1);
        Q13(n, m) = Q(n,or1,m);
            for k = 1:or21;
                if  m~=or1 && k~=or1;
                    C123 = C123 + abs(Q(n,m,k)^2);
                end;
                PART12= PART12 + (m-or1)^2*abs(Q(n,m,k)^2); % computing partial derivatives
                PART13= PART13 + (k-or1)^2*abs(Q(n,m,k)^2);
                TOT = TOT + abs(Q(n,m,k).^2); % computing norm of complete coupling function
            end;
    end;
end
C12=sqrt(sum(sum(abs(Q12).^2))); % computing partial norms
C13=sqrt(sum(sum(abs(Q13).^2)));
C123=sqrt(C123);
TOT=sqrt(TOT);
PART12=sqrt(PART12);
PART13=sqrt(PART13);
end