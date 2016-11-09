function [COUP] = co_nettri(PHI, PHI_dot, N, meth, thresh)
% DAMOCO Toolbox, function CO_NETTRI, version 12.03.14
%
% This function is a tool to analyze the coupling structure of
% a network of N>3 coupled oscillators. 
% Given the phases of all oscillators, PHI, and their derivatives, PHI_dot,
% co_nettri yields the coupling structure in form of a coupling matrix: COUP.
% Form of call:
%               [COUP] = co_nettri(PHI, PHI_dot, N, meth, thresh)
%               [COUP] = co_nettri(PHI, PHI_dot, N, meth): thresh = default
%               [COUP] = co_nettri(PHI, PHI_dot, N): meth=default, thresh = default
% Input:
%               PHI:        A matrix containing all phases of the network, the matrix is of the size
%                           (number of nodes) times (number of time points)
%               PHI_dot:    A matrix containing the corresponding
%                           derivatives of the phases
%               N:          Order of Fourier series of the phase models. Note, that the computational time 
%                           grows  ~N^3. An order of N=5 is already sufficient for quite complex coupling functions.
%               meth:       Choose method of analysis: meth=1  (default) ==> partial norm, 
%                                                      meth=2  ==> partial norm / omega, 
%                                                      meth=3  ==> partial derivative, 
%                                                      meth=4  ==> partial derivative / omega.  
%                                                      For  details see: help to co_tricplfan and Ref. [4]
%               thresh:     threshold for reduction of noise effects (also of computational noise) in terms 
%                           of percent of the maximal absolute value of the coefficients Qcoef: 
%                           thresh = 10 is good choice for dominant components only, the default
%                           value is thresh = 2%; typically this removes numerical artifacts, while keeping 
%                           all other terms. Feel free to vary this parameter [0,..100], to remove noise, etc..
% Output:       COUP:       For a network of M oscillators, this M x M matrix contains measures of coupling strengths.
%                           The n-th row contains measures of the effect ON the n-th oscillator.
%                           The m-th column contains measures of the effect CAUSED BY the m-th oscillator. 
%                           Note, that this is a directional measure:  so generally COUP(n,m) not equal COUP(m,n).
%
if nargin == 4; 
    thresh = 2;                       % defining default value thresh=2%
end;
if nargin == 3;
    meth = 1;                         % default method =1: partial norm 
    thresh = 2;                       % defining default value thresh=2; 
end;
% correct orientation of the input
S=size(PHI);
if S(2)<S(1); PHI=PHI';  end
S=size(PHI_dot);
if S(2)<S(1); PHI_dot=PHI_dot';  end
M = min([S(1) S(2)]);
C = zeros(M,M,M-2);
IN = ones(M,M);
for n = 1:M;
    for m = n+1 : M;
        for k = m+1:  M;
            [~, maxind]=co_maxsync3(PHI(n,:), PHI(m,:), PHI(k,:), N); % check for bad synchrony
            if maxind > 0.5; display('Warning: oscillators are to close to synchrony, results might be not reliable!'); end
            [Qcoef1, Qcoef2, Qcoef3]=co_fcpltri(PHI(n,:), PHI(m,:), PHI(k,:), PHI_dot(n,:), PHI_dot(m,:), PHI_dot(k,:), N); % computing coupling functions of tripletts
            [COUP] = co_tricplfan(Qcoef1, Qcoef2, Qcoef3, meth, thresh); % computing coupling matrix of tripletts
            C(n, m ,IN(n,m)) = COUP(1,2); IN(n,m) = IN(n,m)+1; % sorting all results in a corresponding MxM matrix. Note, for each entry we compute M-2 values, since each pair of oscillators appears in M-2 tripletts.
            C(m, n ,IN(m,n)) = COUP(2,1); IN(m,n) = IN(m,n)+1;
            C(n, k ,IN(n,k)) = COUP(1,3); IN(n,k) = IN(n,k)+1;
            C(k, n ,IN(k,n)) = COUP(3,1); IN(k,n) = IN(k,n)+1;
            C(m, k ,IN(m,k)) = COUP(2,3); IN(m,k) = IN(m,k)+1;
            C(k, m ,IN(k,m)) = COUP(3,2); IN(k,m) = IN(k,m)+1;  
        end
    end
end
for n = 1:M;
    for m = 1:M
       COUP(n,m)=min(C(n,m,:)); % finding minimum over all triplets for each cell of the MxM matrix.
    end
end
end

