function [phi, Cav] = co_avcyc(x, theta, N, PL, alpha)
% DAMOCO Toolbox, function CO_AVCYC, version 18.05.14
%
% Given an oscillatory periodic time series, x, and a first estimate 
% of its protophase, theta, the function yields a new protophase, phi, 
% based on the projetion of the embedded time series on the average cycle.
%
% Form of call: 
%                   [phi, Cav] = co_avcyc(x, theta, N, alpha, PL)
% Input:            x:          oscillatory time series
%                   theta:      first estimate of the protophase
%                   N:          maximal order of the Fourier expansion
%                   alpha:      defines the relative weight in the 
%                               minimization procedure: 
%                               (1) alpha=0: the new protophase, phi, is
%                                   based on minimizing the distance to the
%                                   average cycle only. alpha =0 is default
%                                   value.
%                               (2) weight alpha quantifies how the 
%                                   difference of phi to the first estimate, 
%                                   theta, is used as a second constraint.
%                                   Hence, for very large alpha,  phi=theta.
%                                   This is useful in case of time series
%                                   possessing flat regions where the change
%                                   of the signal is almost on the noise level. 
%                                   In these regions the signal does not 
%                                   contain information on the evolution of 
%                                   the phase. A method to optimize alpha is
%                                   described in the supplementary material 
%                                   in Ref. [5]. 
%                   PL:         If PL==1, some results are plotted:
%                                     (1) average cycle
%                                     (2) the coefficients Cav
%                                     (3) the derivative of new protophase, phi (green), 
%                                         and of the first estimate (red).
% Output:           phi:        the new protophase, based on projection on the average cycle. 
%                   Cav:        the coefficients of Fourier expansion of the average cycle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <5; alpha = 0; end;                      % default value  of alpha
if nargin < 4; PL=0; end;                          % default value of PL
S = size(x); if S(1)> S(2); x=x.'; end             % Check for correct format
S = size(theta); if S(1)> S(2); theta=theta.'; end % Check for correct format
x = x / max(x);                                    % Normalize amplitude of the time series
H = hilbert(x);                                    % Computing 2-dimensional embedding 
%                                                    by means of Hilbert transform 
H = H(1000:end-1000);
theta = theta(1000:end-1000);
[Cav] = av_comp(H, theta, N, PL);               % Computing coefficients of the average 
%                                                    cycle, Cav.
S=size(Cav); if S(1)> S(2); Cav=Cav.'; end         % Check for correct format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Computation of new protophase, phi, by means of
%%%%%%%%%%%%%%%%%%% projection on the average cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(1,length(theta));
DIS=mean(gradient(theta));                          % defining size of search region
for n = 1:length(theta)                      % Estimating new protophase, phi, by minimizing error
    phi(n) = fminbnd(@ERav, theta(n)-DIS, theta(n) + DIS,[], Cav, H(n), theta(n), alpha);   
end
if PL ==1;
    figure;
        plot(mod(phi,2*pi), gradient(phi),'g.','MarkerSize',4);
        hold on
            plot(mod(theta,2*pi), gradient(theta),'r.','MarkerSize',4);
        hold off
end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  SUBFUNCTION   [Cav] = co_av_comp(H, theta, N, PL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cav] = av_comp(H, theta, N, PL)
% DAMOCO Toolbox, function CO_AV_COMP, version 06.05.14
% Given a 2-dimensional complex embedding by means of Hilbert transform, H, and an estimate of the
% protophase, theta, the function yields the coefficients of the average cycle, Cav.
% Form of call: 
%          [Cav] = co_av_comp(H, theta, N, PL)
% Input:   H:       2-dimensional embedding by means of the Hilbert transform
%          theta:   estimate of the protophase
%          N:       maximal order of Fourier expansion,
%          PL:      if PL ==1 the results are plotted
% Output:  Cav      Fourier coefficients of the average cycle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=size(theta); if S(1) > S(2); theta = theta'; end  % Check for correct format
S=size(H); if S(1) > S(2);  H = H.'; end
Cav = zeros(N+1,1);
d_theta = gradient(theta);                          % Computing derivative of theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%     Computing the coefficients Cav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0 : N 
    Cav(n+1) =  sum(H.*exp(-1i*n*theta).*d_theta)/(theta(end)-theta(1)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Plotting results if PL ==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PL ==1;
    Y= zeros(1,10000);
    p2=2*pi*(0:9999)/10000;
    for n = 0 : N
        Y = Y + Cav(n +1)*exp(1i*n*p2);
    end
    figure;
    plot(H,'.','MarkerSize', 4);
    hold on;
        plot(Y,'.r','MarkerSize',4);
    hold off
    figure; plot(abs(Cav))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% ERRORFUNCTION  ERav(phi, Cav, H, theta, alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ER] = ERav(phi, Cav, H, theta, alpha)
N = (0:length(Cav)-1).';
Z = Cav * exp(1i*N*phi);
ER = (abs(H-Z).^2 + alpha* abs(exp(1i*phi)-exp(1i*theta) ).^2);
end