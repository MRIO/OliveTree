function q1= co_kcplfct1(phi1,phi2,phi1_dot,ngrid,fignum,al_x,al_y)
% DAMOCO Toolbox, function CO_KCPLFCT1, version 04.03.14
%
% Given two phases (NOT protophases) and the derivative of the first one,
% the function yields the coupling function q1(phi1,phi2), using
% kernel smoothing technique.
%
% Form of call: q1= co_kcplfct1(phi1,phi2,phi1_dot,ngrid,fignum);
%               q1= co_kcplfct1(phi1,phi2,phi1_dot,ngrid,fignum,al_x);
%               q1= co_kcplfct1(phi1,phi2,phi1_dot,ngrid,fignum,al_x,al_y);
%
% The function is similar to CO_KCPLFCT2, just it computes one
%             function and has an incorporated option for plotting
% Parameters: phi1, phi2 are phases, phi1_dot is time derivative of phi1
%             ngrid is size of the grid for function computation;
%             fignum is the number of the graphic window for output;
%             if fignum=0, then the function is computed, but not
%             plotted. 
%             al_x,al_y are the smoothing parameters for the kernel
%             If they are omitted, the default values are used
%             If only first is given, then al_x=al_y is used
%
pi2=pi+pi;            ng1=ngrid-1; 
if nargin == 5       % smoothing factor, default value, al_x=al_y;
    al_x = ng1/pi2; al_x=al_x*al_x; al_y=al_x;
end
if nargin == 6       % one smoothing factor is given, al_x=al_y;
    al_y=al_x;
end

q1=zeros(ngrid,ngrid);                % function to be computed
Nrm=zeros(ngrid,ngrid);               % denominator of the expression for Q
exp_x=zeros(ng1,1);     exp_y=zeros(ng1,1);    % Kernels K_x and K_y 

dx=pi2/ng1;  % step in the argument of Q1(x,y)  
cdx=cos(dx); sdx=sin(dx);    % cosine and sine of dx are used for efficient 
                             % computation of cosine functions of (phi1-x) 
                             % and of (phi2-y) on the grid
for l=1:length(phi1)         % Main cycle over all points of the time series
    % computation of the kernel functions K_x and K_y at x=0, y=0
    cx=cos(phi1(l)); sx=sin(phi1(l));  cy=cos(phi2(l)); sy=sin(phi2(l));
    exp_x(1)=exp(al_x*(cx-1));   exp_y(1)=exp(al_y*(cy-1)); 
    for k=2:ng1              % Kernel function at other grid points
        cxn=cx*cdx+sx*sdx;  sxn=sx*cdx-cx*sdx; 
        cyn=cy*cdx+sy*sdx;  syn=sy*cdx-cy*sdx;
        cx=cxn; sx=sxn; cy=cyn; sy=syn; 
        exp_x(k)=exp(al_x*(cx-1)); exp_y(k)=exp(al_y*(cy-1)); 
    end
    for k=1:ng1              % Computing coupling function 
        for n=1:ng1
            Kernel=exp_x(k)*exp_y(n);
            q1(k,n)=q1(k,n)+phi1_dot(l)*Kernel;   % Numerator
            Nrm(k,n)=Nrm(k,n)+Kernel;          % Denominator
        end
    end
end                          % End of cycle over all points of the time series
q1(1:ng1,1:ng1)=q1(1:ng1,1:ng1)./Nrm(1:ng1,1:ng1);  %  Normalizing coupling function
for k=1:ng1                  % Filling last points on the grid using periodicity
    q1(ngrid,k)=q1(1,k); q1(k,ngrid)=q1(k,1);
end
q1(ngrid,ngrid)=q1(1,1); 
if fignum >0                 % Plot function, if required
    x=0:dx:pi2; 
    figure(fignum); 
    surf(x,x,q1); axis([0 pi2 0 pi2 -inf inf]); axis square; colormap jet;
    xlabel('\phi_1'); ylabel('\phi_2'); zlabel('q_1')
end
end