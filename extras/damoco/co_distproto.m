function [theta, Start, Stop]=co_distproto(x, NV)
% DAMOCO Toolbox, function CO_DISTPROTO, version 10.04.14
%
% Given a 2-dim embedding x of a time series x, the function computes an
% estimate of the protophase, theta , based on the ratio of the length 
% of the trajectory in the phase plane, see Ref [3] for details.
%       
%                    !!!!   WARNING   !!!!
% Note, that this method is very sensitive to amplitude variations and, 
% thus, yield only rough estimate of the prophase.
% It should be used with caution and only when the amplitude variation 
% is weak and no other methods are available, e.g., if the trajectory has
% self-crossing loops
%
% Form of call: co_distproto(x,NV)
%               
% Input:        x          : 2-dim embedding of the analyzed time series;
%                            x can be an array of complex numbers, e.g. 
%                            produced by the MATLAB function 'hilbert', 
%                            or a 2-dim array of real numbers (two 
%                            coordinates of a point in the phase plane).
%               NV         : vector NV = [a, b] is normal to the Poincare 
%                            section, which defines the beginning of 
%                            each period.
% Output:
%               theta      : estimate of the protophase
%                            while x has N points x(1), .., x(N),
%                            the protophase theta is defined as 
%                            theta(Start), theta(Start+1), .., theta(Stop),
%                            where 1< Start <Stop < N
%               Start      : Start and Stop correspond to the beginning of the 
%               Stop         first complete cycle and to the end of the last one, 
%                            respectively, i.e., to the first and last crossing 
%                            of the Poincare section
%
if isreal(x)==0; % Check input
    y(1,:)=real(x);
    y(2,:)=imag(x);
else
    y=x;
end
x=[];
S=size(y); if S(1)>S(2); y=y'; end
y(1,:) = y(1,:)/std(y(1,:));
y(2,:) = y(2,:)/std(y(2,:));
Pro = zeros(1,length(y)); % Allocating space
Se  = zeros(1,length(y));
dd  = zeros(1,length(y));
theta  = zeros(1,length(y));

for n= 1: length(y);         % Find intersection with Poincare plane to 
  Pro(n)=NV'*y(:,n);         % define the beginning of the periods.
end;
IN=1;
for n = 2:length(Pro);
    if ((Pro(n)>0)&&(Pro(n-1)<0)) ; % Intersection with Poincare plane
        Se(n)=1;
        V(IN) = Pro(n)/ (Pro(n)-Pro(n-1));
        IN = IN+1;
    else
        Se(n)=0;
    end;
end;
dy=gradient(y); % Computing the covered distance in the state space
for n= 1: length(y);
    dd(n)=norm(dy(:,n));
end;

Dis = cumsum(dd);  % Covered distance along the trajectory

Pmin=find(Se==1); % Indices of the beginning of the cycles

for i= 1 : length(Pmin)-1;    % Computing protophase theta
    for j= Pmin(i) : Pmin(i+1)-1;
        R1 =  V(i)*(  Dis(Pmin(i)) - Dis(Pmin(i)-1)  );
        R2 =  (1-V(i+1)) * ( Dis(Pmin(i+1))-Dis(Pmin(i+1)-1) );
        theta(j)= 2*pi* ( Dis(j) - (Dis(Pmin(i))-R1) )   /   ( (Dis(Pmin(i+1)-1)+R2) - (Dis(Pmin(i))-R1)  ) ;
    end;
end;
theta=unwrap(theta);
Start = Pmin(1);
Stop =  Pmin(end)-1;
end