function [theta, START, STOP]=co_mmzproto(x, pl)
% DAMOCO Toolbox, function CO_MMZPROTO, version 06.05.14
%
% Given a time series x, co_mmzproto computes a protophase 
% by means of linear interpolation between marker events. 
% It assumes, that the time series has one maximum, one minimum, 
% and two zero-crossing in each cycle, and uses these 4 events 
% as markers to compute the protophase.
% Since the algorithm uses only 4 markers per cycle, it grasps 
% only slow components of the protophase dynamics. 
% If omega0 is the frequency of the oscillation, then, roughly, 
% the components of the phase dynamics for omega > 2*omega0 are 
% filtered out by simplifying the signal to 4 markers per cycle.
%
%                   WARNING! 
% Do not use the algorithm blindly: it possibly has to be adjusted
% to deal with your particular time series! 
%
%       INPUT
%                       x:       given time series
%                       pl:      if pl ==1 the (1) time series and the
%                                marker-events and (2) the time series as 
%                                function of theta are plotted to check 
%                                whether the markers are chosen properly. 
%                                If the choice is reasonable, plot (2) 
%                                shows a softly smeared periodic pattern. 
%                                Default value: pl=0, no plot.
%
%       OUTPUT
%                   theta:       protophase     
%                   START:       The protophase is defined for complete
%                                cycles only. START is the index of the 
%                                beginning of the first cycle in the time
%                                series x.
%                   STOP:        STOP is the index of the end of the last
%                                complete cycle in the time series x.
%                                Thus,theta is defined for x(START:STOP).
%
if nargin <2; pl=0; end;
S=size(x);if S(1) > S(2); x=x'; end;
gx = diff(x); % gradient(x);
s = sign(gx(1:end-1))-sign(gx(2:end));  % Marking changes of the derivative: s=0: no change, 
%                                         s=2: positive to negative (maximum), 
%                                         s=-2: negative to positive (minimum)
IN1 = find(s ==2 )+1;                     % Index of maxima
IN1c = IN1 + gx(IN1-1) ./ (gx(IN1-1)-gx(IN1)); % Correction of the position of the maxima
IN3 = find(s == -2)+1;                    % Index of minima
IN3c = IN3 + gx(IN3-1) ./ (gx(IN3-1)-gx(IN3)); % Correction of the position of the minima
IN2 = find(x(1:end-1)>0 & x(2:end)<0);  % Index of the zero crossings: positive to negative  
IN2c = IN2 + x(IN2) ./ ( x(IN2) - x(IN2+1)); % Correction of the position of the zero-crossings
IN4 = find(x(1:end-1)<0 & x(2:end)>0);  % Index of the zero crossings: negative to positive  
IN4c = IN4 + x(IN4) ./ ( x(IN4) - x(IN4+1)); % Correction of the position of the zero-crossings

if IN2c(1)<IN1c(1); IN2c(1)=[]; end;     % Cycle starts with first maximum. 
                                         % All others markers before the first maximum are eliminated
if IN2c(end)>IN1c(end); IN2c(end)=[]; end;
if IN3c(1)<IN1c(1); IN3c(1)=[]; end;
if IN3c(end)>IN1c(end); IN3c(end)=[]; end;
if IN4c(1)<IN1c(1); IN4c(1)=[]; end;
if IN4c(end)>IN1c(end); IN4c(end)=[]; end;

%%%%%%%%%%%%%%%%%%  CHECK FOR CONSISTENCE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% The detection of the markers is based on finding of local minima / maxima: 
% they also can be due to noise or other artifacts. To avoid  
% errors die to artifacts the algorithm checks for consistence of the number 
% of the 4 types of marker events. In case of deviating numbers of marker events
% it displays a warning and plots the difference of the index of the
% marker events, diff(IN). A sudden jump of diff(IN) indicates the location
% of the artifact. In case of inconsistence, smoothing the time series typically helps.
%
% Check for consistence of the number of marker events 
if length(IN2c)~=length(IN3c)|length(IN3c)~=length(IN4c)|length(IN2c)~=length(IN4c)|length(IN1c)-1~=length(IN2c); 
    display('WARNING: number of marker events inconsistent');
    theta=NaN; START=NaN; STOP=NaN;
    figure;                         % Plots to locate the source of the inconsistence;
    plot(diff(IN1c)); hold on;      % sudden 'jumps' of diff(IN) indicate irregular events 
    plot(diff(IN2c),'r');
    plot(diff(IN3c),'g');
    plot(diff(IN4c),'k');
    hold off
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:length(IN1c)-1
    D(n) = IN1c(n+1)-IN1c(n);            % Length of each cycle
    R2(n) = (IN2c(n)-IN1c(n)) / D(n); % Position of the 1st zero crossing relative to the cycle length
    R3(n) = (IN3c(n)-IN1c(n)) / D(n); % Position of the minimum relatively to the cycle length
    R4(n) = (IN4c(n)-IN1c(n)) / D(n); % Position of the 2nd zero crossing relative to the cycle length
end
R2 = mean(R2); % Average position of the 1st zero crossing
R3 = mean(R3); % Average position of the minima
R4 = mean(R4); % Average position of the 2nd zero crossings
IN = [];
Pin = [];
for n= 1: length(IN1c)-1;
    IN =[IN IN1c(n) IN2c(n) IN3c(n) IN4c(n)]; % Ordering the positions of markers for interpolation
    Pin = [Pin 2*pi*(n-1)  2*pi*(n-1+R2) 2*pi*(n-1+R3) 2*pi*(n-1+R4)];  % Ordering the values of the 
end                                                 %  Protophase of the corresponding markers
IN=[IN IN1c(end)];
Pin=[Pin 2*pi*(length(IN1c)-1)];
theta= interp1(IN, Pin, (IN1(1)+1:1:IN1(end)),'linear'); % Computing the protophase by 
                                                          % linear interpolation between markers.
START=IN1(1)+1; STOP=IN1(end);
if pl==1;
    figure;
        plot(mod(theta,2*pi),x(START:STOP),'.','markersize',4)
    figure;
        plot(x);
    hold on; 
        plot(IN1,x(IN1),'.r')
        plot(IN2c, 0,'.g')
        plot(IN3,x(IN3),'.k')
        plot(IN4c, 0,'.m')
    hold off
  figure;
  plot(mod(theta,2*pi),gradient(theta),'.','markersize',4)
end
end
end
