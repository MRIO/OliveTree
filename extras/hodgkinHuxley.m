% hodgkinHuxley.m

function Dydt = hodgkinHuxley(t,y)

% This function simulates the generation of action potentials in the brain
% by way of the Hodgkin Huxley model. 

% y(1) = Voltage
% y(2) = m Gating Variable
% y(3) = h Gating Variable
% y(4) = n Gating Variable

C   = 1;     % Capacitance of the cellular membrane

gNa = 120;   % Sodium channel conductance in mS/cm^2
             % If you ate poorly prepared Fugu (puffer blowfish), you would
             % consume vast amounts of tetrodoxin, blocking sodium channels
             % causing the conductance to go to zero (and causing you to die)
             
gK  = 36;    % Potassium channel conductance in mS/cm^2
             % If you were shot by a poisoned dart from the Amazon, it
             % would likely be laced with curare, which blocks potassium
             % channels, causing the conductance to go to zero (and causing
             % you to die, yet again)
             
gL  = .3;    % This is the conductance of the membrane, thankfully unblockable

vNa = 115;   % The resting potential of sodium (in mV). Due to differences in
             % the concentration of sodium inside and outside the cell.
             
vK  = -12;   % The resting potential of potassium (in mV). Due to differences
             % in the resting concentration of potassium inside and outside
             % the cell.
             
vL  = 10.6;  % The resting potential due to various other factors.

I   = 20;    % The stimulating current (from the synaptic input). Varying this
             % value leads to interesting effects. Too low, and action
             % potentials do not occur. Above about 7, and increases in
             % current will actually lead to firing of action potentials
             % more rapidly. There is a maximal limit, though!


             % [t y] = ode45(@hodgkinHuxley,[0 40],[0 .5 .5 .5]);
             
aM = .1*(25-y(1))./(exp((25-y(1))/10)-1);   % Basically curve fitting
aH = .07*exp(-y(1)/20);                     % Basically curve fitting
aN = .01*(10-y(1))./(exp((10-y(1))/10)-1);  % Basically curve fitting
bM  = 4*exp(-y(1)/18);                       % Basically curve fitting
bH  = 1/(exp((30-y(1))/10)+1);               % Basically curve fitting
bN  = .125*exp(-y(1)/80);                    % Basically curve fitting

Dydt(1,1) = (I - gNa*y(2)^3*y(3)*(y(1)-vNa) - gK*y(4)^4*(y(1)-vK) - gL*(y(1)-vL));
Dydt(2,1) = (1-y(2))*aM - y(2)*bM;
Dydt(3,1) = (1-y(3))*aH - y(3)*bH;
Dydt(4,1) = (1-y(4))*aN - y(4)*bN;




