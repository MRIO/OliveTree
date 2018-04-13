% ## Copyright (C)  
% ##
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 2 of the License, or
% ## (at your option) any later version.
% ##
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ##
% ## You should have received a copy of the GNU General Public License
% ## along with this program; If not, see <http://www.gnu.org/licenses/>.

% ## Author:  mnegrello@gmail.com

clea
r
cell_function = 'vanilla'; % unless specified below
cell_function = 'devel'; % unless specified below


availablefieldstosave =  {'V_soma', 'I_cx36', 'Sodium_h', 'Potassium_n', 'Potassium_x_s', 'Calcium_k', 'Calcium_l', 'V_dend', 'Calcium_r',...
 'I_CaH', 'Potassium_s', 'Hcurrent_q', 'Ca2Plus', 'V_axon', 'Sodium_m_a', 'Sodium_h_a', 'Potassium_x_a', 'Ca2_soma'};
to_report =  {'V_soma','V_dend', 'Ca2Plus', 'I_cx36'};
to_report = availablefieldstosave;

pert = []; I_app = [];
dt = .025;


experiment = 'inhibition';
switch experiment

	case 'inhibition'

		dt = 0.01;
		simtime = 1000;
		netsize = [1 1 1];
		gap = 0;

		W = 0;
		cell_parameters = createDefaultNeurons(1);

		
		I_app = zeros(1, simtime*(1/dt));
		I_app(1,(300*(1/dt):400*(1/dt))) = -15;% uA/cm^2 
		gnoise = [0 0 0 0]; % 



end


% [transients] = IOnet('networksize', netsize,'perturbation', pert ,'appCurrent',I_app,'time',simtime,'g_CaL', g_CaL ,'W', W ,'ou_noise', gnoise ,'to_report', to_report,'gpu', gpu);
[transients] = IOnet('cell_function', cell_function, 'networksize', netsize, 'cell_parameters', cell_parameters,  ...
	'perturbation', pert ,'appCurrent',I_app,'time',simtime ,'W', W ,'ou_noise', gnoise ,'to_report', to_report ,'gpu',0 , 'delta', dt);


subplot(3,1,1)
plot(transients.networkHistory.V_soma)
ylabel('mV')
subplot(3,1,2)
plot(transients.networkHistory.V_soma'), hold on
plot(transients.networkHistory.V_dend','b')
ylabel('mV')
subplot(3,1,3)
plot(transients.networkHistory.Ca2Plus')
xlabel('ms')
ylabel('[Ca]')


replayResults(transients, 'plotallfields',1)


