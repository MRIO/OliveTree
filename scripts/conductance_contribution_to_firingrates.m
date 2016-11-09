function D = conductance_contribution_to_firingrates(sim)


	% eval(['load /Users/M/Cabinet/SyncBox/Bench/stim5pa_1hz/sim3D_long_Xcorr_onlyexcitation_' num2str(s) '.mat']);

	netsize = sim.networksize;

	if isfield('triggers', sim.perturbation)
		trigger = 2;
		triggers = sim.perturbation.triggers{trigger};
	end
	spks = spikedetect(sim ,0,0);

	W = sim.networkParameters.connectivityMatrix;

	gapneighborhood = full(sum(triu(W)))+eps;

	meanISI = cell2mat(cellfun(@mean, cellfun(@diff,spks.spikes,'uniformoutput', false),'uniformoutput', 0));
	meanISI(isnan(meanISI))= 0;

	D = dataset(ones(prod(netsize),1), gapneighborhood', sim.cellParameters.g_CaL, sim.cellParameters.g_h, sim.cellParameters.g_ls, spks.spikespercell', meanISI', sim.simulationParameters.noiseApplied);
	D.Properties.VarNames = {'sim', 'gap', 'caL' ,'h' ,'ls', 'no_spks', 'meanISI', 'noise' }




az = 17;
el = 22;



if 1



	subplot(3,2,1)
	scatter3( sim.cellParameters.g_CaL, gapneighborhood, spks.spikespercell+1, 25,sim.cellParameters.g_h,'filled'); 
	% l1 = line( sim.cellParameters.g_CaL, gapneighborhood, spks.spikespercell, 'marker', '+','color', 'r', 'markersize',50,'linewidth',4);
	xlabel('Ttype conductance'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)

	subplot(3,2,2)
	scatter3( sim.cellParameters.g_h, gapneighborhood, spks.spikespercell+1, 25, sim.cellParameters.g_CaL,'filled'); view([az el])
	xlabel('g_h'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)

	subplot(3,2,3)
	scatter3( sim.cellParameters.g_int, gapneighborhood, spks.spikespercell+1, 25, sim.cellParameters.g_CaL,'filled'); view([az el])
	xlabel('g_int'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)

	subplot(3,2,4)
	scatter3( sim.cellParameters.g_ls, gapneighborhood, spks.spikespercell+1, 25, sim.cellParameters.g_CaL,'filled'); view([az el])
	xlabel('g_l_s'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)


	subplot(3,2,5)
	scatter3( sim.cellParameters.g_K_Ca, gapneighborhood, spks.spikespercell+1, 25, sim.cellParameters.g_CaL,'filled'); view([az el])
	xlabel('g_K_Ca'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)

	subplot(3,2,6)
	scatter3( sim.simulationParameters.noiseApplied, gapneighborhood, spks.spikespercell+1, 25, sim.cellParameters.g_CaL,'filled'); view([az el])
	xlabel('noise applied'), ylabel( 'gapneighborhood'),zlabel( 'spikes frequency')
	view(az,el)







% 	figure



% 	cla
% 	scatter3(coord1,coord2, coord3,100,spks.spikespercell+1,'filled'); view([az el]); axis equal
% 	l2 = line(coord1(i), coord2(i), coord3(i),'marker', '.','color', 'r', 'markersize',70);
% 	axis equal
% 	view(-25,28)
% 	title('spikes per cell')

end