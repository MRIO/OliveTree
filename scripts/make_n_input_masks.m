make_n_input_masks.m




% ampa random inputs
numberofmasks = 10; 
stim_interval = 100;
onset_of_stim = 300;
n_of_pulses = 20;
for nm = 1:numberofmasks

	pert.mask  	  {nm} = create_input_mask(netsize, 'all', 'synapseprobability', .25);
	pert.amplitude{nm} = 1;
	pert.triggers {nm} = onset_of_stim + cumsum(poissrnd(stim_interval,n_of_pulses,1)) ;
	pert.duration {nm} = 1;
	pert.type	  {nm} = 'ampa';
	
end