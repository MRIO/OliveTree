% run_sanity_checks.m

experiments = {'one_cell'; 'tycho'; 'three_cells'; '2Dcenter'; '3Dcenter'; 'gaba_soma'; 'gaba_soma_dend'; 'two_cells';'continuations'; 'somatic_ampa'};


for e = 5:numel(experiments)

	experiment = experiments{e};
		
	try
		sanity_check
	catch E
		E.stack
		E.message

		ERROR{e} = E;

	end
end
