% #jochen


thisprefix = 'iso_sigfiles_';
% thisprefix = 'cluster_vs_iso';

	selected_neurons_sel{1} = 'neighbors';	
	selected_neurons_sel{2} = 'nextneighbors'; 
	selected_neurons_sel{3} = 'stimulated';
	selected_neurons_sel{4} = 'highgapconn';
	selected_neurons_sel{5} = 'leastgapconn';
	selected_neurons_sel{6} = 'acrosscluster';
	selected_neurons_sel{7} = 'all';

	clear T

	for fff = [1 2 3 7]
		selected_neurons = selected_neurons_sel{fff}
		
		xcorr_feature_analysis;
		disp('step' )

		% keyboard
		
		pth = [ thisprefix selected_neurons_sel{fff} '/'] 
		eval(['mkdir ' pth])

		saveallfigs('prefix', [pth '2clust_', selected_neurons])
		disp('saved figs')

		T(fff,:) = [length(selectedneurons) sum(cell2mat(sigs)) cell2mat(prop_sig_pairs)];
		
		eval([selected_neurons '_pack= {collectX ; collectAmpl ;collectDelay ;collectAsym ; sig_pairs; prop_sig_pairs }']) ;
		% save([pth selected_neurons '_pack' ])





		

		close all

	end		

T = array2table(T);
 T.Properties.RowNames = selected_neurons_sel;
 T.Properties.VariableNames = {'neurons' 'WT1' 'MT1' 'WTspont' 'MTspont'  'WT1prop' 'MT1prop' 'WTspontprop' 'MTspontprop'}

% 	numsel	 WT1H	   MT1 		WT_spot	  MT_spot	   WT1H	   MT1 		WT_spot	  MT_spot
%     9.0000   36.0000   31.0000   36.0000   31.0000    1.0000    0.8611    1.0000    0.8611
%    30.0000  349.0000  275.0000  354.0000  276.0000    0.8023    0.6322    0.8138    0.6345
%    30.0000  405.0000  404.0000  246.0000  253.0000    0.9310    0.9287    0.5655    0.5816
%     5.0000   10.0000    6.0000   10.0000    6.0000    1.0000    0.6000    1.0000    0.6000
%     3.0000    3.0000    1.0000    3.0000    1.0000    1.0000    0.3333    1.0000    0.3333
%    30.0000  281.0000  205.0000  275.0000  204.0000    0.6460    0.4713    0.6322    0.4690
%    30.0000  282.0000  232.0000  284.0000  228.0000    0.6483    0.5333    0.6529    0.5241



% T1 
%  1.0000    0.8611    1.0000    0.8611
%  0.8023    0.6322    0.8138    0.6345
%  0.9310    0.9287    0.5655    0.5816
%  1.0000    0.6000    1.0000    0.6000
%  1.0000    0.3333    1.0000    0.3333
%  0.6460    0.4713    0.6322    0.4690
%  0.6483    0.5333    0.6529    0.5241
