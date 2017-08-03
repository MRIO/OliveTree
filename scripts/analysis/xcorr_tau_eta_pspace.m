
files = {...
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_100.2_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_100.4_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_100.6_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_100_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_200.2_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_200.4_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_200.6_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_200_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_300.2_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_300.4_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_300.6_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_300_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_400.2_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_400.4_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_400.6_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_400_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_500.2_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_500.4_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_500.6_1_iso_spont_50000_1_.mat'
'/home/titanuser1/Sync/Titan/Bench2/periodic_ampa_tau_eta_pspace_4xcorr/periodic_ampa__tau_eta_pspace_tau_eta_500_1_iso_spont_50000_1_.mat'};



for f = 1:length(files)
	load(files{f})
	S{f} = simresults{1};
	R{f} = profile_sim(simresults{1})
	X{f} = xcorr_summa(S{f}, 'plotme',1);
	freq(f) = sum(X{f}.spikespercell/(50*200));
	params(f,:) = [1/S{f}.perturbation.noise(1) S{f}.simulationParameters.sametoall];
	% saveallfigs('prefix', ['xc_pspace' num2str(params(f,:))]  )
	close all

end



for f = 1:length(files)
	
	XCs(f,:) = X{f}.XcorrNoAc;
	
end

[params_ord ord  ] = sortrows(params,[2 1]);

XCs_ord = XCs(ord,:);



cb = cbrewer('seq', 'Blues', 8);
corder(1:5,:)   = cb(4:8,:);
cb = cbrewer('seq', 'Greens', 8);
corder(6:10,:)  = cb(4:8,:);
cb = cbrewer('seq', 'Oranges', 8);
corder(11:15,:) = cb(4:8,:);
cb = cbrewer('seq', 'Purples', 8);
corder(16:20,:) = cb(4:8,:);

figure
set(0,'defaultAxesColorOrder', corder)

plot([-400:400],XCs_ord')

legen = [params_ord freq(ord)']

legend(num2str(legen))



figure
plot([-400:400],bsxfun(@rdivide, XCs_ord, XCs_ord(:,401)))


figure
set(0,'defaultAxesColorOrder', corder(5:5:end,:))
plot([10:10:50], reshape(freq(ord),5,4))