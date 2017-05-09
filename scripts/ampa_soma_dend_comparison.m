% ampa_soma_dend_comparison.m


gap = 0.04;  noisesig = 0; noiseamp = 0 ; tau = 30; sametoall = 0.15; spont = 0; conntype = 'iso' ;  gapcomp = 0;
simtime = 800;
netsize = [1 1 1]

M1 = zeros(prod(netsize),1); M1(1) = 1;
M2 = zeros(prod(netsize),1); M2(3) = 1;

	pert.mask{1}  	= M1;
	% pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa_dend';
	pert.duration{1}  = 1;
	pert.triggers{1} = [501:4:512  ];

	singlesim

	vdend{1} = simresults.networkHistory.V_dend;
	vsoma{1} = simresults.networkHistory.V_soma;


	pert = [];

	pert.mask{2}  	  = M1;
	% pert.amplitude{2} = 2;
	pert.type{2}	  = 'ampa';
	pert.duration{2}  = 1;
	pert.triggers{2} = [501:4:512];

	singlesim

	vdend{2} = simresults.networkHistory.V_dend;
	vsoma{2} = simresults.networkHistory.V_soma;


figure
plot(vsoma{1}','b')
hold on
plot(vsoma{2}','r')
ylim([-70 -40])
legend({'dendrite' 'soma'})
title('Vm @ soma')

figure
plot(vdend{1}','b')
hold on
plot(vdend{2}','r')
title('Vm @ dendrite')
ylim([-70 -40])
legend({'dendrite' 'soma'})