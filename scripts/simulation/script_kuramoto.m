% script_kuramoto.m

clear % lazy Mario ;D

pspaceoscillatornumber = 0;
	plotoscillatorexamples  = 1;
	num_oscillator = [2 5 10];
	nt = 5; % number of trials
	time  = 1 ; %s;
	cc = 10; % coupling coefficient

pspacefreqcoupling = 1; % for two coupled oscillators
	deltaf = 2;
	no_osc = 2;
	basefrequencies = [1 9];
	coupling_coefficients = [5 15];
	time = 5;


comparetwooscillators = 0; %slow vs fast frequencies
	time = 3; %s
	K = [10 10]*5;
	slowfreq = [.5 ; 9.5];
	fastfreq = [10.5 ; 19.5];

checkexpsync = 0;

% [=================================================================]
%  parameter space for different number of oscillator and 
% [=================================================================]

if pspaceoscillatornumber
	clear HH

	
	set(0,'defaultaxescolororder',linspecer(8))


	gs = 0;
	for no_osc = num_oscillator
		gs = gs+1;
		% basefreq = 3;
		% deltaf = [.45 .8 1];
		% oscillators = linspace(basefreq-deltaf(gs),basefreq+deltaf(gs),no_osc);

		for ii =1:nt
			seed = randi(10000);
			init_cond = linspace(0,2*pi,no_osc+1);
			init_cond = init_cond(1:end-1);
			K{ii} = kuramotoSheet([no_osc, 1],cc,'clusterize', [0 0 0 0],'time', time,'seed', seed,'dt', .001, 'plotme',0,'init_cond', init_cond ,'radius', no_osc, 'oscillators',[] );
			KK(ii,:) = K{ii}.orderparameter;
			OO(ii,:) = K{ii}.oscillators;
			SS(ii,:) = sum(K{ii}.state);

		end
	
		[HH(gs,:) XX] = hist(KK(:),[0:.1:1]);

		if plotoscillatorexamples
			figure
			ax(1) = subplot(2,1,1);
				[val ord] = sort([max(OO')'-min(OO')']);
				plot(KK(ord,:)')
				legend(num2str(val));
				% legend(num2str(std(OO')'));
				set(ax(1), 'position',[0.15    0.5824    0.55   0.3426])


			subplot(2,1,2)
				SS = SS(ord,:);
				SS = cell2mat(arrayfun(@(r) SS(r,:) + r*2*pi , [1:nt]' ,'uniformoutput',0));
				plot(SS','linewidth',2)
				hold on, plot(K{ii}.state')

			axes('Position', [0.7461    0.5824    0.16    0.3426])
				% [HH(gs,:) XX] = hist(KK(:),[0:.1:1]);
				
				% [HH(gs,:) XX] = hist(KK,[0:.1:1]);
				barh(XX,HH(gs,:))
				axis off
				axis tight
		end
		clear K
		clear KK
		clear OO
		clear SS



	end

	figure
	plot(XX,HH')
end

% [=================================================================]
%  compare two oscillators
% [=================================================================]

if comparetwooscillators

	%% for slower oscillators, there's more sync:
	no_osc =2;
	seed = 1
	init_cond = linspace(0,2*pi,no_osc+1);
	init_cond = init_cond(1:end-1);
	KKK = kuramotoSheet([no_osc, 1],K(1),'clusterize', [0 10 1 0],'time', time,'seed', seed,'dt', .001, 'plotme',0,'init_cond', init_cond ,'radius', no_osc,'oscillators', slowfreq*2*pi);

	figure
	plot(sum(KKK.state),'linewidth',2)
	hold on
	plot(KKK.state')



	no_osc =2;
	seed = 1
	init_cond = linspace(0,2*pi,no_osc+1);
	init_cond = init_cond(1:end-1);
	KKK = kuramotoSheet([no_osc, 1],K(2),'clusterize', [0 10 1 0],'time', time,'seed', seed,'dt', .001, 'plotme',0,'init_cond', init_cond ,'radius', no_osc,'oscillators', fastfreq*2*pi);


	figure
	plot(sum(KKK.state),'linewidth',2)
	hold on
	plot(KKK.state')


end

% [=================================================================]
%  parameter space  basefreq x coupling
% [=================================================================]


gs = 0;
if pspacefreqcoupling
	
	
	seed = 1;
	
	init_cond = linspace(0,2*pi,no_osc+1);
	init_cond = init_cond(1:end-1);

	gs = gs + 1;
	for cc = coupling_coefficients
		for ff = basefrequencies

		K{gs} = kuramotoSheet([no_osc, 1],cc,'clusterize', [0 10 1 0],'time', time,'seed', seed,'dt', .001, 'plotme',0,'init_cond', init_cond ,'radius', no_osc,'oscillators', [ff ; ff+deltaf]*2*pi);

			KK(gs,:) = 		K{gs}.orderparameter;
			OO(gs,:) = 		K{gs}.oscillators;
			SS(gs,:) =  sum(K{gs}.state);



		figure
			ax(1) = subplot(2,1,1);
			% [val ord] = sort([max(OO')'-min(OO')']);
			% plot(KK(ord,:)')
			plot(KK')
			% legend(num2str(val));
			set(ax(1), 'position',[0.15    0.5824    0.55   0.3426])


			ax(2) = subplot(2,1,2);
			% SS = SS(ord,:);
			% SS = cell2mat(arrayfun(@(r) SS(r,:) + r*2*pi , [1:nt]' ,'uniformoutput',0));
			plot(SS')
			hold on
			plot(K{gs}.state')


			ax(3) = axes('Position', [0.7461    0.5824    0.16    0.3426]);
				[HH(gs,:) XX] = hist(KK(:),[0:.1:1]);
				% [HH(gs,:) XX] = hist(KK,[0:.1:1]);
				barh(XX,HH(gs,:))
				axis off
				axis tight

		title(['coupling: ' num2str(cc) ' freq:' num2str(ff)])

		figure
			plot(abs(fft(SS', 1024)));
			title(['coupling: ' num2str(cc) ' freq:' num2str(ff)])

		% clear K
		clear KK
		clear OO
		clear SS



		end
	end





end



KKK = kuramotoSheet([10 10],100,'clusterize', [1 10 1 0],'time', 1,'seed', 10,'dt', .001, 'plotme',1, 'radius', 3);

%=============================exponential sync for 2 oscillators==============================%

% if checkexpsync
% it = 0;
% 	phases = linspace(0, pi, 10);
% 	for ph = 
% 		it = it+1;
% 		KKK = kuramotoSheet([1 2],20,'clusterize', [0 10 1 0],'time', 1,'seed', 10,'dt', .001, 'plotme',0, 'init_cond', [-pi ; ph] , 'oscillators', [10 ; 10]*2*pi);


% 		KK(it,:) = KKK.orderparameter;
% 	end

% 	plot(KK')

% 	xlabel('ms')
% 	title('time to sync')
% 	ylabel('order parameter')

% 	legend(num2str(phases'))

% end		




kuramotoSheet([10 10],200,'clusterize', [1 10 1 0],'time', 1,'seed', 0,'dt', .001, 'plotme',1,'radius', 3, 'oscillators',[] );