% summa_xcorr.m
function out = xcorr_summa(varargin)

	p = inputParser;
	p.addRequired('sim')
	p.addParamValue('selectedneurons',[])
	p.addParamValue('nwins',2 )
	p.addParamValue('plotme',1 )

	p.parse(varargin{:});

	simulation = p.Results.sim;
	selectedneurons = p.Results.selectedneurons;
	nwins = p.Results.nwins;
	plotthem = p.Results.plotme;

	multiwindowxcorr = 0;
	applyconvkernel =1;	
	sortbyasym = 1;

	flipasym = 0;

	lag = 400;
	centerwin = 50;
	sidewin = [51:200];
	noneur = 20;
	nselneu = length(selectedneurons);


	try
		duration = simulation.duration;
	catch
		duration = simulation.simulationParameters.duration;
	end

	winsize = floor(duration/(nwins));
	xcorrwin = 1:duration;

	XnoAC  = [];
	
	VSoma = simulation.networkHistory.V_soma;
	
	
	% [================================================]
	%  start
	% [================================================]
	
	



if applyconvkernel
	% CS Kernel
	N = 51;
	mog_filter_2d_CS = zeros(1, N);
	x_2d        = (fix(-N/2):fix(N/2));
	for n = [4 7 10] 
	    mog_filter_2d_CS = mog_filter_2d_CS + exp(-x_2d.^2/(2*n^2));
	end
	mog_filter_2d_CS = mog_filter_2d_CS./sum(mog_filter_2d_CS);

	
end


		

% load /Users/M/Cabinet/SyncBox/Bench/stim5pa_1hz/sim3D_long_Xcorr_onlyexcitation_1.mat
% load /Users/M/Public/BitSync/Bench2/sim3D_long_Xcorr_onlyexcitation_5.mat
% load /Users/M/Public/BitSync/Bench/sim3D_long_Xcorr_spontaneous_5.mat
% load /Users/M/Public/BitSync/Bench2/simresults_ampa30stim.mat


		W = simulation.networkParameters.connectivityMatrix;
		netsize = simulation.networksize;
		depth   = [1:netsize(1)];
        breadth = [1:netsize(2)];
        height  = [1:netsize(3)];


		gapneighborhood = full(sum(triu(W)))+eps;

		
		W  = simulation.networkParameters.connectivityMatrix;
		

%=============================select neurons==============================%
			
			if isfield(simulation, 'spikes')
				spks = simulation.spikes;
				VSB  = spks.binaryspikes;
			else
				spks = spikedetect(simulation);
				VSB = spks.binaryspikes;
			end

			
			takemostactive = 0;
			takeall = 1;
			if not(isempty(selectedneurons))
				% do nothign because
				% selectedneurons = selectedneurons;

			elseif takemostactive

				% take N most active neurons
				[V sorted] = sort(sum(VSB(:,xcorrwin)'),2,'descend');

				selectedneurons = sorted(1:noneur);

			else

				selectedneurons = find(sum(VSB(:,xcorrwin)')>5);

				if length(selectedneurons) <= noneur
					disp('not enough firing neurons')
					keyboard
					length(selectedneurons)
					out = [];
					return
				else
					selectedneurons = selectedneurons(randperm(length(selectedneurons)));
					selectedneurons = selectedneurons(1:noneur);
				end

				
			end

			N = length(selectedneurons);


%=============================compute xcorrs==============================%

pairs = find(triu(ones(N) - eye(N))') ;
		
		xcorrwin = round(linspace(1,duration,nwins+1));


		for nw = 1:length(xcorrwin)-1;


			XC{nw} = xcorr(VSB( selectedneurons ,xcorrwin(nw):xcorrwin(nw+1))',lag,'coeff');

			xcorrset = XC{nw}(:,pairs)';

			if applyconvkernel
				XCs{nw} = conv2(1, mog_filter_2d_CS,xcorrset,'same');
			end

			
			if flipasym
				asym = abs(sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)'));
			else
				asym = sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)');
			end
			
			XnoAC(:,nw) = mean(XCs{nw},'omitnan');


			[ampl(nw,:) delay(nw,:)] = max(XCs{nw}(:,lag-centerwin:lag+centerwin+1), [], 2);

			if flipasym
				delay = abs(delay);
			end


			%=============================title==============================%

			% if flipasym
			% 	sdpk = abs(sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)'));
			% else
			% 	sdpk = sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)');
			% end



			if isfield(simulation,'noiseapplied')
				XCNoise{nw} = xcorr(simulation.noiseapplied(selectedneurons , xcorrwin(nw)-1:1/dt:xcorrwin(nw+1))',lag,'coeff');
			end



		end


	delay = delay - centerwin +1 ;


vsomacorr = 0;
if vsomacorr
	% VS = simulation.networkHistory.V_soma;
end



out.winsize = winsize;
out.asymmetry = asym;
out.xcorr_pairs = XCs;
out.spikespercell = sum(VSB');
out.pairs = pairs;
out.delay = delay;
out.amplitude = ampl;
out.XcorrNoAc = XnoAC;
out.XC = XCs;
out.selectedneurons = selectedneurons;

if isfield(simulation,'noiseapplied')
	out.XCNoise = XCNoise;
end


if plotthem

	if 0
		figure
	        [X Y Z] = meshgrid(depth,breadth,height);
	        X = X(:); Y = Y(:); Z = Z(:);

	        idx = zeros(length(X),1);
	        idx(selectedneurons) = 1;

	        plotnetstruct(W,X,Y,Z,sum(VSB'))
	end

	if 0
		figure
			if sortbyasym
				[V O] = sort(asym);
				imagesc([-lag:lag], [1:length(pairs)], xcorrset(O,:))
			else
				imagesc([-lag:lag], [1:length(pairs)], xcorrset)
			end

			line([0 0], [0 length(pairs)],'color', 'w')
			% title({num2str(Ptable(simcount,:)) ; 'sc gp nc na fr'})
			ylabel('pairs')
			xlabel('lag')
	end
		
		% Ptable = [sc gp nc na fr ];


	figure
		highlighted = pairs(1);
		corder = cbrewer('seq', 'Greys',50);
		set(0,'defaultaxescolororder', flip(corder))
		axi(1) = subplot(5,1,1);
		plot([-lag: lag], XnoAC)
		hold on
		plot([-lag: lag], XCs{1}(highlighted,:),'r')
		xlabel('lag(ms)')
		
		% legend(num2str(sims2p'))
		% legend({'WT' ; 'Cx36' })
		axis tight
		xlabel('ms')
		ylabel('correlation (coeff)')
		xlabel('aggregate correlation (windowed)')	


		axi(2) = subplot(5,1,2);
		% plot([-lag: lag],XCs{1} )
		area([-lag: lag], XCs{1}')
		xlabel('lag(ms)')
		ylabel('aggregate coeff')


		axis tight
		xlabel('ms')
		ylabel('correlation (coeff)')
		xlabel('individual cross-correlations (sample)')	


		axi(3) = subplot(5,1,3);
			imagesc( -lag:lag, 1:length(pairs) , XCs{1})
			line([-lag -lag+10], [highlighted *ones(1,2)],'color', 'r')
			xlabel('lag(ms)')
			ylabel('pairs')

			
		axi(4) = subplot(5,1,4);
			% 
			if nwins == 2
				line([zeros(size(delay(1,:))) ; ones(size(delay(1,:)))  ], [delay(1,:); delay(2,:)])
			else
				[v_  sortord] = sort(XCs{1}(:,lag+1));
				imagesc([-lag: lag], 1:length(pairs), XCs{1}(sortord,:))
			end

		axi_ = subplot(5,1,5);
		
		
			plot(0:1000, VSoma(selectedneurons,5500:6500))
			hold on, 
			plot(0:1000, VSoma(selectedneurons(highlighted),5500:6500),'r','linewidth',2)
			colormap(flip(corder))

		linkaxes(axi,'x')



			if 0
			figure
				subplot(131)
				boxplot([ampl(1,:)'],'notch', 'on'); title('amplitude'); 
				subplot(132)
				boxplot([ delay(1,:)' ],'notch', 'on'); title('delay')
				subplot(133)
				boxplot([ asym(1,:)'],'notch', 'on'); title('asymmetry')
			end	


end