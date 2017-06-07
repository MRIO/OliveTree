
function out = ETR(varargin)
	% 
	% @OliveTree
% p.addRequired('triggers',@isvector) % triggers and spikes must be in ms
% p.addRequired('spikes',@isvector) % triggers and spikes must be in ms
%
% if SortISI = true, we sort by the first interval after the trigger.
% 
% sortISI can take arbitrary values (positive or negative) when 'spikewindow' is set to 'count'. 
% 
% p.addRequired('triggers',@isvector) % triggers and spikes must be in ms
% p.addRequired('spikes',@isvector) % triggers and spikes must be in ms
% p.addOptional('waves',[])
% p.addOptional('triggerlabels',[])
% p.addParamValue('spikewindow', 'time') % time or count
% p.addParamValue('spksLR', [20 20])
% p.addParamValue('sortISI', 1) % ISI to sort by (1 = first interval after the CS), 2 is the second, ....
% p.addParamValue('span', [-300 300]) 
% p.addParamValue('plotQ', true)
% p.addParamValue('Fs', 1000)
% p.addParamValue('fname', [])
% p.addParamValue('bin', 5)
% p.addParamValue('markercolor', [0 0 0])
% p.addParamValue('markersize', 5)
% p.addParamValue('markertype','.')
% p.addParamValue('markercolorEv', [.3 .3 1])
% p.addParamValue('axes', [])


p = inputParser;
p.addRequired('triggers',@isvector) % triggers and spikes must be in ms
p.addRequired('spikes',  @isvector)   % triggers and spikes must be in ms
p.addOptional('waves',[])
p.addOptional('triggerlabels',[])

p.addParamValue('spikewindow', 'time') % time or count (legacy)
p.addParamValue('spksLR', [20 20])
p.addParamValue('sortISI', 1) % ISI to sort by (1 = first interval after the CS), 2 is the second, ....
p.addParamValue('span', [500], @isscalar)
p.addParamValue('plotQ', true)
p.addParamValue('Fs', 1000)
p.addParamValue('fname', [])
p.addParamValue('bin', 5)
p.addParamValue('markercolor', [0 0 0])
p.addParamValue('markersize', 5)
p.addParamValue('markertype','+')
p.addParamValue('markercolorEv', [.3 .3 1])
p.addParamValue('axes', [])
normalizewaveamplitudes = 0;

p.parse(varargin{:});


triggers = p.Results.triggers;
spikes = p.Results.spikes;
LFP = p.Results.waves;
Fs = p.Results.Fs;
preSpikes = p.Results.spksLR(1);
postSpikes = p.Results.spksLR(2);
span = p.Results.span; 
sortISI = p.Results.sortISI;
plotQ = p.Results.plotQ;
triggerslabels = p.Results.triggerlabels;
fname = p.Results.fname;
bin = p.Results.bin;
spikewindow = p.Results.spikewindow;
markercolor = p.Results.markercolor;
markersize = p.Results.markersize;
markercolorEv = p.Results.markercolorEv;
markertype = p.Results.markertype;
axeshandle = p.Results.axes;

win = abs(round(span*Fs/1000)); % in samples

if size(triggers,2) ~= 1; triggers = triggers'; end
if size(spikes,2) ~= 1; spikes = spikes'; end
if size(LFP,2) ~= 1; LFP = LFP'; end


%%====================================================================     Tests

	% padding preemptively


if triggers-win<=0
	triggers(triggers-win<=0) = [];
end

if not(isempty(LFP))
	if sum(triggers>=length(LFP)-win)
		padsize =  triggers(end) - length(LFP) + win;
		LFP = padarray(LFP, padsize);
	else 
		padsize = 0;
	end
end



%%====================================================     Find spikes around triggers

switch spikewindow
	case 'count'
		ssAroundEvent = @(cs)([spikes(find(spikes<cs,preSpikes, 'last'))' cs spikes(find(spikes>cs,postSpikes,'first'))']) ;
		ssPreEvent = ...
		 @(cs)([spikes(find(spikes<cs,preSpikes, 'last'))']) ;
		ssPosEvent = ...
		 @(cs)([spikes(find(spikes>cs,postSpikes,'first'))']) ;
		spikesartriggers = [];
		for cscount = 1:length(triggers)
			ssPre = ssPreEvent(triggers(cscount));
			ssPos = ssPosEvent(triggers(cscount));
			noPre = length(ssPre);
			noPos = length(ssPos);
	
			try
				if  noPre~= preSpikes ;
					ssPre = padarray(ssPre,[0 preSpikes-noPre], NaN, 'pre');
				end
				if noPos ~= postSpikes ;
					ssPos = padarray(ssPos,[0 postSpikes-noPos], NaN, 'pos');
				end
				spikesartriggers(cscount, :) = [ssPre triggers(cscount) ssPos];
			catch E

				keyboard
			end
		end


		% checks
		try
		% center spikes
	
		centspikesartriggers = spikesartriggers - repmat(triggers,[ 1 preSpikes+postSpikes+1] );


		catch E
			keyboard
		end

		% sort by ISI
		if sortISI
			[OETR spikeOrder] = sortrows(centspikesartriggers,[preSpikes+1+sortISI preSpikes+2+sortISI]) ;
		else
			OETR = centspikesartriggers;
			spikeOrder = [1:length(triggers)];
		end




	case 'time'

		% SPKS = sort([spikes ; triggers]);
		SPKS = sort([spikes ]);
		
		spks  = @(cs)(SPKS(find(SPKS> (cs - win) & SPKS < (cs + win)) )- cs); 
		cspks = @(cs)(triggers(find(triggers> cs - win & triggers < cs + win) )- cs); 
		% remzero = @(c) c(find(c~=0));

		% spikes around triggers: spartrig
		spartrig = arrayfun(spks, triggers([1:length(triggers)]'), 'uniformoutput', false);                 
		% spartrig = cellfun(remzero, spartrig, 'uniformoutput', false);
		spartrigtriggers = arrayfun(cspks, triggers([1:length(triggers)]'), 'uniformoutput', false);
		spartrigmat = cell2mat(spartrig);                                                             

		getspkafint = @(l)(spartrig{l}(find(spartrig{l} > 0, 1,'first')));

		spkafint = arrayfun(getspkafint,[1:length(triggers)], 'uniformoutput' , 0);

		% check for pauses longer than window (empty cell in spkafint)
		emptycells = find(cellfun(@isempty, spkafint));
		
		if ~isempty(emptycells)
			% warning('found pauses larger than window in ETR' )
			for ec = emptycells;
				spkafint{ec} = win;
			end
		end
		spkafint = cell2mat(spkafint);
		
		if sortISI
			[OETR_ spikeOrder] = sort(spkafint) ;
		
			OETR = cell(length(triggers),1);
			OETER = cell(length(triggers),1);
			for cc = 1:length(triggers)
				OETR{cc} = spartrig{spikeOrder(cc)};
				OETER{cc} = spartrigtriggers{spikeOrder(cc)};
			end		
		else
			OETR = spartrig;
			OETER = spartrigtriggers;
			spikeOrder = [1:length(triggers)];
		end

end

%%======================================================================     ETWaves

if ~isempty(LFP)


	if normalizewaveamplitudes
			lfpAroundEvent = ...
		@(t)(LFP( round(t)-win:round(t)+win) /norm(LFP(round(t)-win:round(t)+win))) ;
	else
		lfpAroundEvent = ...
		@(t) LFP( round(t)-win:round(t)+win);
	end

		ETW = cell2mat(arrayfun(lfpAroundEvent, triggers'+padsize,'uniformoutput', false));

	OETW = ETW(:,spikeOrder);

end


%%=====================================================================     Plot


if ~isempty(cell2mat(OETR(:)))
		[H X] = ksdensity(cell2mat(OETR(:)), [-win:1:win],'kernel', 'box', 'width', bin);
	else
		disp('Ordered ETR variable is empty')
		return
end



if plotQ & isempty(triggerslabels)
	
		
	
	axcount = 1;
	if ~exist('markercolor'); markercolor  = [.7 0 0]; end
	linecolor =  [.6 .6 1];
	
	
	switch spikewindow
		case 'count'
	
		a(3) = axes;
		axis off
		% text(0,0,fname)
	
		a(1) = axes('color', 'none');
		if ~isempty(LFP)
			plottype = 'pcolor';
			switch plottype
				case 'pcolor'
					try
					pcolor([-win:win]/Fs*1000,[1:length(triggers)],OETW'); shading flat, colormap copper
					set(a(1),'position', [.13 .41 .80 .5])
					catch
						keyboard
					end

				case 'wavelines'
					waveLine = @(l)(line( [-win:win]/Fs*1000,  l+ETW(:,spikeOrder(l))/max(ETW(:,spikeOrder(l))), 'linestyle', 'none', 'marker', '.', 'color', 'k', 'markersize', markersize*5));
					arrayfun(waveLine, [1:length(triggers)])
				end
		end
	

		rasterLine = @(l)(line( OETR(l,:), l*ones(size(OETR(l,:))), 'linestyle', 'none', 'marker', '.', 'color', markercolor, 'markersize', markersize*10));
		arrayfun(rasterLine, [1:length(triggers)])
		xlabel(['Event triggered raster, sorted by ' num2str(sortISI) ' interval(s) after event.'])
		axis tight

		%%================================================================     histogram
	
		a(2) = axes;

		set(a(1), 'position', [.13 .11 .8 .7]);
		set(a(2), 'position', [.13 .81 .8 .13]);


		bar(X,sum(H'),1, 'edgecolor', 'k', 'facecolor', 'k'), hold on
		grid on, box off, axis off

		set(gcf,'color', ones(1,3));
		linkaxes(a, 'x');
		ylim([0 4*median(H(find(H)))]);
		xlim([-span span]);


		
	case 'time'
	
		if ~isempty(axeshandle) & isempty(LFP)
			if ~exist('axeshandle')
				axeshandle = axes;
			end
			axes(axeshandle)
			axis off
			p = get(axeshandle, 'position');
			a(1) = axes('position', [p(1) p(2)         p(3)  .7*p(4)]);
			a(2) = axes('position', [p(1) p(2)+.7*p(4) p(3) .15*p(4)]);
		elseif ~isempty(axeshandle)
			axes(axeshandle)
			axis off
			p = get(axeshandle, 'position');
			a(1) = axes('position', [p(1) p(2)         p(3)  .7*p(4)],'color', 'none');
			a(2) = axes('position', [p(1) p(2)+.7*p(4) p(3) .15*p(4)]);				
		else
			a(3) = axes;
			axis off
			a(1) = axes('position', [.13 .11 .8 .7]); ;
			a(2) = axes('position', [.13 .81 .8 .15]);;
		end
	
	
	%%==========================================================     Continuous Data
	
	if ~isempty(LFP)
		plottype = 'pcolor';
		switch plottype
			case 'pcolor'
				figure
				pcolor([-win:win]/Fs*1000,[1:length(triggers)],OETW'); shading flat, colormap copper


			case 'wavelines'
				axes(a(2))
				waveLine = @(l)(line( [-win:win]/Fs*1000,  l+ETW(:,spikeOrder(l))/max(ETW(:,spikeOrder(l))), 'linestyle', 'none', 'marker', '.', 'color', linecolor, 'markersize', 1));
				arrayfun(waveLine, [1:length(triggers)])
			end
	end
	
	
	%%===========================add=======================================     Rasters
	
		raster = OETR;
		rasterEv = OETER;
	
		rasterspikes = @(l)(line( raster{l}, l*ones(size(raster{l})), 'linestyle', 'none', 'marker', '.', 'color', markercolor, 'markersize', markersize*3));

		rastertriggers = @(l)(line( rasterEv{l}, l*ones(size(rasterEv{l})), 'linestyle', 'none', 'marker', markertype, 'color', markercolor, 'markersize', markersize*1.1));

		rastertriggers_core = @(l)(line( rasterEv{l}, l*ones(size(rasterEv{l})), 'linestyle', 'none', 'marker', markertype, 'color', markercolorEv, 'markersize', markersize));

		axes(a(1))
		arrayfun(rastertriggers, [1:length(triggers)])
		arrayfun(rasterspikes, [1:length(triggers)])
		arrayfun(rastertriggers_core, [1:length(triggers)])

		ylab = get(a(1),'ytick');
		set(a(1),'ytick',ylab(1:end-1))
		axis tight
		
		
		
		%%================================================================     histogram

		axes(a(2))
		% [H X] = hist(cell2mat(OETR(:)), [-win:bin:win]);
		colorscheme = 'blackonwhite';
		switch colorscheme
			case 'whiteonblack'
				area(X,H,1, 'edgecolor', [0 0 0], 'facecolor', [0 0 0])		

				box off
				set(a(2),'xtick', [])
				ylim([0 1.2*max(H(find(H)))]);
				set(gca,'color', ones(1,3));	

				line([-120 -20], [30 30], 'linewidth', 10,'color', 'k')		
		 		line([-120 -20], [30 30], 'linewidth', 8,'color', 'w')
				line([-70 -20], [20 20], 'linewidth', 10,'color', 'k')		
		 		line([-70 -20], [20 20], 'linewidth', 8,'color', 'w')
				line([-30 -20], [10 10], 'linewidth', 10,'color', 'k')		
		 		line([-30 -20], [10 10], 'linewidth', 8,'color', 'w')

				linkaxes(a, 'x');

				grid on
			case 'blackonwhite'
				bar(X,H,1, 'edgecolor', [0 0 0], 'facecolor', [0 0 0])		

				box off
				set(a(2),'xtick', [])
				ylim([0 1.2*max(H(find(H)))]);
				set(gca,'color', ones(1,3));	

				line([-120 -20], [30 30], 'linewidth', 10,'color', 'k')		
		 		line([-120 -20], [30 30], 'linewidth', 8,'color', 'w')
				line([-70 -20], [20 20], 'linewidth', 10,'color', 'k')		
		 		line([-70 -20], [20 20], 'linewidth', 8,'color', 'w')
				line([-30 -20], [10 10], 'linewidth', 10,'color', 'k')		
		 		line([-30 -20], [10 10], 'linewidth', 8,'color', 'w')

				linkaxes(a, 'x');

				grid on
		end
								
		
	end


	
end


out.eventTriggeredRaster = OETR;
out.triggers = triggers;
out.events = spikes;
out.eventOrder = spikeOrder;
out.Fs =  Fs;
if ~isempty(LFP)
out.eventTriggeredWaveforms = OETW;
end
out.CS = triggers(spikeOrder);
out.histogram = {X,H};

if exist('a')
	out.axeshandle = a;
end


function S = errorFunc(S, varargin)

