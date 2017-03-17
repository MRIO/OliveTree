function R = profile_sim(varargin)
% plot partial correlations between cell parameters and cell behavior




	p = inputParser;
	p.addRequired('sim')  
	% p.addRequired('input')  
	% p.addRequired('input')  
	% p.addRequired('input')  


	p.addParamValue('tslice', []) 
	p.addParamValue('plotme', 0)
	% p.addParamValue('input', default) 
	% p.addParamValue('input', default) 
	% p.addParamValue('input', default) 

	p.parse(varargin{:});

	sim = p.Results.sim;
	% input = p.Results.input;
	% input = p.Results.input;
	% input = p.Results.input;
	tslice = p.Results.tslice;
	plotme = p.Results.plotme;
	% input = p.Results.input;
	% input = p.Results.input;
	% input = p.Results.input;




   	T = sim;
   	noneurons = prod(sim.networksize);
   	
   	if isempty(tslice); tslice = [1:min(1000,sim.duration)]; end

	V = double(T.networkHistory.V_soma(:,tslice));

	if length(tslice) < 400
		disp('profile_sim requires at least 400ms of network activity')
		return
	end

	minV = min(V,[],2);
	maxV = max(V,[],2);

	spikes  = spikedetect(sim);

	% ampl = max(V, [], 2) - min(V, [], 2);

	for ii = 1:noneurons
		
		ampl(1,ii) = mean(findpeaks(V(ii,:),'minpeakdist',50)) -(-(mean(findpeaks(-V(ii,:),'minpeakdist',50))));
	

	end

	ampl(1,isnan(ampl))=0;

	ampl = ampl';

	meanVm = mean(V, 2);
	spks = spikes.spikespercell'/length(tslice)*1e3;

	midH = bsxfun(@minus, V, [ (max(V')-min(V'))/2+ min(V') ]');
	pos = sum(midH>1, 2);
	neg = sum(midH<-1, 2);
	supth = pos/length(tslice);


	K = measureGlobalSync(sim,'plotme', 0,'duration',tslice); %, 'duration', tslice

	freq_each = K.frequency'*1e3;

	freq = median(freq_each);

	noparams = length(fields(T.cellParameters));

	str = [];
	bla = 0;
	for ff = fields(T.cellParameters)'
		if strcmp(ff{1}, 'Plist') | strcmp(ff{1}, 'Pnames'); 
			continue
		end
		eval([ff{1} '= T.cellParameters.' ff{1} ';'])
				str = [str ff{1} ','];
	end


	eval(['R.allneurons = table(' str '  freq_each, ampl, meanVm, spks, supth, minV, maxV );'])

	% ColumnsOfInterst = T.Pnames; 
	
	if isfield(T, 'Pnames')
		Pnames = T.Pnames';
		ColumnsOfInterst = [Pnames  'freq_each', 'ampl', 'meanVm', 'spks', 'supth' ];
		[rho pval] = partialcorr(table2array(R.allneurons(:,ColumnsOfInterst)),'rows', 'complete');
	else
		ColumnsOfInterst = R.allneurons.Properties.VariableNames;
		[rho pval] = partialcorr(table2array(R.allneurons),'rows', 'complete');
	end

	% 

	% [rho pval] = partialcorr(table2array(R.allneurons,'rows', 'complete');


	R.partialcorr = rho;
	R.pval		  = pval;



if plotme
	try
		CB = flipud(cbrewer('div', 'RdBu',20));
		set(0, 'defaultfigurecolormap', CB)
	catch
		warning('no color brewer')
	end

	no_vars=  length(ColumnsOfInterst);

	figure 
		subplot(121)
			imagesc(rho,[-1 1]); colorbar
			set(gca,'xtick', [1:no_vars], 'xticklabel', ColumnsOfInterst)
			set(gca,'ytick', [1:no_vars], 'yticklabel', ColumnsOfInterst)
			title('partial correlation')

		subplot(122)
			imagesc(pval<0.05); colorbar

			set(gca,'xtick', [1:no_vars], 'xticklabel', ColumnsOfInterst)
			set(gca,'ytick', [1:no_vars], 'yticklabel', ColumnsOfInterst)
			title('p<0.05?')

	if isfield(sim,'Plist')
		sim.Plist = round(sim.Plist*100)/100;
		figure
			ca = axis;
			set(0,'defaultaxescolororder', linspecer(length(sim.Plist)))
			% p = plot(tslice,   sim.networkHistory.V_soma');
			p = plot(  sim.networkHistory.V_soma');
			legend(num2str(sim.Plist))

		try; maximize_fig; catch; end

		figure
			imagesc(sim.networkHistory.V_soma,[-80 -20]), colorbar
			set(gca,'ytick', [1:length(sim.Plist)],'yticklabel', num2str(sim.Plist),'fontsize',8)
			try
			tit = sim.Pnames';
			title(tit)
		catch
			disp('bla')
		end

		try; maximize_fig; catch; end

	end

end








% from mathworks

% X = [ones(size(x1)) x1 x2 x1.*x2];
% b = regress(y,X)    % Removes NaN data

% scatter3(x1,x2,y,'filled')
% hold on
% x1fit = min(x1):100:max(x1);
% x2fit = min(x2):10:max(x2);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
% mesh(X1FIT,X2FIT,YFIT)



% scatter3(x1,x2,y,'filled')
% hold on
% x1fit = min(x1):100:max(x1);
% x2fit = min(x2):10:max(x2);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
% mesh(X1FIT,X2FIT,YFIT)