function R = profile_sim(sim)
% plot partial correlations between cell parameters and cell behavior


	plotme = 1;

   	T = sim;
   	tslice = [1:min(1000,sim.duration)];
	V = T.networkHistory.V_soma(:,tslice);

	minV = min(V,[],2);
	maxV = max(V,[],2);

	spikes  = spikedetect(sim);

	ampl = max(V, [], 2) - min(V, [], 2);
	meanVm = mean(V, 2);
	spks = spikes.spikespercell'/length(tslice)*1e3;

	midH = bsxfun(@minus, V, mean(V,2));
	pos = sum(sign(midH)>0, 2);
	neg = sum(sign(midH)<0, 2);
	supth = pos/length(tslice);

	K = measureGlobalSync(sim,'plotme', 0,'duration',tslice); %, 'duration', tslice
	freq_each = K.frequency'*1e3;

	freq = median(freq_each);
	

	noparams = length(fields(T.cellParameters));

	str = [];
	for ff = fields(T.cellParameters)'
		eval([ff{1} '= T.cellParameters.' ff{1} ';'])
				str = [str ff{1} ','];

	end

	eval(['R.allneurons = table(' str '  freq_each, ampl, meanVm, spks, supth, minV, maxV );'])


	[rho pval] = partialcorr(table2array(R.allneurons));

	R.partialcorr = rho;
	R.pval		  = pval;



if plotme
	try
		CB = flipud(cbrewer('div', 'RdBu',20));
		set(0, 'defaultfigurecolormap', CB)
	catch
		warning('no color brewer')
	end

	no_vars=  length(R.allneurons.Properties.VariableNames);

	figure 
		subplot(121)
			imagesc(rho,[-1 1]); colorbar
			set(gca,'xtick', [1:no_vars], 'xticklabel', R.allneurons.Properties.VariableNames)
			set(gca,'ytick', [1:no_vars], 'yticklabel', R.allneurons.Properties.VariableNames)
			title('partial correlation')

		subplot(122)
			imagesc(pval<0.05); colorbar

			set(gca,'xtick', [1:no_vars], 'xticklabel', R.allneurons.Properties.VariableNames)
			set(gca,'ytick', [1:no_vars], 'yticklabel', R.allneurons.Properties.VariableNames)
			title('p values')

	if isfield(sim,'Plist')
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
			legend(num2str(sim.Plist))

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