
function varargout = NDscatter(varargin)



p = inputParser;
p.addRequired('intable')
p.addOptional('groupby',[]) 
p.addParamValue('fields',[])

p.parse(varargin{:}) 

intable = p.Results.intable;
groupby = p.Results.groupby;
selectedfields  = p.Results.fields;
plotdistributions = 0;

VarNames = intable.Properties.VariableNames;
data = table2array(intable);
N = length(VarNames);

figure

if ~isempty(groupby)
	if length(groupby)>1
		groups = groupby;
	else
		[C_ r_ groups] = unique(data(:,groupby));
	end
else
		groups = ones(size(data,2));
end
	
	
	
	bgcolor = [1 1 1];
	axiscolor = [.1 .1 .1];
	markercolor = [0 0 0];
	
	ngroups = max(groups);
	try
	groupcolors = linspecer(ngroups);
	catch
		groupcolors  = [0 0 0];
	end



	l1 = linspace(0, .88, N+1); l1 = l1(1:end -1)+.08; 
	l2 = linspace(0, .88, N+1); l2 = l2(1:end -1)+.08; l2 = fliplr(l2);
	w  = .80/N;
	h  = .80/N;

	r = 0; c = 0;
	for rr = l1
		r = r+1;
		for cc = l2
			c = c+1;
			
			hand(r, c) = axes('position', [cc rr w h]);
			if r > c
				hold on
				for g = 1:ngroups

					line(data(groups==g,c), data(groups==g,r), ...
						 'markersize', 5, 'color', groupcolors(g,:), 'marker','.','linestyle','none','markersize',15)
				end
				line(data(:,c), data(:,r) , 'color', markercolor, 'marker','.','linestyle','none','markersize', 5)
				hold off
			elseif r == c
				if min(data(:,r)) == max(data(:,r))
					continue
				end
				xxx = linspace(min(data(:,r)), max(data(:,r)), 30);

				if ngroups>1 & ngroups <3
					for g = 1:ngroups
						[hhh] = hist(data(groups==g,r),xxx);
						plot(xxx, hhh,'color', groupcolors(g,:) )
						hold on

					end
				else
					
					[hhh xxx] = hist(data(:,r), xxx);
					bar(xxx, hhh, 'facecolor', markercolor)
				end

				set(gca,'ytick', [0 max(hhh)])

			end

			if c > r % & plotdistributions
				edges = {linspace(min(data(:,c)), max(data(:,c)) , 20); linspace(min(data(:,r)), max(data(:,r)) , 20)};
				line(data(:,c), data(:,r) , 'color', markercolor, 'marker','.','linestyle','none','markersize', 5)
				% imagesc(edges{1}, edges{2}, hist3([data(:,c) data(:,r) ], 'Edges', edges ))
				axis xy
			% elseif c > r & ~plotdistributions
			% 	delete(hand(r,c))
			end


			if r == N
				title(VarNames{c})
				% ylabel(VarNames{c})
			end
			if c == N
				ylabel(VarNames{r})
			end
			if r ~= 1
				set(gca, 'xtick' , [])
			end
			if c ~= N & c ~= 1
				set(gca, 'ytick' , [])
			end
			if c == r
				set(gca,'ytick', [0 max(hhh)])
			end
			if c == 1
				set(gca,'yaxislocation', 'right')
			end


			try
			xlim([min(data(:,c)) max(data(:,c))])
			catch
			end


			set(gca,'color',bgcolor,'xcolor', axiscolor, 'ycolor', axiscolor,'tickdir', 'out','box','on');

		end

		c = 0;

	end

	colormap(jet(5))


	pl_i = [1:length(l1)];
	for r = pl_i

		linkaxes(hand(r,pl_i(find(pl_i~=r))),'y')
		% linkaxes(hand(pl_i(find(pl_i~=r)),r),'x')
		linkaxes(hand(pl_i,r),'x')

	end

	for c = pl_i

		linkaxes(hand(pl_i(find(pl_i~=c)),c),'x')
		linkaxes(hand(c, pl_i(find(pl_i~=c))),'y')

	end




	if nargout>0
		varargout{1} = hand;
	end


