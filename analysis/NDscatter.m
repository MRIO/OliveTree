
function varargout = NDscatter(varargin)



p = inputParser;
p.addRequired('intable')
p.addOptional('groupby',[]) 
p.addParamValue('fields',[])

p.parse(varargin{:}) 

intable = p.Results.intable;
groupby = p.Results.groupby;
selectedfields  = p.Results.fields;

VarNames = intable.Properties.VariableNames;
data = table2array(intable);
N = length(VarNames);



if ~isempty(groupby)
	[C_ r_ groups] = unique(data(:,groupby));
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



	l1 = linspace(0, .9, N+1); l1 = l1(1:end -1)+.05; 
	l2 = linspace(0, .9, N+1); l2 = l2(1:end -1)+.05; l2 = fliplr(l2);
	w  = .9/N;
	h  = .9/N;

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
				xxx = linspace(min(data(:,r)), max(data(:,r)), 30);
				[hhh xxx] = hist(data(:,r), xxx);
				bar(xxx, hhh, 'facecolor', markercolor)
				hold on
				% for g = 1:ngroups
				% 	[hhh] = hist(data(groups==g,r),xxx);
				% 	plot(xxx, hhh,'color', groupcolors(g,:) )

				% end
				% axis tight

			end

			if c > r
				edges = {linspace(min(data(:,c)), max(data(:,c)) , 20); linspace(min(data(:,r)), max(data(:,r)) , 20)};
				line(data(:,c), data(:,r) , 'color', markercolor, 'marker','.','linestyle','none','markersize', 5)
				% imagesc(edges{1}, edges{2}, hist3([data(:,c) data(:,r) ], 'Edges', edges ))
				axis xy
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
			if c ~= N
				set(gca, 'ytick' , [])
			end
		
			set(gca,'color',bgcolor,'xcolor', axiscolor, 'ycolor', axiscolor);
		end

		c = 0;
	end

	colormap(jet(5))


	pl_i = [1:length(l1)];
	for r = pl_i

		linkaxes(hand(r,pl_i(find(pl_i~=r))),'y')
		linkaxes(hand(pl_i(find(pl_i~=r)),r),'x')

	end


	if nargout>0
		varargout{1} = hand;
	end


