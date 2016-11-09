
function hand = NDscatter(intable, groupby)

	
	
	
	bgcolor = [1 1 1];
	axiscolor = [.1 .1 .1];
	markercolor = [0 0 0];


	VarNames = intable.Properties.VariableNames;
	data = table2array(intable);
	N = length(VarNames);

	if isempty(groupby)
		groups = ones(size(data,2));
	else
		[C_ r_ groups] = unique(data(:,groupby));
	end
	
	ngroups = max(groups);
	groupcolors = linspecer(ngroups);




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
				[hhh xxx] = hist(data(c,:),30);
				bar(xxx, hhh, 'facecolor', markercolor)
				hold on
				for g = 1:ngroups
					[hhh] = hist(data(groups==g,c),xxx);
					plot(xxx, hhh,'color', groupcolors(g,:) )
				end

			end

			if c > r
				edges = {linspace(min(data(:,c)), max(data(:,c)) , 20); linspace(min(data(:,r)), max(data(:,r)) , 20)};
				% line(data(:,c), data(:,r) , 'color', markercolor, 'marker','.','linestyle','none','markersize', 5)
				imagesc(edges{1}, edges{2}, hist3([data(:,c) data(:,r) ], 'Edges', edges ))
				axis xy
			end


			if r == 1
				xlabel(VarNames{c})
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
