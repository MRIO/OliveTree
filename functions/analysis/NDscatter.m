
function varargout = NDscatter(varargin)



p = inputParser;
p.addRequired('intable')
p.addOptional('groupby',[]) 
p.addParamValue('fields',[])
p.addParameter('colors',[])
p.addParameter('hist2d',0)
p.addParameter('stackhist',0)

p.parse(varargin{:}) 

intable = p.Results.intable;
groupby = p.Results.groupby;
selectedfields  = p.Results.fields;
plotdistributions = 0;
colors = p.Results.colors;
hist2d = p.Results.hist2d;
stackhist = p.Results.stackhist;

VarNames = intable.Properties.VariableNames;
data = table2array(intable);
N = length(VarNames);

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
	
	ngroups = max(groups)

	try
		groupcolors = linspecer(ngroups);
	catch
		groupcolors  = [0 0 0];
	end

	if ngroups ==2
		% groupcolors = [.9 0 .2 ; .2 .7 .2];
		% groupcolors = [.9 0 .2 ; .2 .7 .2];
		groupcolors = [.2 .7 .2 ; .2 .2 .7];
	end

	if not(isempty(colors))
		groupcolors = colors;
	end

	figure('color', [1 1 1], 'colormap', groupcolors)


	l1 = linspace(0, .88, N+1); l1 = l1(1:end -1)+.08; 
	l2 = linspace(0, .88, N+1); l2 = l2(1:end -1)+.08; l2 = fliplr(l2);
	w  = .80/N;
	h  = .80/N;

	r = 0; c = 0;
	for rr = l1
		r = r+1;
		for cc = l2
			c = c+1;
			
			hand(r, c) = axes('position', [cc rr w h],'colororder', groupcolors);
			if r > c
				hold on
				for g = 1:ngroups

					line(data(groups==g,c), data(groups==g,r), ...
						 'markersize', 5, 'color', groupcolors(g,:), 'marker','.','linestyle','none','markersize',15)

				end
				% line(data(:,c), data(:,r) , 'color', markercolor, 'marker','.','linestyle','none','markersize', 5)
				hold off
			elseif r == c
				if min(data(:,r)) == max(data(:,r))
					continue
				end
				xxx = linspace(min(data(:,r)), max(data(:,r)), 30);

				if ngroups==2
					for g = 1:ngroups
						[hhh(:,g)] = hist(data(groups==g,r),xxx);
						% plot(xxx, hhh,'color', groupcolors(g,:) )
						hold on
					end

					if stackhist
						h_area = area(xxx, hhh);
					else
						for g = 1:ngroups
							h_area = patch([xxx(1) xxx xxx(end) xxx(1)] , [0; hhh(:,g) ; 0; 0], groupcolors(g,:))
						end

					end
					alpha(.7)
					colormap(groupcolors)
					
				else
					
					[hhh xxx] = hist(data(:,r), xxx);
					bar(xxx, hhh, 'facecolor', markercolor)
				end

				set(gca,'ytick', [0 max(hhh(:))])

			end

			if c > r % & plotdistributions
				
				if ~hist2d
					if ngroups ==2
						for g = 1:ngroups
							% line(data(groups==g,c), data(groups==g,r), ...
								 % 'markersize', 5, 'color', groupcolors(g,:), 'marker','.','linestyle','none','markersize',15)

							[groupcorr groupval] = corr(data(groups==g,c), data(groups==g,r));
							groupregress = regress(data(groups==g,c), data(groups==g,r));
							text_regression{g,1} = ['r=' num2str(groupcorr) '; p=' num2str(groupval) ];
							text_regression{g,2} = ['R^2=' num2str(groupregress)];

						end
					else
						
							[groupcorr groupval] = corr(data(:,c), data(:,r));
							groupregress = regress(data(:,c), data(:,r));
							text_regression{1,1} = ['r=' num2str(groupcorr) '; p=' num2str(groupval) ];
							text_regression{1,2} = ['R^2=' num2str(groupregress)];
					end

					line(data(:,c), data(:,r) , 'color', [.8 .8 .8], 'marker','.','linestyle','none','markersize', 10)
					bla = axes('position', get(gca, 'position'));
					axis off	
					t = text(double(min(data(:,c))), double(min(data(:,r))), text_regression);
					t.VerticalAlignment = 'bottom';

					% axis off


					

				else
					nedges = 20;
					edges = {linspace(min(data(:,c)), max(data(:,c)) , nedges); linspace(min(data(:,r)), max(data(:,r)) , nedges)};
					HH = hist3([data(:,c) data(:,r) ]);
					imagesc(edges{1}, edges{2}, HH ); colormap(bone)
					axis xy

				end
			% elseif c > r & ~plotdistributions
			% 	delete(hand(r,c))
			end


			if r == N
				title(hand(r,c), VarNames{c},'interpreter', 'none')
				% ylabel(VarNames{c})
			end
			if c == N
				ylabel(hand(r,c), VarNames{r},'interpreter', 'none')
			end
			if r ~= c
				set(hand(r,c), 'xtick' , [])
			end
			if c ~= N & c ~= 1
				set(hand(r,c), 'ytick' , [])
			end
			if c == r
				set(hand(r,c),'ytick', [0 max(hhh(:))/2],'yaxislocation', 'left')
				% set(hand(r,c),'xtick', [linspace(min(data(:,c)), max(data(:,c)), 5)])

			end
			if c == 1 & (c~=r)
				set(hand(r,c),'yaxislocation', 'right')
			end


			try
			xlim([min(data(:,c)) max(data(:,c))])
			catch
			end


			set(hand(r,c),'color',bgcolor,'xcolor', axiscolor, 'ycolor', axiscolor,'tickdir', 'out','box','on');

		end

		c = 0;

	end

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


