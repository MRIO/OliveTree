function out = plot_mean_and_std(varargin)

p = inputParser;
	p.addOptional('x',[])
	p.addRequired('data')
	p.addParamValue('color', [1 0 0])
	p.addParamValue('quantiles', [.1 .9])

p.parse(varargin{:});

	x 			= p.Results.x;
	data 		= p.Results.data;
	color 		= p.Results.color;
	quantiles	= p.Results.quantiles;


if isempty(x)
	x = [1:size(data,2)];
end

	sigm = @(x) 1./(1+exp(-x));


	M = mean(data);
	S = std(data);
	Q = quantile(data, quantiles);


	fill_between_lines = @(X,Y1,Y2, color) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], color ,'edgecolor','none');


	fill_between_lines(x,Q(1,:),Q(2,:), sigm([.8 .8 .8]+color)) 



	line(x, M,'color', color, 'linewidth',3)


out.mean 		= M; 
out.std  		= S;
out.quantiles 	= Q;