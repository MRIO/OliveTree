function colorlines_gradient(axh,color, varargin)

nlines = length(get(axh));

switch color
	case 'b'
		clrs = ones(nlines,3);
		clrs(:,1) = linspace(0,.9,nlines);
		clrs(:,2) = linspace(0,.9,nlines);
	case 'r'
		clrs = ones(nlines,3);
		clrs(:,2) = linspace(0,.9,nlines);
		clrs(:,3) = linspace(0,.9,nlines);
	case 'g'
		clrs = ones(nlines,3);
		clrs(:,1) = linspace(0,.9,nlines);
		clrs(:,3) = linspace(0,.9,nlines);
	case 'y'
		clrs = ones(nlines,3);
		clrs(:,3) = linspace(0,.9,nlines);
	case 'p'
		clrs = ones(nlines,3);
		clrs(:,2) = linspace(0,.9,nlines);
	case 'x'
		clrs = ones(nlines,3);
		clrs(:,1) = linspace(0,.9,nlines);
end

if nargin ==3
	clrs = varargin{1};
end
	


arrayfun(@(n) set(axh(n),'color', clrs(n,:)), [1:nlines]);

