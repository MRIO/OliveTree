
% openfig_savepng.m

function saveallfigs(varargin)


	p = inputParser;
	
	p.addParameter('prefix', []) 
	p.addParameter('path', []) 
	p.addParameter('format', 'pdf')
	p.addParameter('style', '12x12')
	p.addParameter('savefig', 1)

	
	p.parse(varargin{:});
	prefix = p.Results.prefix;
	path = p.Results.path;
	format = p.Results.format;
	style = p.Results.style;
	savefig = p.Results.savefig;

	h = get(0,'children');

for i=1:length(h)
	h(i).Color = [1 1 1];
	if savefig
	   saveas(h(i), [prefix '_' num2str(i)], 'fig');
	end

	fname = [prefix '_' num2str(i) '.' format];
	snam=style;
	s=hgexport('readstyle',snam);
    s.Format = format;
    hgexport(h(i),fname,s);

end


