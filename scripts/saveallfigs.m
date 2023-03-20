
% openfig_savepng.m

function saveallfigs(varargin)


	p = inputParser;
	
	p.addParameter('prefix', []) 
	p.addParameter('path', []) 
	p.addParameter('format', 'pdf')
    p.addParameter('figsize', [6 6])
	p.addParameter('style', 'default')
	p.addParameter('savefig', 0)
	p.addParameter('exportfig', 1)

	
	p.parse(varargin{:});
	prefix = p.Results.prefix;
	path = p.Results.path;
	format = p.Results.format;
	style = p.Results.style;
	savefigs = p.Results.savefig;
	exportfig = p.Results.exportfig;
    figsize = p.Results.figsize;

	h = get(0,'children');

for i=1:length(h)
	h(i).Color = [1 1 1];
	fname = [prefix '_' num2str(i) ]
	if savefigs
		figure(h(i))
	   savefig(h(i), [fname '.fig']);

	end

	if exportfig
		fname = [prefix '_' num2str(i) '.' format];

            h(i).Units = 'centimeters';
            h(i).PaperUnits = 'centimeters';
            h(i).OuterPosition = [0 0 figsize(1) figsize(2)];
            h(i).Children.FontSize = 7;
	        exportgraphics(f100,[name '.pdf'],'Resolution',300)



		snam=style;
		s=hgexport('readstyle',snam);
	    s.Format = format;
	    hgexport(h(i),fname,s);
	end

end


