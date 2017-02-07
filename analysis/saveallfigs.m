
% openfig_savepng.m

function saveallfigs(varargin)


	p = inputParser;
	
	p.addParameter('prefix', []) 
	p.addParameter('path', []) 
	p.addParameter('formats', 'style')
	p.addParameter('style', '9x6')

	
	p.parse(varargin{:});
	prefix = p.Results.prefix;
	path = p.Results.path;
	formats = p.Results.formats;
	style = p.Results.style;

	mag = '-m4'

	h = get(0,'children');

for i=1:length(h)
	h(i).Color = [1 1 1];
   saveas(h(i), [prefix num2str(i)], 'fig');


   switch formats
   case 'style'
   		fname = [prefix num2str(i) '.png'];
   		snam=style;
   		s=hgexport('readstyle',snam);
	    s.Format = 'png';
	    hgexport(h(i),fname,s);


   	case 'png'
   		export_fig([prefix num2str(i)], '-png',mag,h(i))
   	case 'svg'
   		print(h(i), '-dsvg', [prefix num2str(i)] )
	case 'pdf'
   		print(h(i), '-dpdf', [prefix num2str(i)] )
   	otherwise

   	end
   % plot2svg([num2str(i) '.svg'],h(i));
end


