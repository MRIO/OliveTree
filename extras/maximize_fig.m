% maximize_fig.m
function maximize_fig(varargin)

if isempty(varargin)

	MS = getMonitorSize;
	set(gcf, 'position', MS);
else
	MS = getMonitorSize;
	set(varargin{1}, 'position', MS);
end	