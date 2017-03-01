function cells =  jitter_cell_parameters(cells, varargin)

fieldstojitter = {'g_CaL'};
noneurons = length(cells.g_h);


if nargin ==2
	spread = varargin{1}
else
	spread = .1;
end


for ff = fields(cells)'
	if ismember(fieldstojitter, ff)
			
			maxP = max(eval(['cells.' ff{1}]));
			minP = min(eval(['cells.' ff{1}]));
			param_interv = maxP - minP;
			if param_interv ~=0
				str = ['cells.' ff{1} ' = rand(noneurons,1) * param_interv + minP;' ];
				eval(str)
			else
				str = ['cells.' ff{1} ' = cells.' ff{1} '+ spread * randn(noneurons,1) .* mean(cells.' ff{1} ');' ];
				eval(str)
			end
	end
end

