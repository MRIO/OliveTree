function out = retrieve_data_M(varargin)

p = inputParser;
p.addRequired('sim')  % a matrix with two columns or a cell array with two cells;
p.addParamValue('FOI', []) % stdandard deviation criterion for offset threshold

p.parse(varargin{:});

sim = p.Results.sim;
FOI = p.Results.FOI;


networkHistory  = sim.networkHistory;
rows 			= sim.rows;
columns			= sim.columns;
lastState 		= sim.lastState;
duration 		= sim.duration;
dt 				= sim.dt;
g_CaL 			= sim.g_CaL;
g_Gap 			= sim.g_Gap;
time 			= sim.time;

if isempty(FOI)
	FOI= 'V_soma'
end

samples = duration*(1/dt);

traces = zeros(rows*columns,samples);

eval(['rd = @(i,k) networkHistory(i,k).' FOI  ';' ])

for i = 1:rows*columns

		for k = 1:samples

		traces(i,k) = rd(i,k);

		end

end


out = sim;

eval(['out.' FOI ' = traces(1:rows*columns,1:end)']);


V_soma = traces(1:rows*columns,1:end);
% make this a function that operates on the structure fields
V_soma_traces = reshape(V_soma,rows*columns,[])';


out.traces.V_soma = V_soma;
out.traces.V_soma_unwrapped = V_soma_traces;
% out.traces.V_soma_unwrapper = @reshape(V_soma,rows*columns,[])';





