function out = spikedetect(varargin) 
  % function out = spikedetect(sim, plotme, printme) 
  % Detect spikes with threshold
  % 
  % written by Elisabetta Iavaronne

  p = inputParser;
  p.addRequired('sim')
  p.addParamValue('plotme', 0)
  p.addParamValue('printme', 0)
  p.addParamValue('thresh', -30)
  p.addParamValue('returnbinarymatrix', 1)
  p.addParamValue('time_slice', [])

  p.parse(varargin{:});

  sim = p.Results.sim;
  plotme = p.Results.plotme;
  printme = p.Results.printme;
  spkthresh = p.Results.thresh;
  returnbinarymatrix = p.Results.returnbinarymatrix;
  tslice = p.Results.time_slice;

try
  netsize = sim.networkParameters.networksize;
catch
  netsize = sim.networksize;
end

noneurons = prod(netsize);

if isempty(tslice)
  tslice = [1:size(sim.networkHistory.V_soma,2)];
end

if isfield(sim, 'fname')
    fname = sim.fname;
else
    fname = '';
end

if printme && ~plotme
    plotme = 1;
end

if printme
    visible = 'off';
else
    visible = 'on';
end

A = zeros(noneurons,21);


VS = sim.networkHistory.V_soma(:,tslice);

simtime = size(VS,2);

for i = 1:noneurons
   abovethreshold{i} = find(VS(i,:)> spkthresh);
   b{i} = diff([-2 abovethreshold{i}]);
   spikeonset = find(b{i}>1);
   d{i} = abovethreshold{i}(spikeonset) + tslice(1);
   ISI{i} = diff(d{i});
   if ISI{i}>0

    freq{i} = 1./ISI{i}*1000; % in Hz
  else
    freq{i} = 0;
  end
  
   med{i} = median(freq{i});
    
end




if returnbinarymatrix
  binaryspikes = zeros(noneurons, simtime);
  for i = 1:noneurons
    binaryspikes(i,d{i}) = 1;
  end
end




% Find proportion of spiking cells
spkneurons = [find(~cellfun(@isempty,d))];
propspkneurons = numel(spkneurons)/noneurons;

 
for i = 1:noneurons;
  spkspercell(i) = numel(d{i});
end

popmedian = cell2mat(med);
popfreq = sum(spkspercell)/noneurons/simtime*1e3;

if plotme
  fighandle = figure('visible', visible);
  out.fhandle = fighandle;

  edges1 = [-1 0:.25:max(popmedian)];

  histmed = histc(popmedian,edges1);

  bar(edges1,histmed,'edgecolor', 'none', 'facecolor', [0 0 0]);
  xlabel('Median frequency (Hz)'); ylabel(('Num. neurons'));
  xlim([-1 10]);
  text(5,max(histmed)*0.9, ['Spk. neurons:' num2str(propspkneurons)]);

  figure, boxplot(popmedian);
  text(1.2, 1, ['Spk. neurons:' num2str(propspkneurons)]);
  ylabel('ISF (Hz)');

  figure, bar3(A);
  ylim([0 noneurons+1]);
  set(gca,'xTick',[0:1:21]);
  set(gca,'xTickLabel',[-1 0:0.5:10]);
  xlabel('Hz'); ylabel('neurons'); zlabel('spike count');


    if printme
        plot2svg([fname '_spiking_hist.svg'], fighandle)
        close(fighandle)
    end
end
    


out.spikespercell = spkspercell;
out.spikes = d;
out.medfreq = popmedian;
out.popfrequency = popfreq;
out.propspkneurons = propspkneurons;
out.cellISI = ISI;

if returnbinarymatrix
  out.binaryspikes = binaryspikes;
end




