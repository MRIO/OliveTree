% joinsim.m
function out = joinsim(sims, simstojoin)


spikespercell = 0;
binspikes = [];


numn = prod(sims{1}.networksize);
spikeset = cell(1, numn);


d = sims{simstojoin(1)}.duration;

if isfield(sims{1}.perturbation, 'triggers');
	p = length(sims{1}.perturbation.triggers);
	triggers = cell(p);
end


VSomaHist  = [];
i = 1;
for s = simstojoin

	thesespikes = sims{s}.spikes;

	binspikes = [binspikes thesespikes.binaryspikes];
	
	spikespercell = spikespercell + thesespikes.spikespercell;

	for n = 1:numn
		spikeset{n} = union(spikeset{n}, thesespikes.spikes{n}+d*(i-1));
	end

	if isfield(sims{1}.perturbation, 'triggers');
		for t = 1:p
			triggers{t} = [triggers{t} ; ((i - 1)*d + sims{s}.perturbation.triggers{t})];
		end
	end

	d = sims{s}.duration;

	if isfield(sims{s}.networkHistory, 'V_soma');
		VSomaHist = [VSomaHist sims{s}.networkHistory.V_soma];
	end

	i = i +1;




end

if isfield(sims{simstojoin(1)}.perturbation, 'triggers');
	% pruning only the first pulse of the triggers
	for t = 1:p
		triggers{t} = triggers{t}((diff([1; triggers{t}]))>4);
	end
end


out = sims{simstojoin(1)};
out.duration = d*length(simstojoin);
out.spikes.spikespercell = spikespercell;
out.spikes.spikes  = spikeset;
out.spikes.binaryspikes = binspikes;

if isfield(sims{simstojoin(1)}.perturbation, 'triggers');
	out.perturbation.triggers = triggers;
end

if isfield(sims{simstojoin(1)}.networkHistory, 'V_soma');
	out.networkHistory.V_soma = VSomaHist;
end

