function out = sync_triggered_X(simulation)


% id = 3;
% simulation = simresults{id};

vsoma = simulation.networkHistory.V_soma;
noise = simulation.networkHistory.backgroundnoise;
% noise = resample(noise', 1,1/simulation.simulationParameters.dt)';
spks = spikedetect(simulation);
simtime = length(vsoma);
binsize = 30;
synctrignoise = 1;

win = 300;
noneurons = size(vsoma,1);

% calcula a probabilidade de N eventos simulationultaneos


% encontra aqueles com probabilidade menor que 1 %


% spike triggered noise average

perispikevsoma = [];
perispikenoise = [];

for c = 1:noneurons
	somasnip = @(spk) vsoma(c,spk-win: spk+win)';
	perispikevsoma = [perispikevsoma ; cell2mat(arrayfun(somasnip, spks.spikes{c}(3:end-3),'uniformoutput', 0))'];

	noissnip = @(spk) noise(c,spk-win: spk+win)';
	perispikenoise = [perispikenoise ; cell2mat(arrayfun(noissnip, spks.spikes{c}(3:end-3),'uniformoutput', 0))'];
end


allspks = cell2mat(spks.spikes); 
allspks = allspks(:);
spkcount = accumarray(allspks,1);
[count time] = hist(allspks, [1:binsize:simtime]);
spkdensity = ksdensity(allspks,time,'kernel', 'epanechnikov', 'bandwidth',binsize);


if synctrignoise

	% i  = 0; 
	for s = 1:max(count) % S = COUNT IS NUMBER OF SYNC spikes
			% i = i+1;
			t = time(find(count==s)); % find the bins where N sync spikes happened
			t(t<win+1) = [];
			t(t>(simtime-win)) = [];
			n = length(t); % how many times N sync happened

			perispikevsoma_perspikecount = [];
			perispikenoise_perspikecount = [];

			s

		if n >= 1

			for c = 1:noneurons

					% somasnip =  @(spk) vsoma(c,spk-win: spk+win)';
					% perispikevsoma_perspikecount = [perispikevsoma_perspikecount ; cell2mat(arrayfun(somasnip, t,'uniformoutput', 0))'];

					noissnip = @(spk) noise(c,spk-win: spk+win)';
					perispikenoise_perspikecount = [perispikenoise_perspikecount ; cell2mat(arrayfun(noissnip, t,'uniformoutput', 0))'];

					vsomasnip = @(spk) vsoma(c,spk-win: spk+win)';
					perispikenvsoma_perspikecount = [perispikevsoma_perspikecount ; cell2mat(arrayfun(vsomasnip, t,'uniformoutput', 0))'];
			end

			perispikenoise_syncspk(s,:) = mean(perispikenoise_perspikecount);

		end

		syncspks(s) = n;

	end

	% [================================================]
	%  sync triggered coherence
	% [================================================]
	% 
	% how much of spike synchrony is explained by subthreshold coherence?
	% how much is explained simply by noise amplitude?


	plotthis = 0;
	if plotthis

		figure
			subplot(2,1,1)
			hold on
			plot_mean_and_std(-win:win, perispikenoise)
			legend({'noise' })
			axis tight

			subplot(2,1,2)
			hold on
			plot_mean_and_std(-win:win, perispikevsoma)
			legend({'soma' })
			axis tight

	end

	out.perispikenoise_syncspk = perispikenoise_syncspk;

end


figure
bar(time, count)


%simultaneous events


out.perispikenoise = perispikenoise;
out.perispikevsoma = perispikevsoma;
out.num_syncspks = syncspks;

out.hist.count = count;
out.hist.time  = time;

