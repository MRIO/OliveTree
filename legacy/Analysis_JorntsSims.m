

set(0, 'defaultfigurecolor', ones(3,1))
set(0, 'defaultfigurecolormap', linspecer(64))
set(0, 'defaultfaxescolororder', linspecer(100))


%
	% time = 200; %ms 
	% I_app = zeros(rows,columns,fs*time);

	% % perturbation = zeros(rows, columns); perturbation(ceil(rows/2),ceil(columns/2))=1;
	% perturbation = rand(rows,columns);
	% stim_on = [20 50 100]*fs;
	% % I_app(:,:,stim_on) = rand(rows,columns,length(stim_on))*-100;

	% I_app(:,:,stim_on) = repmat(perturbation,[1,1,length(stim_on)])*1000;
	% [sim] = IOtoplevelCa_M('rows',rows,'columns', columns,'appCurrent',I_app,'time',time,'delta',dt,'g_CaL',g_CaL);
	% sim_r = retrieve_data(sim);
	% replayResults(sim_r)


	% time = 100; %ms 
	% I_app = zeros(rows,columns,1/dt*time);
	% stim_on = [0:200*fs:time*fs];
	% % I_app(:,:,stim_on) = rand(rows,columns,length(stim_on))*-100;
	% I_app(:,:,stim_on) = repmat(rand(rows,columns),[1,1,length(stim_on)])*100;
	% [sim] = IOtoplevelCa_M('rows',rows,'columns', columns,'appCurrent',I_app,'time',time,'delta',dt);


	% [hist,steps] = IOtoplevelCa(dt,time,rows,columns,gap_c,1,zeros(rows,columns,time/dt),nan(rows,columns,time/dt), g_CaL);

	% imagesc(V_soma_traces);



% load sim_results_with_stim_1
% load sim_results_with_stim_2
% load sim_results_with_stim_3
% load sim_results_with_stim_4
% load sim_results_with_stim_5
% load sim_results_with_stim_6
% load sim_results_with_stim_7
% load sim_results_with_stim_8

path_to_sims = ['/Users/M/Cogsci/Experiments/Olive/JorntsOlive/Simulations for Mario/Synchrony sims/'];
folder = ['1/'];
fname = 'simresults1';
fs  = 20; %KHz
data_path = [path_to_sims folder fname];
load(data_path)
duration = size(traces,3)/fs;
traces_ms = traces(1:10,1:10,1:fs:end);
unwrapped_traces_1 = reshape(traces_ms,100,duration)';

path_to_sims = ['/Users/M/Cogsci/Experiments/Olive/JorntsOlive/Simulations for Mario/Synchrony sims/'];
folder = ['1/'];
fname = 'simresults3';
data_path = [path_to_sims folder fname];
load(data_path)
duration = size(traces,3)/fs;
traces_ms = traces(1:10,1:10,1:fs:end);
unwrapped_traces_2 = reshape(traces_ms,100,duration)';

duration = size(traces,3);


subplot(131); imagesc(unwrapped_traces_1)
subplot(132); imagesc(unwrapped_traces_2)
subplot(133); imagesc(unwrapped_traces_1-unwrapped_traces_2)
xlabel('neurons')
ylabel('simulated t (ms)')


%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%


results_folder = 1;
folder = [];
fnumber = 2;
path_to_sims = ['/Users/M/Cogsci/Experiments/Olive/JorntsOlive/Simulations for Mario/sim_results_' num2str(results_folder) '/'];
fname = ['simresults_with_stim_'  num2str(fnumber) '.mat'];
data_path = [path_to_sims folder fname];
load(data_path)


fs  = 20;
duration = size(traces,3)/fs;

traces_ms = traces(1:10,1:10,1:fs:end);
unwrapped_traces = reshape(traces_ms,100,duration)';

appCurrent_ms = appCurrent(1:10, 1:10,1:fs:end);
stimulus_ms = reshape(appCurrent_ms,100,duration)';
[sumstim_t stim_n] = find(sum(stimulus_ms,2));

stim_t = sumstim_t(find(diff([0 ; stim_t])>1));

subplot(2,3,1:3]), plot(unwrapped_traces)
subplot(2,3,4), plot(unwrapped_traces)        
subplot(2,3,5), plot(unwrapped_traces)
subplot(2,3,6), plot(unwrapped_traces)
ylim([-70 20])
ylim([-70 30])                        
ylim([-70 30])
ylim([-70 30])
xlabel('ms')
ylabel('mV')


