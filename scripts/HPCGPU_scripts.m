% HPCGPU_scripts.m

% NOTES:
% nameprefix must not contain spaces

%============================= 4 Jochen ==============================%
% X_README = 'XCORR sims 7/6/2016';
% gaps = [0.04] ; gapcomp = 0   ; tau = 20; noisesig = .7; noisemu = -1; sametoall = 0.3; simtype = '1Hz'   ; conntype = 'cluster' ; HPCGPU_periodic_ampa;
% gaps = [0.04] ; gapcomp = 0   ; tau = 20; noisesig = .7; noisemu = -1; sametoall = 0.3; simtype = 'spont' ; conntype = 'cluster' ; HPCGPU_periodic_ampa;
% gaps = eps; gapcomp = -1.6	; tau = 20; noisesig = .7; noisemu = -1; sametoall = 0.3; simtype = '1Hz'   ; conntype = 'cluster' ; HPCGPU_periodic_ampa;	
% gaps = eps; gapcomp = -1.6	; tau = 20;	noisesig = .7; noisemu = -1; sametoall = 0.3; simtype = 'spont' ; conntype = 'cluster' ; HPCGPU_periodic_ampa;	

%============================= 4 Pascal ==============================%

% if strcmp(pwd,'/home/titanuser1/Titan/Bench')
% 	X_README = 'sims for pascal_'
% 	nameprefix = 'moreoscillations_nocorr_'
% 	gaps = 0.04; tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = '1Hz'   ;  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 3: 
% 	gaps = 0.04; tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = 'spont' ;  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 1: 
% 	gaps = 0.04; tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = 'gallop';  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 1: 

% 	gaps = 0;	 tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = '1Hz'   ;  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 3: 
% 	gaps = 0;	 tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = 'spont' ;  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 1: 
% 	gaps = 0;	 tau = 30; noisesig = .4; noisemu = -.4; sametoall = 0; simtype = 'gallop';  conntype = 'iso'     ; HPCGPU_periodic_ampa		;% 4Pascal 1: 


% end

%=============================recomputing to generate all vsoma output, also, add random seed to gauge variability ==============================%
% if strcmp(pwd,'/home/titanuser1/Sync/Titan/Bench')
% 	X_README = 'replay_06_12_16 for pascal_'
% 	nameprefix = 'replay_06_12_16_'


% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 4: 
% 	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 2: 

% 	seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.2; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
% 	seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.2; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa		;% 4Pascal 4: 
% 	seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.2; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa		;% 4Pascal 2: 

% 	seed = 1; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
% 	seed = 1; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 4: 
% 	seed = 1; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.2; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 2: 

	

% end



if strcmp(pwd,'/home/titanuser1/Sync/Titan/Bench4')
	X_README = 'replay_06_12_16 for pascal_with_spont'
	nameprefix = 'replay_06_12_16_with_spont'

	seed = 0; tau = 20; noisesig = 0;  noisemu = 0;	  sametoall = 0.1; simtype = 'spont' ; conntype = 'iso' ; numruns = 1;  HPCGPU_periodic_ampa		;% 4Pascal 2: 

	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.1; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.1; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 4: 
	seed = 0; tau = 20; noisesig = .6; noisemu = -.6; sametoall = 0.1; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  HPCGPU_periodic_ampa		;% 4Pascal 2: 

	% seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.1; simtype = 'gallop'; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa	 	;% 4Pascal 2: 
	% seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.1; simtype = '1Hz'   ; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa		;% 4Pascal 4: 
	% seed = 0; tau = 20; noisesig = .3; noisemu = -.6; sametoall = 0.1; simtype = 'spont' ; conntype = 'iso' ; numruns = 4;  gaps = eps; HPCGPU_periodic_ampa		;% 4Pascal 2: 


end




% [=================================================================]
% clustersize
% [=================================================================]

% X_README = 'CLUSTER SIZE SIMULATION'
% nameprefix = 'clustersize';
% 	p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0 .1 .2 .3];			% sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [0];	% noiseamp of ornstein uhlenbeck process
% 	p4  = [0.04 0.001]; % average gap leak per connection
% 	p5  = [3];			% single cell connection radius 
% 	p6  = [10 25 50 100 200];			% clustersize (iff conntype='cluster')
% 	p7  = [1];			% intraclusterP (iff conntype='cluster')
% 	p8  = [.05];		% extraclusterP (iff conntype='cluster')
% 	p9  = [8];			% mean number of connections
% 	p10 = [10];			% N, size of N x N network
% 	p11 = [2];			% depth in Z		
% 	p12 = [1];			% whether using gap compensation in cells
% 	p13 = [.7];			% variability of noise injected

% 	 simtime = 5000;
% 	 conntype = 'cluster';

% 	 NetPspace

% [=================================================================]
%  sig x tau
% [=================================================================]
% Lastrun: 7/6/2016
% X_README = 'THIS IS A SIMULATION FOR TAU - gapfactor introduced'
% nameprefix = 'tauXncorr';
% 	p1  = [5 10 20 30 40 50];	% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0 .3];	    % sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [0];			% noiseamp of ornstein uhlenbeck process
% 	p4  = [0.04 eps]; 	% average gap leak per connection
% 	p5  = [3];			% single cell connection radius 
% 	p6  = [36];			% clustersize (iff conntype='cluster')
% 	p7  = [1];			% intraclusterP (iff conntype='cluster')
% 	p8  = [0];			% extraclusterP (iff conntype='cluster')
% 	p9  = [8];			% mean number of connections
% 	p10 = [6];			% N, size of N x N network
% 	p11 = [2];			% depth in Z		
% 	p12 = [1];	% whether using gap compensation in cells
% 	p13 = [.5];	    % variability of noise injected (sig) - 1Hz when noise_amplitude = 0, gap = 0.04 and n = 8;

% Lastrun: 20/6/2016
% if strcmp(pwd,'/home/titanuser1/Titan/Bench2')
% 	X_README = 'THIS IS A SIMULATION FOR TAU - gapfactor introduced'
% 	p1  = [5 10 20 40 60 80];		% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0:0.05:.3];	    % sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [-.6];			% noiseamp of ornstein uhlenbeck process
% 	p4  = [0.04 eps]; 	% average gap leak per connection
% 	p5  = [2];			% single cell connection radius 
% 	p6  = [36];			% clustersize (iff conntype='cluster')
% 	p7  = [1];			% intraclusterP (iff conntype='cluster')
% 	p8  = [0];			% extraclusterP (iff conntype='cluster')
% 	p9  = [8];			% mean number of connections
% 	p10 = [6];			% N, size of N x N network
% 	p11 = [2];			% depth in Z		
% 	p12 = [1];			% whether using gap compensation in cells
% 	p13 = [.6];	    % variability of noise injected (sig) - 1Hz when noise_amplitude = 0, gap = 0.04 and n = 8;

% 	conntype = 'iso';

% 	simtime = 3000;
% 	conntype = 'iso';

% 	NetPspace;
% end

%============================= gap compensation tests ==============================%

 % conntype = 'iso' 
 % X_README = 'THIS IS A SIMULATION FOR GAP COMPENSATION'
 % nameprefix = 'gapcomp';
	% p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
	% p2  = [0 .1 .2];	% sametoall (noise correlation) of ornstein uhlenbeck process
	% p3  = [0];			% noiseamp of ornstein uhlenbeck process
	% p4  = [0.04 eps]; 	% average gap leak per connection
	% p5  = [3];			% single cell connection radius 
	% p6  = [36];			% clustersize (iff conntype='cluster')
	% p7  = [1];			% intraclusterP (iff conntype='cluster')
	% p8  = [0];			% extraclusterP (iff conntype='cluster')
	% p9  = [8];			% mean number of connections
	% p10 = [6];			% N, size of N x N network
	% p11 = [2];			% depth in Z		
	% p12 = [-1:.5:2];			% whether using gap compensation in cells
	% p13 = [.5];		% variability of noise injected

	% simtime = 5000;
	% conntype = 'iso';

	% % NetPspace;


% [=================================================================]
%  noise tests
% [=================================================================]

	% p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
	% p2  = [0];			% sametoall (noise correlation) of ornstein uhlenbeck process
	% p3  = [-1:.25:1];	% noiseamp of ornstein uhlenbeck process
	% p4  = [0.04 0.001]; % average gap leak per connection
	% p5  = [3];			% single cell connection radius 
	% p6  = [20];			% clustersize (iff conntype='cluster')
	% p7  = [.9];			% intraclusterP (iff conntype='cluster')
	% p8  = [.1];			% extraclusterP (iff conntype='cluster')
	% p9  = [8];			% mean number of connections
	% p10 = [6];			% N, size of N x N network
	% p11 = [2];			% depth in Z		
	% p12 = [0];			% whether using gap compensation in cells
	% p13 = [0:.2:1];		% variability of noise injected

	% NetPspace;


% [================================================]
%  sig x mu (corr =0)
% [================================================]


% if strcmp(pwd,'/home/titanuser1/Titan/Bench2')
% 	X_README = 'SIG X MU Corr = 0'
% 	nameprefix = 'sigXmu_0corr';
% 	p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0];			% sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [-.6 :.1: 0];	% noiseamp of ornstein uhlenbeck process
% 	p4  = [0.04 eps]; 	% average gap leak per connection
% 	p5  = [3];		% single cell connection radius 
% 	p6  = [36];			% clustersize (iff conntype='cluster')
% 	p7  = [1];			% intraclusterP (iff conntype='cluster')
% 	p8  = [0];			% extraclusterP (iff conntype='cluster')
% 	p9  = [6];			% mean number of connections
% 	p10 = [15];			% N, size of N x N network
% 	p11 = [2];			% depth in Z		
% 	p12 = [1];			% whether using gap compensation in cells
% 	p13 = [0:.1:.6];	    	% variability of noise injected (sig) - 1Hz when noise_amplitude = 0, gap = 0.04 and n = 8;

% 	simtime = 3000;
% 	conntype = 'iso';

% 	NetPspace;
% 	clear;
% end


% [================================================]
%  sig x mu (corr =0.1)
% [================================================]


if strcmp(pwd,'/home/titanuser1/Sync/Titan/Bench2')
	X_README = 'SIG X MU Corr = 0.1 R2'
	nameprefix = 'sigXmu_01corr_R2';
	p1  = [20];			% time constant of ornstein uhlenbeck process (tau)
	p2  = [0.1];			% sametoall (noise correlation) of ornstein uhlenbeck process
	p3  = [-.6 :.1: 0];	% noiseamp of ornstein uhlenbeck process
	p4  = [0.04 eps]; 	% average gap leak per connection
	p5  = [2];		% single cell connection radius 
	p6  = [36];			% clustersize (iff conntype='cluster')
	p7  = [1];			% intraclusterP (iff conntype='cluster')
	p8  = [0];			% extraclusterP (iff conntype='cluster')
	p9  = [6];			% mean number of connections
	p10 = [15];			% N, size of N x N network
	p11 = [2];			% depth in Z		
	p12 = [1];			% whether using gap compensation in cells
	p13 = [0:.1:.6];	    	% variability of noise injected (sig) - 1Hz when noise_amplitude = 0, gap = 0.04 and n = 8;

	simtime = 3000;
	conntype = 'iso';

	NetPspace;
	clear;
end




% [================================================]
%  connectivity radius
% [================================================]



% if strcmp(pwd,'/home/titanuser1/Sync/Titan/Bench2')
% 	X_README = 'THIS IS A SIMULATION FOR RADIUS X NCORR'
% 	nameprefix = 'RADIUS_X_NCORR';
% 	p1  = [20];		% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0.1];	% sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [-.6];		% noiseamp of ornstein uhlenbeck process
% 	p4  = [0.04 eps]; 	% average gap leak per connection
% 	p5  = [1:5];		% single cell connection radius 
% 	p6  = [36];			% clustersize (iff conntype='cluster')
% 	p7  = [1];			% intraclusterP (iff conntype='cluster')
% 	p8  = [0];			% extraclusterP (iff conntype='cluster')
% 	p9  = [4:2:12];			% mean number of connections
% 	p10 = [15];			% N, size of N x N network
% 	p11 = [3];			% depth in Z		
% 	p12 = [1];			% whether using gap compensation in cells
% 	p13 = [.6];	    	% variability of noise injected (sig) - 1Hz when noise_amplitude = 0, gap = 0.04 and n = 8;

% 	conntype = 'iso';

% 	simtime = 3000;
% 	conntype = 'iso';

% 	NetPspace;
% end


