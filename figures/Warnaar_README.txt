***********************
Warnaar_README.txt
***********************

To run the model and simulations the whole tree of the model must be in the path

addpath(genpath('OliveTree'))

To reproduce the analysis and panels from Warnaar et al. (2018), one can either use the precomputed data provided in dataWarnaar\ or recompute the simulations. Remember to add the data folder to the path.


NOTE: Running the model requires a CUDA enabled GPU.

All analysis is performed in these scripts:

- Warnaar_all_simulations_analyses_and_figures.m
	Produces thorough network analysis including PSTH's and analysis of phase resets.

- Warnaar_all_prc_figures
	Produces phase response analysis

- Warnaar_model_figures
	Produces the panels for the model introduction figure


- analysis_NetPspace.m
	Returns the analysis of the parameter space for tau x eta.
	Must have the results of the snippet below in the workspace.

% [=================================================================]
%  tau x eta NetPspace
% [=================================================================]

% Lastrun: 20/6/2016
% if strcmp(pwd,'/home/titanuser1/Titan/Bench2')
% 	X_README = 'THIS IS A SIMULATION FOR TAU - gapfactor introduced'
% 	p1  = [5 10 20 40 60 80];		% time constant of ornstein uhlenbeck process (tau)
% 	p2  = [0:0.05:.3];	% sametoall (noise correlation) of ornstein uhlenbeck process
% 	p3  = [-.6];		% noiseamp of ornstein uhlenbeck process
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

% 	NetPspace;
% end




