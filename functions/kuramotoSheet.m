function out = kuramotoSheet(varargin)
% function out = kuramotoSheet([N,M], K, 'Param', 'Value')]
% 
% Calculates the sheet of kuramoto oscillators with different connectivity schemes.
% User can select the oscillator frequencies, initial conditions, and 
% time base (dt for the forward euler solver)
% The script also outputs the synchronization (kuramoto parameter), i.e., the 
% centroid of the phase of the group of oscillators.
% 
% inputs:
% 	networksize, [N M]
% 	scaling for coupling - K
% outputs:
% 	out.state = sin ( theta_t) ;
% 	out.phase = phase;
% 	out.parameters = parameters
% 	out.oscillators = intrinsic frequency of oscillators;
% 	out.orderparameter = synchrony measure (kuramoto parameter)
% 	out.meanphase = mean phase of oscillators at t
% 	out.seed = random seed
% 
% parameter value pairs:
% 
% ('radius', 4)
% ('dt', 1e-3) 	
% ('simtime',1) % in seconds
% ('omega_mean', 10) 
% ('plotme', 1) 
% ('noise', 0) 
% ('connectivity', []) 
% connectivity  = 'inverse dist';
% connectivity  = 'euclidean';
% connectivity = 'inverse dist';
% connectivity = 'chebychev';
 % 
 % kuramotoSheet([50 50],105)
 % 
 % author: m@negrello.org
 % all rights to kuramoto and mathworks, all wrongs are mine ;D

plasticity = 0;
gpu = 0;

anim = 1; makemovie = 1;

% [=================================================================]
%  parse inputs
% [=================================================================]

	p = inputParser;
	p.addRequired('networksize')  
	p.addRequired('scaling')  

	p.addParameter('radius', 3)
	p.addParameter('dt', 1e-3) 	
	p.addParameter('time',1) % in seconds
	p.addParameter('omega_mean', 7) % in Hz 
	p.addParameter('init_cond', []) % in Hz 
	p.addParameter('omega_std', 2) % in Hz 
	p.addParameter('oscillators', []) 
	p.addParameter('plotme', 1) 
	p.addParameter('noise', 0) 
	p.addParameter('connectivity', 'euclidean')  % adjacency matrix
	p.addParameter('clusterize', [0 0 0 0],@isvector)
	p.addParameter('seed', 0)

	p.parse(varargin{:});

	netsize = p.Results.networksize;
	scaling = p.Results.scaling;
	radius  = p.Results.radius;
	connectivity = p.Results.connectivity;
	dt = p.Results.dt;
	simtime = p.Results.time;
	omega_mean = p.Results.omega_mean;
	omega_std = p.Results.omega_std;
	plotme = p.Results.plotme;
	noise = p.Results.noise;
	clusterize = p.Results.clusterize;
	seed = p.Results.seed;
	init_cond = p.Results.init_cond;
	oscillators = p.Results.oscillators;

	N = netsize(1);
	M = netsize(2);
	NO = prod(netsize);
	idx = ones(1,NO);

% [=================================================================]
%  randomize oscillator intrinsic frequencies
% [=================================================================]

rng(seed,'twister')

if oscillators
	omega_i = oscillators;
else
	omega_i = (randn(N*M,1)*omega_std+omega_mean)*2*pi;
end

scale_to_intrinsic_freq = 0;




% [=================================================================]
%  connectivity
% [=================================================================]


switch connectivity
	case 'all to all'
		connectivity = ones(N*M) - eye(N*M);

	case 'chebychev'
        
        [X Y] = meshgrid([1:N],[1:M]);
        X = X(:); Y = Y(:); 

        % # compute adjacency matrix
		connectivity = squareform( pdist([X Y], 'chebychev') <= radius );

	case 'euclidean'

		[X Y] = meshgrid([1:N],[1:M]);
        X = X(:); Y = Y(:); 

        % # compute adjacency matrix
		connectivity = squareform( pdist([X Y], 'euclidean') <= radius );


	case 'inverse dist'
		depth   = [1:N];
        breadth = [1:M];
        
        [X Y] = meshgrid(depth,breadth);
        X = X(:); Y = Y(:); 

        % # compute adjacency matrix
		connectivity = 1./squareform( pdist([X Y], 'euclidean') );



	case 'random'
		% W = (ones(N*M)-eye(N*M) ) .* rand(N*M);


	
	otherwise 
		% we use the matrix passed as a parameter
		connectivity = connectivity

end


% [=================================================================]
%  this
% [=================================================================]


if clusterize(1)
	W = connectivity;

    groupsize = clusterize(2);

    k = round(NO/groupsize);
    [idx] = kmeans([X Y], k);

    ProbCluster = clusterize(3);
    ProbOriginal = clusterize(4);

    cW = zeros(NO);
    for ii = 1:NO
        for jj = 1:NO
            if idx(ii)==idx(jj)
                cW(ii,jj) = 1;
            end
        end
    end

    W = cW.* ( (rand(NO)+(eye(NO))) <= ProbCluster) + W.* ( rand(NO) <= ProbOriginal) ;
    

	    out.stats.clusters = idx;


connectivity = W;
end


% ensure that there are no self connections
connectivity(find(eye(M*N))) = 0;


% [=================================================================]
%  Scale coupling parameter?
% [=================================================================]

% Hu, X., Boccaletti, S., Huang, W., Zhang, X., Liu, Z., Guan, S., & Lai, C.-H. (2014). Exact solution for first-order synchronization transition in a generalized Kuramoto model. Scientific Reports, 4, 7262–6. http://doi.org/10.1038/srep07262


if scale_to_intrinsic_freq
	connectivity = bsxfun(@times, omega_i, connectivity) * scaling / NO;
else
	connectivity = scaling*connectivity;
end


% [=================================================================]
%  create noise
% [=================================================================]

ou_noise = zeros(NO, simtime/dt);
if noise

	mixalpha = .3;
	sig = .05;
	th = 10;
	mu = 0;

	for t = 2:simtime/dt
		ou_noise(:,t) =  ou_noise(:,t-1) +th*(mu-ou_noise(:,t-1))*dt + ...
	        (1-mixalpha)*sig*sqrt(dt)*randn(NO,1) + ...
	         mixalpha*sig*sqrt(dt)*randn*ones(NO,1);
	end
end


% [=================================================================]
%  randomize initial condition
% [=================================================================]
%initial condition (initial phase)
if isempty(init_cond)
	theta_t(:,1) = rand(N*M,1)*2*pi;
else
	% disp('using initial condition (assuming between 0 and 2pi)')
	theta_t(:,1) = init_cond;
end

k = zeros(1,simtime*(1/dt)); 
MP = zeros(simtime*(1/dt));

if gpu
	theta_t = gpuArray(theta_t);
	connectivity = gpuArray(connectivity);
	ou_noise = gpuArray(ou_noise);
	k = gpuArray(k);
end

% [=================================================================]
%  simulate
% [=================================================================]

if plotme; f = figure(100); a(1) = subplot(121);a(2) = subplot(122); end
for t = 2:simtime/dt

	phasedifferences = bsxfun(@minus, theta_t(:,t-1)',theta_t(:,t-1));

	phasedifferences_W = connectivity.*sin(phasedifferences);
	
	summed_sin_diffs = mean(phasedifferences_W,2);

	theta_t(:,t) = theta_t(:,t-1) + dt*( omega_i + summed_sin_diffs  ) + ou_noise(:,t);

	PP(:,t) = sin(mod(theta_t(:,t),2*pi));

	if plasticity
		% potentiation of connections by phase coincidence of pairs

		delta_W = connectivity./connectivity;

	end

	% [=================================================================]
	%  order parameter
	% [=================================================================]
	if ~clusterize(1)
		GP = theta_t(:,t);
		MP = circ_mean(GP)+pi;
		
		k(t) = mean( exp(i*(bsxfun(@minus, GP, MP))));
	else
		for ui = unique(idx)'
			GP = theta_t(find(idx==ui),t);
			MP = circ_mean(GP)+pi;
			
			k(ui,t) = mean( exp(i*(bsxfun(@minus, GP, MP))));
		end
	end

end


% [=================================================================]
%  plots
% [=================================================================]



if plotme
	
	ffff = figure


	subplot(2,2,1)
	plot(linspace(0,simtime, simtime*dt^-1), mod(theta_t,2*pi)')
	ylabel('phase (theta)')
	subplot(2,2,2)
	imagesc(connectivity), colorbar
	title('connectivity')
	subplot(2,2,3)
	plot(linspace(0,simtime, simtime*dt^-1), sin(theta_t'))
	ylabel('sin(theta)')
	xlabel('seconds')
	subplot(2,2,4)
	hist(omega_i/(2*pi),N)
	ylabel('#')
	xlabel('frequency (Hz)')
	
	figure(f)
	axes(a(2))
		line(repmat(linspace(0,simtime, simtime*dt^-1), length(unique(idx)), 1)', [abs(k)]')


		title('kuramoto parameter')
		xlabel('time (ms)')
		ylim([0 1])
	
	[XX YY] = meshgrid([1:M],[1:N]);
	XX = XX(:); YY = YY(:);
	

	if exist('linspecer')
		LSpec = linspecer(length(unique(idx)))
		set(ffff, 'colormap',   LSpec);
		set(a(1), 'colororder', LSpec);
		set(a(2), 'colororder', LSpec);
	else
		set(ffff, 'colormap',   jet(length(unique(idx))));
		set(a(1), 'colororder', jet(length(unique(idx))));
		set(a(2), 'colororder', jet(length(unique(idx))));
	end

	
	while anim
		
		for t = 2:simtime/dt
			
			SS = reshape(PP(:,t),N,M);
			axes(a(1)); 			cla
			% imagesc(reshape(theta_t(:,t),N,M))
			mesh(SS); hold on
			scatter3(XX, YY ,PP(:,t),60,LSpec(idx,:),'filled')
			title('phase')
			caxis([0 2*pi])
			zlim([-3 3])

			drawnow
			if makemovie
				MOV(t) = getframe(ffff);
			end

		end
		anim = input(['repeat? ; 0 - no, 1 - yes \n'  ])
	end

end


% [================================================]
%  Write out
% [================================================]


if gpu
	theta_t = gather(theta_t);
	connectivity = gather(connectivity);
	ou_noise = gather(ou_noise);
	k = gather(k);
end

out.state = PP;
out.phase = theta_t;
out.parameters = p.Results;
out.oscillators = omega_i/(2*pi);
out.orderparameter = abs(k);
out.meanphase = MP;
out.seed = seed;
 if makemovie && plotme ; out.movie = MOV; end
% out.all =  sin(mod(theta_t,2*pi));
