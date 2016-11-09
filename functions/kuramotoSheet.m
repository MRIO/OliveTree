function out = kuramotoSheet(varargin)
% inputs:
% 	networksize, [N M]
% 	scaling for coupling - K
% 
% parameter value pairs
% 
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

plasticity = 1;


% [=================================================================]
%  parse inputs
% [=================================================================]

	p = inputParser;
	p.addRequired('networksize')  
	p.addRequired('scaling')  

	p.addParamValue('radius', 3)
	p.addParamValue('dt', 1e-3) 	
	p.addParamValue('simtime',1) % in seconds
	p.addParamValue('omega_mean', 10) % in Hz 
	p.addParamValue('omega_std', 2) % in Hz 
	p.addParamValue('plotme', 1) 
	p.addParamValue('noise', 0) 
	p.addParamValue('connectivity', 'euclidean')  % adjacency matrix

	p.parse(varargin{:});

	netsize = p.Results.networksize;
	scaling = p.Results.scaling;
	radius  = p.Results.radius;
	connectivity = p.Results.connectivity;
	dt = p.Results.dt;
	simtime = p.Results.simtime;
	omega_mean = p.Results.omega_mean;
	omega_std = p.Results.omega_std;
	plotme = p.Results.plotme;
	noise = p.Results.noise;
	
	N = netsize(1);
	M = netsize(2);

% [=================================================================]
%  randomize oscillator intrinsic frequencies
% [=================================================================]

omega_i = (randn(N*M,1)*omega_std+omega_mean)*2*pi;



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

% ensure that there are no self connections
connectivity(find(eye(M*N))) = 0;

NO = prod(netsize);
connectivity = (scaling)*(connectivity);



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
theta_t(:,1) = rand(N*M,1)*2*pi;
k = zeros(simtime*(1/dt)); 
MP = zeros(simtime*(1/dt));

% [=================================================================]
%  simulate
% [=================================================================]

if plotme; f = figure(100); a(1) = subplot(121);a(2) = subplot(122); end
for t = 2:simtime/dt

	phasedifferences = bsxfun(@minus, theta_t(:,t-1)',theta_t(:,t-1));

	phasedifferences_W = connectivity.*sin(phasedifferences);
	
	summed_sin_diffs = mean(phasedifferences_W,2);

	theta_t(:,t) = theta_t(:,t-1) + dt*( omega_i + summed_sin_diffs  ) + ou_noise(:,t);

	if plasticity
		% potentiation of connections by phase coincidence of pairs

		delta_W = connectivity./connectivity;

	end

	% [=================================================================]
	%  order parameter
	% [=================================================================]
		GP = theta_t(:,t);
		MP = circ_mean(GP)+pi;

		k(t) = mean( exp(i*(bsxfun(@minus, GP, MP))));
	
	
	if 0


		axes(a(1))
			% imagesc(reshape(theta_t(:,t),N,M))
			mesh(reshape(sin(mod(theta_t(:,t),2*pi)),N,M))
			caxis([0 2*pi])
			zlim([-pi pi])
		axes(a(2))
			line([t-1 t] , [abs(k(t-1)) abs(k(t))])
			

			drawnow


			
	end

end


% [=================================================================]
%  plots
% [=================================================================]


out.phase = theta_t;
out.parameters = p.Results;
out.oscillators = omega_i/(2*pi);
out.orderparameter = abs(k);
out.meanphase = MP;


if plotme
	
	figure

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
		line([1:t] , [abs(k)])
	repeat = 1;
	while 1
		
		for t = 2:simtime/dt
			axes(a(1))
			% imagesc(reshape(theta_t(:,t),N,M))
			mesh(reshape(sin(mod(theta_t(:,t),2*pi)),N,M))
			caxis([0 2*pi])
			zlim([-pi pi])
			
			

			drawnow


		end
		repeat = input(['repeat? ; 0 - no, 1 - yes \n'  ])
	end


end


