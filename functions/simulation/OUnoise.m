function O = OUnoise(varargin)
% function O = OUnoise(varargin)
% 
% produces an OU noise trace with the given parameters]
%   p.addParamValue('mu', 0)
% 	p.addParamValue('sig', .7)
% 	p.addParamValue('thetas', [1/20])
% 	p.addParamValue('simtime', 2000)
% 	p.addParamValue('dt', 0.02)
% 	p.addParamValue('plotme', 1)
% 	p.addParamValue('seed', 0)

	p = inputParser;
	
	p.addParamValue('mu', 0)
	p.addParamValue('sig', .7)
	p.addParamValue('thetas', [1/20])
	p.addParamValue('simtime', 2000)
	p.addParamValue('dt', 0.02)
	p.addParamValue('plotme', 1)
	p.addParamValue('seed', 0)
	
	p.parse(varargin{:});

	mu = p.Results.mu;
	sig = p.Results.sig;
	thetas = p.Results.thetas;
	simtime = p.Results.simtime;
	dt = p.Results.dt;
	plotme = p.Results.plotme;
	seed = p.Results.seed;
	

newstate = @(state, tau, mu, sig) state + (1/tau)*(mu-state)*dt + sig*sqrt(dt)*randn;

taus = 1./thetas;
state=zeros(length(taus),simtime);

state(length(thetas),1) = mu; n =0;

set(0, 'defaultaxescolororder', linspecer(length(taus)))



for tau = taus
	n = n+1;
	rng(seed)
	for t = 1:simtime*(1/dt)
	% curr_noise(n,t+1) =	curr_noise(n,t)*exp(-delta/tau) + sqrt((1*tau*0.5)*(1-(exp(-delta/tau))^2))*randn;
	
	state(n,t+1) = newstate(state(n,t), tau, mu, sig);

	end
end


if plotme
	figure
		subplot(2,2,[1 2])
		plot([0:simtime*(1/dt)]*dt, state')
		legend(num2str(taus'))
		title('ohrstein uhlenbeck process with 0 mean and 1 std')


		subplot(2,2,3)
		for i = 1:n
			[xc(i,:) tau] = xcorr(state(i,:),'coeff');
			noiseintensity_empirical(i) = sum(xc(i,:));
			noiseintensity_calculated(i) = taus(i)*sig^2;

		end
		plot(tau*dt, xc)
		xlim([-1000 1000])
		title('ohrstein uhlenbeck process with 0 mean and 1 std')

		subplot(2,2,4)
		sup = [-30:1:30];
		for i = 1:n
			nhist(:,i) = hist(state(i,:),sup);
		end
		plot(sup, nhist/length(state))
		title('probability of current applied /ms')
		xlabel('pA')
		ylabel('probability')

end


O = state;


% [=================================================================]
%   summed processes
% [=================================================================]
% nosources = 10;
% curr_noise=zeros(nosources,simtime);
% th = 2; n = 0; tau = 5; th = 1/5; tau = 10 ,th = 1/10; 
% for noisesource = 1:nosources
% 	n = n+1;
% 	for t = 1:simtime
% 	curr_noise(n,t+1) = curr_noise(n,t) + (1/tau)*(mu-curr_noise(n,t))*delta + sig*sqrt(delta)*randn;
% 	end
% end
