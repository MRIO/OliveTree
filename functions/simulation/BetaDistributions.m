function out = BetaDistributions(varargin)
% Family of beta distributions
% 
% Parameter 'alpha', 5
% Parameter 'beta', 2
% Parameter 'no_draws', 100
% Parameter 'scale', 1
% Parameter 'support',[0.4 1.2]
% Parameter 'resolution', 4



p = inputParser;
p.addParamValue('alpha', 5)
p.addParamValue('beta', 2)
p.addParamValue('no_draws', 100)

p.addParamValue('scale', 1)
p.addParamValue('support',[0.4 1.2])


p.addParamValue('plot_distributions', true)

p.parse(varargin{:});

p1 = p.Results.alpha;
p2 = p.Results.beta;

no_draws = p.Results.no_draws;
plot_distributions = p.Results.plot_distributions';
scale = p.Results.scale;

sup = p.Results.support;



int_v1 = [1:length(p1)];
int_v2 = [1:length(p2)];

nbins = linspace(sup(1), sup(end),20);

for a = int_v1
    for b = int_v2
        rng default;

        D =  betarnd(p1(a),p2(b),no_draws,1)';
        
        draws{a,b} = (D*diff(sup) + sup(1))*scale;
        


        parameters{a,b} = [p1(a) p2(b)];
        hist_distr{a,b} = histc(draws{a,b}, nbins);  
            if plot_distributions
                subplot(length(p1),length(p2), sub2ind([length(p1) length(p2)], a,b))
                bar(nbins, hist_distr{a,b},'edgecolor', [0 0 0], 'facecolor', [0 0 0])%;,'edgecolor', [0 0 0], 'facecolor', [0 0 0])
                t1 = text(min(nbins), max(hist_distr{a,b}*0.95), {['alpha: ' num2str(p1(a))]},'verticalalignment', 'top');
                t2 = text(min(nbins), max(hist_distr{a,b}*0.75), {['beta: '  num2str(p2(b))]},'verticalalignment', 'top');
                t3 = text(min(nbins), max(hist_distr{a,b}*0.55), {['scale: ' num2str(scale)]},'verticalalignment', 'top');
                
                line([sup(1):0.01:sup(2)], betapdf([sup(1):0.01:sup(2)]*max(hist_distr{a,b}), p1(a),p2(b)),'color','g');


                
                xlim([sup(1) sup(2)]) 
                ylim([0  max(hist_distr{a,b})])
        end
    end
end

out.sampleDraws = draws;
out.parameters = parameters;
