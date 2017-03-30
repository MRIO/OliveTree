function out  = hilbert_of_membranepotential(V)
% function out  = hilbert_of_membranepotential(V)
% 
%  calculates hilbert transform of membrane potential using DAMOCO toolbox (required)
% 
% 

plotme = 0;

if exist('co_hilbproto')
	V(V>-30) = -30; % crop spikes

	for n = 1:size(V,1) % num neurons
		warning off
		out.protophase(n,:) = co_hilbproto(V(n,:),0,0,0,0);
		warning on
        try
            [out.phi(n,:) out.arg(n,:) out.sig(n,:)] = co_fbtrT(out.protophase(n,:));
        catch


           out.phi(n,:) = zeros(1,size(V,2));
           out.arg(n,:) = 0;
           out.sig(n,:) = 0;

           continue

        end
        
		
	end
	out.phidot = diff(unwrap(out.phi)') /(2*pi) * 1e3; 

	mean_phase = circ_mean(out.phi);
	out.order_parameter = mean( exp(i*(bsxfun(@minus, out.phi, mean_phase))));
	out.hilbert = out.phi;
	out.phi = [];


else	% if damoco toolbox not present, use old method
		% crop complex spikes

		disp('not using DAMOCO method for instantaneous frequency estimation')
		V(V>-30) = -30;
		N = -V;

		deno = (max(N,[],2)-min(N,[],2))';
		scaled = bsxfun(@minus, N, min(N,[],2));
		scaled = bsxfun(@rdivide, scaled, deno');
		scaled = scaled *2 -1;


		for n = 1:size(V,1) % num neurons
			out.hilbert(n,:) = angle(hilbert(scaled(n,:)));
		end

		mean_phase = circ_mean(out.hilbert);

		out.order_parameter = mean( exp(i*(bsxfun(@minus, out.hilbert, mean_phase))));

end






if plotme
	imagesc(hist(out.phi,100))
end