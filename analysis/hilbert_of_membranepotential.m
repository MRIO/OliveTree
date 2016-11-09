% calc_order_param_from_V.m
function out  = hilbert_of_membranepotential(V)

plotme = 0;

if exist('co_hilbproto')
	V(V>-20) = -20; % crop spikes

	for n = 1:size(V,1) % num neurons
		warning off
		out.protophase(n,:) = co_hilbproto(V(n,:),0,0,0,0);
		[out.phi(n,:) out.arg(n,:) out.sig(n,:)] = co_fbtrT(out.protophase(n,:));
		warning on
	end

	mean_phase = circ_mean(out.phi);
	out.order_parameter = mean( exp(i*(bsxfun(@minus, out.phi, mean_phase))));
	out.hilbert = out.phi;
	out.phi = [];

else	% if damoco toolbox not present, use old method
		% crop complex spikes
		V(V>-40) = -40;
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