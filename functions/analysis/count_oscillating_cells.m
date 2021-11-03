function out = count_oscillating_cells(sims, tt, log_amp_threshold)

transients = 500;


z = sims.networkHistory.V_soma;
ttt = size(sims.networkHistory.V_soma,2);

if isempty(tt)
    tt = transients:size(sims.networkHistory.V_soma,2);
end
    

% Remove the beginning for better processing
z1 = double(z(:, tt)');
% Data without a dc component
z2 = bsxfun(@minus, z1, mean(z1,1));

% Here we first remove the non-oscllating cells.
z3 = hilbert(z2);
logamp = log(mean(abs(z3)));
histog = histogram(logamp, [-7:.1:3]);

% Two peaks at the end are from non-oscillating cells. We remove them from our data...

 
igood = find(logamp>= log_amp_threshold);
ibad = find(logamp< log_amp_threshold);
% 
% figure
% z4 = z2(:,igood);
% subplot(2,1,1)
% plot(z2(:,igood))
% subplot(2,1,1)
% plot(z2(:,ibad))


out.histogram = histog;
out.oscillating= igood;
out.proportion = length(igood) / (length(igood)+length(ibad));