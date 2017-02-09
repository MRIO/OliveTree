function out = connectivity_statistics(in)

W = in.W;


out.stats.stdW = std(W(W~=0));
out.stats.connections = sum(W>0);
out.stats.meanweight  = mean(W(find(W>0)));

if exist('clustering_coef_wd')
    out.stats.clustercoeff.wd = clustering_coef_wd(W);
    out.stats.clustercoeff.bu = clustering_coef_bu(W);
    out.stats.degree = degrees_und(W);
    % out.stats.
else
    out.stats.clustercoeff = 'WARNING: did not find connectivity toolbox';
end

try 
	out.stats.clusters = in.stats.clusters;
catch
	disp('no clusters in struct')
end