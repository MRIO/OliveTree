% connectivity_tests.m

conntype = '3d_chebychev';
netsize = [4 5 1];
radius = 2;
plotthis = 1;
meanconn = 10;
scaling = 1;
maxiter = [];
clusterize = 0;
somatapositions = [];

% wo/randomize
out = createW(conntype, netsize, radius, scaling, 0, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)

% w/randomize
out = createW(conntype, netsize, radius, scaling, 1, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)

% w/randomize w/symmetrize
out = createW(conntype, netsize, radius, scaling, 1, plotthis, maxiter, meanconn, somatapositions, 1, clusterize)

% w/randomize w/symmetrize w/normalize leaks
out = createW(conntype, netsize, radius, scaling, 1, plotthis, maxiter, meanconn, somatapositions, 1, clusterize,1)
out = createW(conntype, netsize, radius, .1, 1, plotthis, maxiter, meanconn, somatapositions, 1, clusterize,1)


imagesc(out.W./out.W')
colorbar


Wbrick 		= createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], 0);
Wcluster50  = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], 0, [1 50 .5 0], 1);
Wallcluster = createW('3d_chebychev', netsize, 3, 1, 1, 1, [], 8, [], 0, [1 10 .9 .01]);
