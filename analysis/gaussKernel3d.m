function g3d = gaussKernel3d(a, n, mu)
	[xx yy zz] = meshgrid([1:n],[1:n],[1:n]);
	gausskernel3d = @(x, y, z)(exp(-((x-mu)^2 + (y-mu)^2 + (z-mu)^2)*a));
	g3d = arrayfun(gausskernel3d, xx, yy, zz);
	g3d = g3d/sum(sum(sum(g3d)));

