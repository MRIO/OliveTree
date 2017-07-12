function addScalebar(fhandle, barlength)

	xl = get(fhandle, 'xlim');
	yl = get(fhandle, 'ylim');

	line([xl(2) xl(2)-barlength(1)],  [yl(1) yl(1)], 'linewidth', 2 )
	line([xl(2) xl(2)],  [yl(1) yl(1)+barlength(2)], 'linewidth', 2 )

