function add_x_scalebar(fhandle, barlength)

	xl = get(fhandle, 'xlim');
	yl = get(fhandle, 'ylim');

	line([xl(2) xl(2)-barlength],  [yl(1) yl(1)], 'linewidth', 2 )

