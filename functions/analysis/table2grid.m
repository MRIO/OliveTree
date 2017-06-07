function out = table2grid(dimnames, measured, STable)

% [================================================================================================]
%  					function out = table2grid(dimnames, measured, table)
% [================================================================================================]
% 
% ============================description=========================
% accumulates a 3d grid from a 2d table such that for every unique 
% value combination of 'dimnames' a 3d grid array entry is equal 
% to the 'measured' field in the input table.
% 
% 
% 
% ============================output==============================
% out.dimnames = dimnames;
% out.measured = measured;
% out.results  = X;
% out.ticks.X = ui;
% out.ticks.Y = uj;
% out.ticks.Z = uk;
% ============================mnegrello==============================


ui = table2array(unique(STable(:,dimnames{1})));
uj = table2array(unique(STable(:,dimnames{2})));
uk = table2array(unique(STable(:,dimnames{3})));

D1 = table2array(STable(:,dimnames{1}));
D2 = table2array(STable(:,dimnames{2}));
D3 = table2array(STable(:,dimnames{3}));

for ii = 1:length(ui)

	for jj = 1:length(uj)

		for kk = 1:length(uk)

			row =    D1== ui(ii) ...
				  &	 D2== uj(jj) ...
				  &  D3== uk(kk);

				V = table2array(STable( row , measured ));

				try
						R(ii,jj,kk) = V;
				catch
					warning('check dimensions in analysis_script')
						R(ii,jj,kk) = V(2);
					dimnames
					ui
					uj
					uk
					% keyboard
					
			end

		end
	end
end

out.dimnames = dimnames;
out.measured = measured;
out.results  = R;
out.ticks.X = ui;
out.ticks.Y = uj;
out.ticks.Z = uk;


