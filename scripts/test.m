
multiwindowxcorr = 0;
computeselectedxcorr = 1; sims2p = 1;
% 

simulation = simresults;


netsize = simulation{1}.networksize;

if computeselectedxcorr
	
	XnoAC  = [];
	lag = 300;
	centerwin = 25;
	noneur = 30;
	xcorrwin = 1:50000;
	nwins = 5; windur = 10000;



	c = 0;
	for simcount = sims2p
		c = c+1;
		

		VS = simulation{simcount}.networkHistory.V_soma;
		W  = simulation{simcount}.networkParameters.connectivityMatrix;
		% NS = simulation{simcount}.noiseapplied;

			spks = spikedetect(simulation{simcount});

			VSB = single(VS);
			VSB(VS<-10)=0;
			VSB(VS>=-10)=1;

			VSB = spks.binaryspikes;


			% take N most active neurons
			[V sorted] = sort(sum(VSB(:,xcorrwin)'),2,'descend');

			selectedneurons = sorted(1:noneur);

			N(c) = length(selectedneurons);

			XC{c} = xcorr(VSB( selectedneurons ,xcorrwin)',lag);

			pairs = setdiff([1:N(c)^2],1:N(c):N(c)^2);
			pairs = find(triu(ones(N(c))) - eye(N(c)) ) ;

			xcorrset = XC{c}(:,pairs)';
			asym{c} = sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)');


			XnoAC(:,c) = sum(xcorrset);





		% sorted by asymmetry
		[V O] = sort(asym{c});

		figure
		imagesc([-lag:lag], [1:length(pairs)], xcorrset(O,:))
		line([0 0], [0 length(pairs)],'color', 'w')
		% title({num2str(Ptable(simcount,:)) ; 'sc gp nc na fr'})
		ylabel('pairs')
		xlabel('lag')
		
		% Ptable = [sc gp nc na fr ];


        depth   = [1:netsize(1)];
        breadth = [1:netsize(2)];
        height  = [1:netsize(3)];

        % # compute adjacency matrix
        [X Y Z] = meshgrid(depth,breadth,height);
        X = X(:); Y = Y(:); Z = Z(:);

        idx = zeros(length(X),1);
        idx(selectedneurons) = 1;

        % plotnetstruct(W,X,Y,Z,idx)


        plotnetstruct(W,X,Y,Z,sum(VSB'))

	end

	figure
	plot([-lag: lag],XnoAC/ ((N(c)^2-N(c))) )
	legend(num2str(sims2p'))
	legend({'WT' ; 'Cx36' })
	axis tight
	xlabel('ms')
	ylabel('correlation (coeff)')
	xlabel('aggregate correlation')
	

end


	




if multiwindowxcorr 
% [12 16 19 24]
	
	XnoAC  = [];
	lag = 300;
	centerwin = 25;
	noneur = 30;

	xcorrwin = 1:10000:100000;

	c = 0;
	sims2p = 11;
	for simcount = sims2p
		c = c+1;
		

		VS = simulation{simcount}.networkHistory.V_soma;
		W  = simulation{simcount}.networkParameters.connectivityMatrix;
		% NS = simulation{simcount}.noiseapplied;

			spks = spikedetect(simulation{simcount});

			VSB = single(VS);
			VSB(VS<-10)=0;
			VSB(VS>=-10)=1;

			VSB = spks.binaryspikes;

			% take N most active neurons
			[V sorted] = sort(sum(VSB(:,xcorrwin)'),2,'descend');

			selectedneurons = sorted(1:noneur);

			N(c) = length(selectedneurons);

			nw = 0;
			for nw = 1:length(xcorrwin)-1;

				XC{c} = xcorr(VSB( selectedneurons ,xcorrwin(nw):xcorrwin(nw+1))',lag);

				pairs = setdiff([1:N(c)^2],1:N(c):N(c)^2);
				pairs = find(triu(ones(N(c))) - eye(N(c)) ) ;

				xcorrset = XC{c}(:,pairs)';
				asym{c} = sum(xcorrset(:,lag-centerwin:lag)')-sum(xcorrset(:,lag+1:lag+1+centerwin)');


				XnoAC(:,nw) = sum(xcorrset);

			end

	figure
	plot([-lag: lag],XnoAC/ ((N(c)^2-N(c))) )
	legend(num2str(sims2p'))
	% legend({'WT' ; 'Cx36' })
	axis tight
	xlabel('ms')
	ylabel('correlation (coeff)')
	xlabel('aggregate correlation')
	



		% % sorted by asymmetry
		% [V O] = sort(asym{c});

		% figure
		% imagesc([-lag:lag], [1:length(pairs)], xcorrset(O,:))
		% line([0 0], [0 length(pairs)],'color', 'w')
		% title({num2str(Ptable(simcount,:)) ; 'sc gp nc na fr'})
		% ylabel('pairs')
		% xlabel('lag')
		
		% % Ptable = [sc gp nc na fr ];


  %       depth   = [1:netsize(1)];
  %       breadth = [1:netsize(2)];
  %       height  = [1:netsize(3)];

  %       % # compute adjacency matrix
  %       [X Y Z] = meshgrid(depth,breadth,height);
  %       X = X(:); Y = Y(:); Z = Z(:);

  %       idx = zeros(length(X),1);
  %       idx(selectedneurons) = 1;

  %       plotnetstruct(W,X,Y,Z,idx)

	end



end


	





% legend({'.25 .1' ; '.25 .4' ; '0 .1' ; '0 .4'})
% legend(num2str([12 16 19 24]'))



	% Ptable =

	%     1.0000    0.0500         0    1.0000    0.0155
	%     2.0000    0.0500         0    1.2000    0.1165
	%     3.0000    0.0500    0.1000    1.0000    0.0114
	%     4.0000    0.0500    0.1000    1.2000    0.0659
	%     5.0000    0.0500    0.2000    1.0000    0.0234
	%     6.0000    0.0500    0.2000    1.2000    0.0885
	%     7.0000    0.0500    0.4000    1.0000    0.1510
	%     8.0000    0.0500    0.4000    1.2000    0.3759
	%     9.0000    0.0250         0    1.0000    0.4264
	%    10.0000    0.0250         0    1.2000    1.1095
	%    11.0000    0.0250    0.1000    1.0000    0.2401
	%    12.0000    0.0250    0.1000    1.2000    0.7203
	%    13.0000    0.0250    0.2000    1.0000    0.1792
	%%%    14.0000    0.0250    0.2000    1.2000    0.5219
	%    15.0000    0.0250    0.4000    1.0000    0.2953
	%    16.0000    0.0250    0.4000    1.2000    0.6309
	%    17.0000    0.0010         0    1.0000    2.3214
	%    18.0000    0.0010         0    1.2000    3.0339
	%    19.0000    0.0010    0.1000    1.0000    1.9656
	%    20.0000    0.0010    0.1000    1.2000    2.6470
	%    21.0000    0.0010    0.2000    1.0000    1.6602
	%%%    22.0000    0.0010    0.2000    1.2000    2.2909
	%    23.0000    0.0010    0.4000    1.0000    1.2583
	%    24.0000    0.0010    0.4000    1.2000    1.8228
