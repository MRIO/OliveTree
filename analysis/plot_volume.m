function plot_volume(V_soma_unwrapped, coords, ttt)

fname = ['volume_' num2str(now)];
savemovie = 1;
frames2file = 0;
plotdiffact = 0;
if plotdiffact
	 clims = [0 2*pi];
else
	clims = [-68 -40];
end

if isempty(ttt)
	ttt = 1:size(V_soma_unwrapped,2)
end

gksz = 7;
g3d3 = gaussKernel3d(.2,  gksz, ceil(gksz/2)); g3d3 = g3d3/sum(g3d3(:)); g3d3(g3d3<.005) = 0; g3d3(g3d3>0.005)=1; %g3d3 = g3d3/sum(g3d3(:)); 
% figure
% cla
% vol3d('cdata',g3d3)
% colorbar

coarseness = 4;

% [=================================================================]
%  if saving movie
% [=================================================================]

if savemovie
	
		vidObj = VideoWriter(fname,'MPEG-4');
		
		set(vidObj,'FrameRate', 25)

	open(vidObj);
end




fig_volume = figure('color', [1 1 1]);
ax_volume = axes;
colormap(ax_volume, jet(128));
colorbar
for tt = 1:size(V_soma_unwrapped,2)
	VVVV = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), V_soma_unwrapped(:,tt));
	NNNN = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), 1);
	NNNN(NNNN==0)=1;
	VVVV = VVVV./NNNN;
	
	% CCCC = interp3(VVVV, 3);
	% CCCC = convn(VVVV, g3d3, 'same');
	% VVVV = CCCC;
	CCCC = imerode(VVVV,g3d3);
	% CCCC = imdilate(VVVV,g3d3);

	set(0,'CurrentFigure',fig_volume);
    % set(fig1,'CurrentAxes',a(3));

    if plotdiffact
    	AAA = abs(CCCC);
		AAA(AAA<0)  = 0;
		AAA(AAA==0) = 0;
		AAA(AAA>=0)  = 1;
		AAA =AAA*.25;

    else
	    
	    AAA = zeros(size(CCCC));
		AAA(CCCC>-35) = 0;
		AAA(CCCC<=-35 & CCCC>=-70) = 1;
		AAA(CCCC<=-70) = .8;
		AAA = AAA*.25;
		% AAA = ~AAA*.25

	end

	cla(ax_volume)
	set(ax_volume, 'clim',clims)
	
	vol3d('cdata',CCCC, 'Alpha', AAA , 'texture','3D');


	if tt==1;view(3) ; view(-22,-56.4);; axis tight;  daspect([1 1 1]); end %axis off; end

	title([num2str(ttt(tt)) 'ms'])
	drawnow

	if savemovie
		writeVideo(vidObj, getframe(fig_volume))
	end

	if frames2file
		saveallfigs('prefix', ['volact@' num2str(tt) ],'format', 'png','savefig',0,'style','12x6')
	end

end