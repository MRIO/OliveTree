% animate_volume.m

function animate_volume_hilbert(sim,frames, savemovie, varargin)

backgroundcolor = 'dark';
backgroundcolor = 'light';

plotvol = 1; 
trigger = 1;
showcolorbar = 0;
scatterit = 1;
printframes = true;

if nargin==4
	interpv = varargin{1};
else
	interpv = 1;
end


try;simtime = sim.time;catch; simtime = length(sim.networkHistory.V_soma);end



if savemovie
	printframes = false;
	disp('warning: cant print frames when saving movie')
end

if ~isempty(sim.perturbation) % backward compatibility

	try
		pert_mask			= sim.perturbation.mask{trigger};
		if iscell(sim.perturbation.triggers)
			perturbation_triggers = cell2mat(sim.perturbation.triggers);
		else
			perturbation_triggers = sim.perturbation.triggers;
		end
	catch
		pert_mask = [];
		perturbation_triggers = [];
	end

else	
	pert_mask = [];
	perturbation_triggers = [];
end

try
	if isfield(sim.W, 'W')
		coords = sim.W.coords;
	end
catch
	disp('could not read W')
end


if ~isfield(sim, 'networksize')
	depth = 5;
	breadth = 10;
	height = 20;
else
	% depth = sim.networksize(1); 
	% breadth = sim.networksize(2);
	% height = sim.networksize(3);

	% why is this mapping WEIRD?
	depth = sim.networksize(2); 
	breadth = sim.networksize(1);
	height = sim.networksize(3);
end



% SIM3D.     networkHistory: [1x1 struct]
%                rows: 50v
%             columns: 20
%           lastState: [1000x23 double]
%            duration: 1000
%                  dt: 0.0500
%               g_CaL: [1000x1 double]
%               g_Gap: [1000x1000 double]
%                time: 1000
%         timeElapsed: 310.9454
%        initialState: [1000x23 double]
%     adjacencyMatrix: [1000x1000 double]
%           condition: [1x1 struct]

%               __                     __       _                    __  _       _ __       
%  _   ______  / /_  ______ ___  ___  / /______(_)____   ____ ______/ /_(_)   __(_) /___  __
% | | / / __ \/ / / / / __ `__ \/ _ \/ __/ ___/ / ___/  / __ `/ ___/ __/ / | / / / __/ / / /
% | |/ / /_/ / / /_/ / / / / / /  __/ /_/ /  / / /__   / /_/ / /__/ /_/ /| |/ / / /_/ /_/ / 
% |___/\____/_/\__,_/_/ /_/ /_/\___/\__/_/  /_/\___/   \__,_/\___/\__/_/ |___/_/\__/\__, /  
%                                                                                  /____/   



% f = figure('position', [440 37 481 761]);



% switch backgroundcolor
% 	case 'light'
% 		try
% 			% load activity_cmap_hot
% 			load activity_cmap_hot
% 			cm = cmap;
% 		catch
% 			cm = hot(40); 
% 		end
% 		try
% 		cm  = flipud(cbrewer('div', 'RdBu', 40));
% 		catch
% 		end

% 		% !cm = bone(40);
% 		set(f,'colormap', cm,'color', [1 1 1])

% 	case 'dark'
% 		try
% 			% load activity_cmap_hot
% 			load activity_cmap_jet
% 			cm = cmap;
% 		catch
% 			cm = jet(40); 
% 		end
		
% 		set(f,'colormap', cm,'color', [.2 .2 .2])
% end

if savemovie
	% fname = [num2str(rows) 'x' num2str(columns) '_.avi']
		try
			vidObj = VideoWriter('volume','MPEG-4');
		catch
			warning('No MPEG-4 encoder in this crappy OS.')
			vidObj = VideoWriter('volume');
		end

		vidObj.FrameRate = 100;
	open(vidObj);
end

	


% interpolate and plot the volume


if showcolorbar
	bla = colorbar;
	switch backgroundcolor
		case 'light'
			set(bla,'color',[0 0 0])
		case 'dark'
			set(bla,'color',[1 1 1])
	end
end


switch backgroundcolor
	case 'light'
		if interpv
			% alphamap([0.1:0.01:.1 .2: .05:.9 ]);
			alphamap([.2: .05:.9 ]);
		else
			alphamap([0.2:0.01:.3 0.5:0.05:.9 ]);
		end
		
		% alphamap([0:0.01:0.7 0.3:0.05:.7 ]);
		% alphamap(0.1);
		% alphamap([0.5:-0.01:0 0:0.01:.5 ]);
		% alphamap([0.3:0.01:.7 ]);
		% alphamap(1+ 1./(1-exp(linspace(-2,2,100))))

	case 'dark'
		if interpv
			alphamap([0.05:0.01:0.1 0.5:.1:.7 ]);
			% alphamap([0:0.01:.1 .2: .05:.9 ]);
		else
			alphamap([0.2:0.01:.3 0.5:0.05:.9 ]);
		end
		


		% alphamap(0.2)
end


% do it for frames
if isempty(frames)
	frames = 1:simtime;
end

HH = hilbert_of_membranepotential(sim.networkHistory.V_soma(:,frames(1):frames(end)));
HHH = HH.hilbert;

figure; imagesc(HHH)


			
			gksz = 7;
			g3d3 = gaussKernel3d(.2,  gksz, ceil(gksz/2)); g3d3 = g3d3/sum(g3d3(:)); g3d3(g3d3<.005) = 0; g3d3(g3d3>0.005)=1; %g3d3 = g3d3/sum(g3d3(:)); 
			% figure
			% cla
			% vol3d('cdata',g3d3)
			% colorbar

			coarseness = 4;


			fig_volume = figure('color', [1 1 1]);
			ax_volume = axes;
			colormap(ax_volume, linspecer(128));
			% colorbar

			for tt = frames-frames(1)+1;
				VVVV = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), HHH(:,tt));
				NNNN = accumarray( round([coords(:,1), coords(:,2), coords(:,3)]/coarseness+1), 1);
				NNNN(NNNN==0)=1;
				VVVV = VVVV./NNNN;
				

				% CCCC = interp3(VVVV, 3);
				% CCCC = convn(VVVV, g3d3, 'same');
				% VVVV = CCCC;
				CCCC = imdilate(VVVV,g3d3);

				set(0,'CurrentFigure',fig_volume);
			    % set(fig1,'CurrentAxes',a(3));


			    AAA = CCCC>.1;

				cla(ax_volume)
				set(ax_volume, 'clim',[0 2*pi])
				
				vol3d('cdata',CCCC, 'Alpha', AAA*.25 , 'texture','3D');


				if tt==1;view(3) ; view(-22,-56.4);axis off; axis tight;  daspect([1 1 1]); axis equal; end
				title([num2str(frames(1)+ tt) 'ms'])
				drawnow
				if savemovie
					writeVideo(vidObj, getframe(fig_volume))
				end

				if printframes
						prefix = 'netstate@';
						style = '4x4';
						fname = [prefix '_' num2str(frames(1)+tt) '.png'];
						snam=style;
						s=hgexport('readstyle',snam);
					    s.Format = 'png';
					    hgexport(fig_volume,fname,s);
				end

				

			end
			if savemovie
				close(vidObj)
			end



