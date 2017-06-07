% animate_volume.m

function animate_volume(sim,frames, savemovie, varargin)

backgroundcolor = 'dark';
backgroundcolor = 'light';

plotvol = 1; 
trigger = 1;
showcolorbar = 1;
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



f = figure('position', [440 37 481 761]);



switch backgroundcolor
	case 'light'
		try
			% load activity_cmap_hot
			load activity_cmap_hot
			cm = cmap;
		catch
			cm = hot(40); 
		end
		try
		cm  = flipud(cbrewer('div', 'RdBu', 40));
		catch
		end

		% !cm = bone(40);
		set(f,'colormap', cm,'color', [1 1 1])

	case 'dark'
		try
			% load activity_cmap_hot
			load activity_cmap_jet
			cm = cmap;
		catch
			cm = jet(40); 
		end
		
		set(f,'colormap', cm,'color', [.2 .2 .2])
	end

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

	

v = reshape(sim.networkHistory.V_soma(:,end), [depth breadth height]);

maxv = 0;
minv = min(min(min(min(sim.networkHistory.V_soma))));

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


for t = frames %1: simtime
	
	v = reshape(sim.networkHistory.V_soma(:,t), [ depth breadth height]);
	% v = permute(v,[3 2 1]);
	
	cla

	spkind =find(v>-20);

	% for visibility of correlations
	thresholdactivity = false;
	if thresholdactivity
		v(v>=-20) = -20;
		v(spkind) = 10;

	end
	
	
	if interpv
		
		Vq = interp3(v, 5);
		% alphamap(1+ 1./(1-exp(linspace(-2,2,100))))
		cdata =Vq;

	elseif ~interpv & plotvol
		
		cdata = v;
		
	elseif scatterit

		v = reshape(sim.networkHistory.V_soma(:,t), 1, depth*breadth*height);
		ac =v;
		ac(1:10) =linspace(minv,0,10);
		
		scatter3(coords(:,1),coords(:,2), coords(:,3), (v+70)*10, ac','filled');
		axis equal
		colorbar
		
	end

		adata = (cdata+70)/80;
		adata2 = adata;
		adata2(find(cdata<-55)) = .3;
		adata2(find(cdata>=-55)) = 0;
		adata2(find(cdata>-53)) = .3;
		cdata= cdata;

		vol3d('cdata',cdata,'Alpha',adata2);

		view(38,34)
		% view(60,10)
		axis off
		axis tight; 
		daspect([1 1 1])
		colorbar
		caxis([-63 -40])
		


	if ismember(t, perturbation_triggers)
		title([num2str(t) 'ms'],'backgroundcolor',[1 0 0], 'color', ones(1,3),'fontsize',30)
	else
		title([num2str(t) 'ms'],'backgroundcolor',[0 0 0], 'color', ones(1,3),'fontsize',30)
	end
	drawnow


	if printframes
		export_fig(['netstate@' num2str(t) ],'-png','-m2')
	end


	if savemovie
		% set(f,'Visible','off')
		currFrame = getframe(f);
		writeVideo(vidObj,currFrame);
	end
end


if savemovie
	close(vidObj)
end

% this is how the olive works: by creating appropriate input, we can generate static phase differences between different muscle groups. 
% These will produce complex spikes in their appropriate 




