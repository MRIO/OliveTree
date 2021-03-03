function replayResults(sim,time_slice, savemovie)

static = 0;

networkHistory  	= sim.networkHistory;
rows 				= sim.rows;
columns				= sim.columns;
lastState 			= sim.lastState;
dt 					= sim.dt;
g_CaL 				= sim.g_CaL;
g_Gap 				= sim.g_Gap;
V_soma_unwrapped 	= sim.networkHistory.V_soma;

simtime				= sim.condition.simtime;
pert_map			= sim.condition.perturbation_map;
perturb_amplitude	= sim.condition.perturb_amplitude;
perturb_onsets      = sim.condition.perturb_onsets;
noise_level			= sim.condition.noise_level;



if ~isempty(time_slice)
	simtime = length(time_slice);
else
	time_slice = [1:simtime];
end


% simtime = [1:length(time_slice)];

V_soma_unwrapped = V_soma_unwrapped(:,time_slice);

V_soma = reshape(V_soma_unwrapped,rows, columns, []);

spks = getSpks(V_soma_unwrapped);

if ishandle(gcf)
	clf
end


i=1; center = [ceil(rows/2) ceil(columns/2)];
for x = 1:rows
	for y =1:columns

		D(i) = (x-center(1))^2+(y-center(2))^2;
		i = i+1;
	end
end
[V O] = sort(D);


if savemovie
	fname = [num2str(rows) 'x' num2str(columns) '_']
	vidObj = VideoWriter('sim.avi','FrameRate', 100);
	open(vidObj);
end

V_soma_unwrapped = V_soma_unwrapped(O,:);

subplot(2,2,[3:4])
	if ~isempty(perturb_onsets)
		if size(perturb_onsets,2)==1; perturb_onsets = perturb_onsets';end
		line([perturb_onsets ; perturb_onsets] , [ones(size(perturb_onsets))*min(V_soma_unwrapped(:)) ;ones(size(perturb_onsets))*max(V_soma_unwrapped(:))],'color','r')
		hold on
		plot(V_soma_unwrapped([1:length(find(pert_map))],:)','g')
		plot(V_soma_unwrapped([1:length(find(~pert_map))],:)')
	else
		plot(V_soma_unwrapped');
	end

	
	xlabel('ms');ylabel('mV'); 
	keyboard
	xlim([time_slice(1) time_slice(end)])
	axis tight

subplot(222)
imagesc(V_soma_unwrapped',[-65 -30]);
 colorbar; 
 hold on;
  ylabel('ms');
  xlabel('neurons');



l = line([1 rows*columns], [0 0],'color', 'r')


if static, simtime=1;end
for t =1:simtime; 
	subplot(221);
	imagesc(V_soma(:,:,t),[-68 -30]);

	title(num2str(t));

	subplot(222)
	set(l,'ydata', [t t])

	
	if savemovie
   % Create an animation.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
	else
		drawnow;
		% keyboard
	end

end               

if savemovie
	% Close the file.
    close(vidObj);
end


%write video
   % Prepare the new file.

    


function spks = getSpks(vsoma)


	[no_neurons time] = size(vsoma)
	for n = 1:no_neurons
		xssings = find(vsoma(n,:)>-30);
		spks.time{n} = xssings(find(diff(xssings)>1));
		spks.count(n) = length(spks.time{n});

	end



	











