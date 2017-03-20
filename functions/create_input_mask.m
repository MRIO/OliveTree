% create_input_mask.m
function mask = create_input_mask(varargin)
	% 
	% requires: cluster toolbox 
	% 
% p = inputParser;
% p.addRequired('netsize')  [x y z]
% p.addRequired('inputmasktype') % mid_layer, random, gauss, all, dist_to_center, 
% p.addParamValue('synapseprobability', .1)
% p.addParamValue('radius', 2)
% p.addParamValue('cell_coordinates', [])
% p.addParamValue('projection_center', [])
% p.addParamValue('plotme', 0)


p = inputParser;
p.addRequired ('netsize')       % a matrix with two columns or a cell array with two cells;
p.addRequired ('inputmasktype') % stdandard deviation criterion for offset threshold
p.addParamValue('synapseprobability', 1)
p.addParamValue('radius', 2)
p.addParamValue('cell_coordinates', [])
p.addParamValue('projection_center', [])
p.addParamValue('plotme', 0)


p.parse(varargin{:});

netsize = p.Results.netsize;
inputmasktype = p.Results.inputmasktype;
synapseprobability = p.Results.synapseprobability;
radius = p.Results.radius;
coords = p.Results.cell_coordinates;
projcenter = p.Results.projection_center;
plotme = p.Results.plotme;


depth = netsize(1);
breadth = netsize(2);
height = netsize(3);

switch inputmasktype

	case 'reconstruction'

        
        coords = netsize;

			center = mean(coords);
			
        	distmat = squareform( pdist([center ; coords], 'euclidean') );

        	mask = (distmat(2:end,1) <=radius);

        	mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;


	case 'mid_layer' % misnomer: one layer uniform excitation
			prob = synapseprobability; 

			mask = zeros(depth, breadth, height);
			mask(1:depth, 1:breadth, ceil(height/2)) = 1;
			mask(1:depth, 1:breadth, floor(height/2)) = 1;
			mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;

			mask = reshape(mask, breadth*height*depth,1); 

	case 'random'
			prob = synapseprobability;
			% A mask that stimluates neurons with 'synapseprobability'
			mask = zeros(depth, breadth, height);
			% HALF RECTIFIED GAUSSIAN
			mask(1:depth, 1:breadth, ceil(height/2)) =  1;

			mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;

			mask = reshape(mask, breadth*height*depth,1); 

	case 'gauss'

			G = gausswin(breadth, gaussalpha(1))*gausswin(depth, gaussalpha(2))';
			g(1,1,:) = gausswin(height, gaussalpha(3));
			mask_gaus = bsxfun(@times,G, g );
			mask_norm = mask_gaus/sum(mask_gaus(:));
			mask_shift= circshift(mask_norm,offset);
			mask = reshape(mask_shift, breadth*height*depth,1); % scale to deliver a total of A pA;

	case 'all'
			mask = ones(depth * breadth*  height,1);

			mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;

	case 'column'
			mask = zeros(depth, breadth, height);
			mask(2:3, 2:3 ,:) = 1;

			mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;

			mask = reshape(mask, breadth*height*depth,1); 

	case 'dist_to_center'

        depth   = [1:netsize(1)];
        breadth = [1:netsize(2)];
        height  = [1:netsize(3)];

        % compute distance matrix
        [X Y Z] = meshgrid(depth,breadth,height);
        X = X(:); Y = Y(:); Z = Z(:);
        coords = [X Y Z];

			center = [mean(X) mean(Y) mean(Z)];
			
        	distmat = squareform( pdist([center ; [X Y Z]], 'euclidean') );

        	mask = (distmat(2:end,1) <=radius);

        	mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;


	case 'dist_to_point'
		
		if isempty(coords)
	        depth   = [1:netsize(1)];
	        breadth = [1:netsize(2)];
	        height  = [1:netsize(3)];

	        [X Y Z] = meshgrid(depth,breadth,height);
	        X = X(:); Y = Y(:); Z = Z(:);

			coords = [X Y Z];
		end

    	distmat = squareform( pdist([projcenter ; coords], 'euclidean') );
  	
    	mask = (distmat(2:end,1) <=radius);

    	mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;


	case 'centerneuron'

			mask = zeros(depth, breadth, height);
			mask(ceil(depth/2), ceil(breadth/2), ceil(height/2)) = 1;

			mask = reshape(mask, breadth*height*depth,1); 

			mask(find(mask)) = rand(size(mask(find(mask)))) <synapseprobability;
			

	otherwise
		disp('FAIL: input mask type not recognized.')
		return


end



% plotme = 1;
% if plotme
% 	scatter3(coords(1), coords(2), coors(3), 50, mask+eps)
% end



    	if plotme
    		scatter3(coords(:,1), coords(:,2), coords(:,3), 100, mask+1,'filled'), axis equal


			cm = [150 147 130 ; 250 35 29]/255;
			colormap(cm)
		end