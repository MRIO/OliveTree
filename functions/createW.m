function out = createW(varargin)
% plotnetstruct(W.W, W.coords(:,1), W.coords(:,2), W.coords(:,3), [0 1 0])
% out = createW('type', netsize, radius, scaling, randomize, plotthis, maxiter, meanconn, somatapositions, symmetrize, clusterize)
% 
%  creates connectivity according to rules
% 
% 
% [=================================================================]
%    Creates a connectivity matrix
% [=================================================================]
% 
% Required('connections')  % type of connectivity (string)
% Required('netsize')      % size of network (in 3D: [X Y Z])
% 
% [=================================================================]
%  Parameter Value Pairs
% [=================================================================]
% 
% Optional('radius',1)     % radius of connectivity
% Optional('scaling',1)    % the unitary weight
% Optional('randomize',0)  % post processing: randomization of weights of 10%
% Optional('maxiter', 3)   % this property only applies if connection type is 'rndwalk' 
% Optional('meanconn', 0)  % mean number of connections per cell (enforced)
% Optional('somatapositions', []) % 
% Optional('symmetrize',1)        % symmetrize connections (also used in the case of randomized)
% Optional('clusterize', [0 0 0 1]); % cluster is applied as multiplicative mask
%       four values, clusterize(1) = 0/1 (on or off)
%                    clusterize(2) = groupsize
%                    clusterize(3) = clusterconn
%                    clusterize(4) = otherconn
% 
% Optional('normleak', 0);        % normalization of the weights by total leak (per weight sum of W = scaling)
% Optional('seed', 0);            % random seed for
% Optional('plotthis',1)          % plot the connectivity of the network
% 
% [=================================================================]
%  Output
% [=================================================================]
% 
% out.W = W; %   the actual weight matrix 
% 
% out.params.connections     = connections;
% out.params.netsize         = netsize;
% out.params.radius          = radius;
% out.params.scaling         = scaling;
% out.params.randomize       = randomize;
% out.params.plotthis        = plotthis;
% out.params.maxiter         = maxiter;
% out.params.meanconn        = meanconn;
% out.params.somatapositions = somatapositions;
% out.params.symmetrize      = symmetrize;
% out.params.clusterize      = clusterize;
% out.coords = [X Y Z];
% out.stats.clusters = idx;
% out.params.clusterize = clusterize;
debugging = 1;


    p = inputParser;
    p.addRequired('connections')  % type of connectivity
    p.addRequired('netsize')      % size of network (in 3D: [X Y Z])
    
    p.addOptional('radius',1)     % radius of connectivity
    p.addOptional('scaling',1)    % the unitary weight
    p.addOptional('randomize',0)  % post processing: randomization of weights of 10%
    p.addOptional('plotthis',1)
    p.addOptional('maxiter', 3)   % this property only applies if connection type is 'rndwalk' 
    p.addOptional('meanconn', 0)  % mean number of connections per cell (enforced)
    p.addOptional('somatapositions', []) % 
    p.addOptional('symmetrize',1)        % symmetrize connections (also used in the case of randomized)
    p.addOptional('clusterize', [0 0 0 1]); % four values, [0/1 groupsize clusterconn otherconn]
    p.addOptional('normleak', 0);           % [0/1 groupsize clusterconn otherconn]
    p.addOptional('seed', 0);


    
    p.parse(varargin{:});

    connections     = p.Results.connections;
    netsize         = p.Results.netsize;
    radius          = p.Results.radius;
    scaling         = p.Results.scaling;
    randomize       = p.Results.randomize;
    plotthis        = p.Results.plotthis;
    maxiter         = p.Results.maxiter;
    meanconn        = p.Results.meanconn;
    somatapositions = p.Results.somatapositions;
    symmetrize      = p.Results.symmetrize;
    clusterize      = p.Results.clusterize;
    normleak        = p.Results.normleak;
    seed            = p.Results.seed;
    

noNeurons = prod(netsize);
rng(seed,'twister')

% [=================================================================]
%  coords from a grid
% [=================================================================]
        if not(isempty(somatapositions))
            netsize = [1 length(somatapositions) 1];
        end

        depth   = [1:netsize(1)];
        breadth = [1:netsize(2)];
        height  = [1:netsize(3)];

        % # compute adjacency matrix
        [X Y Z] = meshgrid(depth,breadth,height);
        X = X(:); Y = Y(:); Z = Z(:);


% [=================================================================]
%  connectivity rules
% [=================================================================]

switch connections

    case '3d_reconstruction'
        
        distFunc = 'euclidean';

        % overwrite X,Y,Z with function argument
        X = somatapositions(:,1);
        Y = somatapositions(:,2);
        Z = somatapositions(:,3);

        noNeurons = length(X);

        W = squareform( pdist([X Y Z], distFunc) <= radius );

    case 'all to all'

            W = ones(noNeurons)-eye(noNeurons);
            plotthis = 0;
            
    case 'random'

            W = ones(noNeurons)-eye(noNeurons);
            randomize = 1;

    case '3d_chebychev'
        distFunc = 'chebychev';


        W = squareform( pdist([X Y Z], distFunc ) <= radius );

    case '3d'

        W = squareform( pdist([X Y Z], 'euclidean') <= radius );         

    case '3d_euclidean' % redundant for backward compatibility

        W = squareform( pdist([X Y Z], 'euclidean') <= radius );


    case 'one_cluster'

        clusterradius = 3;
        ProbCluster = clusterize(3);
        ProbOriginal = clusterize(4);

        depth   = [1:netsize(1)];
        breadth = [1:netsize(2)];
        height  = [1:netsize(3)];

        % compute distance matrix
        [X Y Z] = meshgrid(depth,breadth,height);
        X = X(:); Y = Y(:); Z = Z(:);
        coords = [X Y Z];

        center = [mean(X) mean(Y) mean(Z)];
        
        distmat = squareform( pdist([center ; [X Y Z]], 'euclidean') );

        mask = (distmat(2:end,1) <=clusterradius);

        idx(find(mask))=1;
        idx(find(not(mask)))=2;

        W1 = zeros(noNeurons);
        W1(find(mask==0),:)   = 1;
        W1 = and(W1,W1');

        W2 = zeros(noNeurons);
        W2(find(mask==1),:)   = 1;
        W2 = and(W2,W2');

        W3 = squareform( pdist([X Y Z], 'euclidean') <= radius );

        W = ((W1.*W3 + W2.*W3 ).* rand(noNeurons) ) > (1 - ProbCluster);


    case '3d_euclidean_rndwalk'
        % deforms a perfect latice according to connectivity criteria
     
        % maximum step for randomizing positions
        rndsteplength =   .1*radius;
        
        distFunc = 'euclidean';

        W = squareform( pdist([X Y Z], distFunc) <= radius );

        cprob = quantile(sum(W),[.25 .75]);

        % global CG
        CG = [mean(X) mean(Y) mean(Z)];

        % local CG
        CGx = arrayfun(@(l) mean(X(find(W(:,l)))) , 1:noNeurons)';
        CGy = arrayfun(@(l) mean(Y(find(W(:,l)))) , 1:noNeurons)';
        CGz = arrayfun(@(l) mean(Z(find(W(:,l)))) , 1:noNeurons)';


        FEW = find(sum(W)<=cprob(1));
        MANY = find(sum(W)>=cprob(2));

        it = 1; sts = rndsteplength; maxiter = 20;
        while it < maxiter

             % approach or avoid local CG
            CGx = arrayfun(@(l) mean(X(find(W(:,l)))) , 1:noNeurons)';
            CGy = arrayfun(@(l) mean(Y(find(W(:,l)))) , 1:noNeurons)';
            CGz = arrayfun(@(l) mean(Z(find(W(:,l)))) , 1:noNeurons)';


            dolocal = 1;
            if dolocal

                X(FEW)  =   sts*rand*( CGx(FEW) - X(FEW)  ) + X(FEW);
                Y(FEW)  =   sts*rand*( CGy(FEW) - Y(FEW)  ) + Y(FEW);
                Z(FEW)  =   sts*rand*( CGz(FEW) - Z(FEW)  ) + Z(FEW);
                X(MANY) = - sts*rand*( CGx(MANY) - X(MANY)) + X(MANY);
                Y(MANY) = - sts*rand*( CGy(MANY) - Y(MANY)) + Y(MANY);
                Z(MANY) = - sts*rand*( CGz(MANY) - Z(MANY)) + Z(MANY);


            else % global CG -- in the limit transforms data to spherical surface
            
                X(FEW)  =   sts*rand*(CG(1) - X(FEW)) + X(FEW);
                Y(FEW)  =   sts*rand*(CG(2) - Y(FEW)) + Y(FEW);
                Z(FEW)  =   sts*rand*(CG(3) - Z(FEW)) + Z(FEW);
                X(MANY) = - sts*rand*(CG(1) - X(MANY)) + X(MANY);
                Y(MANY) = - sts*rand*(CG(2) - Y(MANY)) + Y(MANY);
                Z(MANY) = - sts*rand*(CG(3) - Z(MANY)) + Z(MANY);
            end

            Wnew = squareform( pdist([X Y Z], distFunc) <= radius ) ;

            FEW = find(sum(Wnew)<=cprob(1));
            MANY = find(sum(Wnew)>=cprob(2));

            % pause(1)
            % clf
            % plotconnectionmatrix(Wnew,X,Y,Z)
            
            it = it + 1;

        end

        W = Wnew;

    otherwise
        disp('case not recognized')
        return

end

Wsteps{1} = W; 
length(find(W.*eye(size(W))))


% [=================================================================]
%  postprocessing
% [=================================================================]

% remove self connections
W = W.*not(eye(size(W)));


% #TODO: recheck if self connections are handled elegantly

if clusterize(1)

    groupsize = clusterize(2);

    k = round(noNeurons/groupsize);
    [idx] = kmeans([X Y Z], k);

    ProbCluster = clusterize(3);
    ProbOriginal = clusterize(4);

    cW = zeros(noNeurons);
    for i = 1:noNeurons
        for j = 1:noNeurons
            if idx(i)==idx(j)
                cW(i,j) = 1;
            end
        end
    end

    % W = cW;
    W = cW.* ( (rand(noNeurons)+(eye(noNeurons))) <= ProbCluster) + (W.*~cW).* ( rand(noNeurons) <= ProbOriginal) ;

    out.stats.clusters = idx;

end

Wsteps{2} = W;

% prune to mean number of connections 
if meanconn

    R = rand(size(W));
    R = (triu(R) + triu(R)')+(eye(size(W)));
    W = bsxfun(@rdivide, W, sum(W>0,2)/meanconn)>= R;

end
Wsteps{3} = W;

% ensure_meanconn = 1;

% if ensure_meanconn
%     D = squareform( pdist([X Y Z], distFunc));
%     W = ensure_connectivity(W, D)
% end



% symmetric or asymmetric gap junctions
if symmetrize
    W = triu(W) + triu(W)';
end
Wsteps{4} = W;


if normleak
    deno = sum(W');
        deno(deno==0)=1;
        W = bsxfun(@rdivide, W', deno);

end
Wsteps{5} = W;


% add a random multiplier for confusion
if randomize 
    R = rand(size(W))*randomize*scaling*.1; % jitter weights  
    R(find(eye(size(R)))) = 0; % no self connections, naturally.
    W = sparse(W .* R + W);
    if symmetrize
        W = triu(W) + triu(W)';
    end

end
Wsteps{6} = W;

% scale single gaps such that average leak per cell is according to scale
if scaling
    W = W / (sum(sum(W))/sum(sum(W~=0)))*scaling;
end
Wsteps{7} = W;

if not(exist('idx')) 
        idx = sum(W);
end

% [=================================================================]
%  finalize
% [=================================================================]

W = full(W);

% get rid of nan's from potential division by zero
W(isnan(W)) = 0;



% [=================================================================]
%  outputs
% [=================================================================]

% stats = connectivity_statistics(W);

out.stats.stdW = std(W(W~=0));
out.stats.connections = sum(W>0);
out.stats.meanweight  = mean(W(find(W>0)));

if exist('clustering_coef_wd')
    out.stats.clustercoeff.wd = clustering_coef_wd(W);
    out.stats.clustercoeff.bu = clustering_coef_bu(W);
else
    out.stats.clustercoeff = 'WARNING: did not find connectivity toolbox';
end

out.W = W;

out.coords = [X Y Z];
out.stats.clusters = idx;
out.params.clusterize = clusterize;

out.params.connections     = connections;
out.params.netsize         = netsize;
out.params.radius          = radius;
out.params.scaling         = scaling;
out.params.randomize       = randomize;
out.params.maxiter         = maxiter;
out.params.meanconn        = meanconn;
out.params.somatapositions = somatapositions;
out.params.symmetrize      = symmetrize;
out.params.clusterize      = clusterize;

% out.stats.betweeness_wei = betweenness_wei(W);
% out.stats.dbscanclusters


if debugging
out.collection = Wsteps;
end


if plotthis==1
    warning('off') % colorbrewer throws an annoyance
           % plotconnectionmatrix(W,X,Y,Z,idx,out)
 
           % plotnetstruct(W,X,Y,Z,idx,out)
           plotnetstruct(W,X,Y,Z,idx)
    warning('on')
elseif plotthis==2
    figure
    imagesc(W)
    title('all to all connections')
    colorbar
end



% function plotconnectionmatrix(W,X,Y,Z,idx,out)


%     noNeurons = size(W,1);
%     plotconnections = 1;
%     plotneurons = 1;
%     onlynetstruct = 1;

%     ncols = length(unique(idx));

%     if strcmp(out.params.connections, 'one_cluster') | out.params.clusterize(1)==1
%         cmaptype = 'qual';
%     else
%         cmaptype = 'div';
%     end


%     try 
%         switch cmaptype
%             case 'qual'
%             cmap = cbrewer('qual', 'Set1', max(ncols,3));
%             case 'div'
%             cmap = cbrewer('div', 'RdBu', ncols);
%         end
%     catch
%         switch cmaptype
%             case 'qual'
%             cmap = jet(ncols);
%             case 'div'
%             cmap = colorcube(ncols);
%         end
%     end
    
%     set(0,'defaultaxescolororder', cmap)
%     set(0,'defaultfigurecolormap', cmap)


%     [ii jj vv] = find(triu(double(W)));
%     asym = abs(triu(W)-tril(W)');
%     normw = vv/max(vv);
%     [vvv iii] = sort(vv);
%     try
%         stdW = quantile(W(:),.15);
%     catch
%         stdW = 0;
%     end


%     if ~onlynetstruct
%         figure

%          subplot(2,2,[1 3])            


%         if plotconnections
%             for li = 1:length(ii)



%                 clusterized = 1;
%                 if clusterized
%                     if idx(ii) == idx(jj)
%                     line([X(ii(li))  X(jj(li))]', ...
%                          [Y(ii(li))  Y(jj(li))]',...
%                          [Z(ii(li))  Z(jj(li))]',...
%                          'linewidth',normw(li) ,'color', [ 1 .2 .2] * normw(li) )  ;
%                     else
%                     line([X(ii(li))  X(jj(li))]', ...
%                          [Y(ii(li))  Y(jj(li))]',...
%                          [Z(ii(li))  Z(jj(li))]',...
%                          'linewidth',normw(li) ,'color', [ .2 .2 .2] )  ;
%                     end
%                 else
%                     if vv(li) > stdW

%                     line([X(ii(li))  X(jj(li))]', ...
%                          [Y(ii(li))  Y(jj(li))]',...
%                          [Z(ii(li))  Z(jj(li))]',...
%                          'linewidth',normw(li) ,'color', [ 1 .2 .2] * normw(li) )  ;
                
%                     end
%                 end
            
%             end
%         end


%         if plotneurons
%             hold on
            
%             if ~isempty(idx)
%                 scatter3(X,Y,Z,sum(W>0)*5+eps,idx,'filled') 
%             else
%                 scatter3(X,Y,Z,sum(W>0)*5+eps,sum(W),'filled') 
%             end
%             colorbar
%             title('connections and gap neighborhood')
             
%             axis equal
%             axis tight
%         end
                
%         view(-20,10)

%             subplot(2,2,2)
%             scatter(sum(W~=0) , sum(W))
%             title('connections x gap leak')
%             xlabel('connections')
%             ylabel('total gap leak')

%             subplot(2,2,4)
%             hist(sum(W))
%             title('gap leak to neighbors')
%             xlabel('leak')

%             drawnow
% end


% if onlynetstruct
%     figure

%      if plotconnections
%                 for li = 1:length(ii)
%                     if vv(li) > stdW

%                         line([X(ii(li))  X(jj(li))]', ...
%                              [Y(ii(li))  Y(jj(li))]',...
%                              [Z(ii(li))  Z(jj(li))]',...
%                              'linewidth',normw(li)*2 ,'color', [ 1 .2 .2] * normw(li) )  ;
%                     end
                
%                 end
%      end


%     if plotneurons
%         hold on
%         % scatter3(X,Y,Z,sum(logical(W))*5+eps,sum(logical(W)),'filled') 
%         if ~isempty(idx)
%             scatter3(X,Y,Z,sum(W>0)*7+eps,idx,'filled') 
%         else
%             scatter3(X,Y,Z,sum(W>0)*7+eps,sum(W),'filled') 
%         end
%         colorbar
%         title('connections and gap neighborhood')
         
%         axis equal
%         axis tight
%     end
                    
%         view(-120,20)
%         drawnow


% end



