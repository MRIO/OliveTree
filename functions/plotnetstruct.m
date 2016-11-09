% plotnetstruct.m

function plotnetstruct(W,X,Y,Z,idx, varargin)
% plotnetstruct(W,X,Y,Z,idx, varargin)

    if isempty(varargin)
        plotconnections = 1;
        plotneurons = 1;
        onlynetstruct = 1;
        plotsingleneuronconnections = 1;
        plotsecondorderneighbors = 0;
        plotconnweights = 0;

    else
        plotconnections = varargin{1}(1);
        plotneurons = varargin{1}(2);
        onlynetstruct = varargin{1}(3);
        plotsingleneuronconnections = 1;
        plotsingleneuronconnections = 1;
        plotsecondorderneighbors = 0;
        plotconnweights = 0;
    end

    noNeurons = size(W,1);
    ncols = length(unique(idx));
    conncolor = [ .6 .6 .6];
    

    try
        % cmap = cbrewer('qual', 'Set1', max(ncols,3));
        % cmap = flipud(cbrewer('seq', 'YlOrRd', max(ncols,3)));
        cmap = cbrewer('seq', 'Purples', 15);
        cmap = cbrewer('div', 'Spectral', length(unique(idx)));
        
    catch
        disp('did not find color brewer')
        cmap = jet(max(ncols,3));
    end

    [ii jj vv] = find(triu(double(W)));
    asym = abs(triu(W)-tril(W)');

if plotconnweights
    normw = vv/max(vv);
    normw = zeros(size(normw));
else
    normw = ones(size(vv));
end

    [vvv iii] = sort(vv);
    
    try
        stdW = quantile(W(:),.15);
    catch
        stdW = 0;
    end

 % index of center neuron
cni = (X==round(max(X)/2)) .* (Y==round(max(Y)/2)) .* (Z==round(max(Z)/2));
cni = find(cni)


    if ~onlynetstruct
        figure
        set(gcf,'defaultaxescolororder', cmap)
        set(gcf,'defaultfigurecolormap', cmap)

         subplot(2,2,[1 3])            


        if plotconnections
            for li = 1:length(ii)

                if vv(li) > stdW & cni~=li

                    line([X(ii(li))  X(jj(li))]', ...
                         [Y(ii(li))  Y(jj(li))]',...
                         [Z(ii(li))  Z(jj(li))]',...
                         'linewidth',(normw(li)+1)/4 ,'color', conncolor * normw(li) )  ;
                
                end
            
            end
        end




        if plotneurons
            hold on
            % scatter3(X,Y,Z,sum(logical(W))*5+eps,sum(logical(W)),'filled') 
            if ~isempty(idx)
                scatter3(X,Y,Z,sum(W>0)*5+eps,idx,'filled') 
            else
                scatter3(X,Y,Z,sum(W>0)*5+eps,sum(W),'filled') 
            end
            colorbar
            title('connections and gap neighborhood')
             
            axis equal
            axis tight
        end
                
        view(-20,10)

            subplot(2,2,2)
            scatter(sum(W~=0) , sum(W))
            title('connections x gap leak')
            xlabel('connections')
            ylabel('total gap leak')

            subplot(2,2,4)
            hist(sum(W))
            title('gap leak to neighbors')
            xlabel('leak')

            drawnow
end


if onlynetstruct
    figure
    set(gca,'colororder', cmap)
    set(gcf,'colormap', cmap)



     if plotconnections
                for li = 1:length(ii)
                    if vv(li) > stdW
                        if  jj(li) ~= cni & ii(li) ~= cni
                        

                        line([X(ii(li))  X(jj(li))]', ...
                             [Y(ii(li))  Y(jj(li))]',...
                             [Z(ii(li))  Z(jj(li))]',...
                             'linewidth',(normw(li)+1)/4 ,'color', conncolor * normw(li) )  ;
                        end
                    end
                
                end
     end


    if plotneurons
        hold on
        % scatter3(X,Y,Z,sum(logical(W))*5+eps,sum(logical(W)),'filled') 
        if ~isempty(idx)
            scatter3(X,Y,Z,sum(W>0)*10+eps,idx,'filled') 
        else
            scatter3(X,Y,Z,sum(W>0)*10+eps,sum(W),'filled') 
        end
        colorbar
        title('connections and gap neighborhood')
         
        axis equal
        axis tight
    end
                    
        view(-120,20)
        drawnow



        
        if plotsingleneuronconnections
            [ii jj vv] = find(triu(double(W)));

            
            

            neighF = @(i) find(W(i,:));

            if plotsecondorderneighbors
                for n = neighF(cni)
                    for nn = neighF(n)

                        if W(n, nn) 

                            line([X(n)  X(nn)]', ...
                                 [Y(n)  Y(nn)]',...
                                 [Z(n)  Z(nn)]',...
                                 'linewidth',3 ,'color', [0 .7 0  ] )  ;
                        end
                    end
                end
            end

            for nn = neighF(cni)
            
                 if W(cni, nn) 

                    line([X(cni)  X(nn)]', ...
                         [Y(cni)  Y(nn)]',...
                         [Z(cni)  Z(nn)]',...
                         'linewidth',5 ,'color', [0 1 0] )  ;
                end
            end
        
        end



end


