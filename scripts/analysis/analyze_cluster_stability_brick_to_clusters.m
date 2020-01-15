% analyze_cluster_stability.m

max_no_clusters = 50;
plot_all = 0;
Nshuffle = 50;


% Load the data
load('/Users/M/Synced/Titan/Bench4/brick_to_clusters_06-May-2019.mat')

for nsims = 1:5

    z = sims{nsims}.networkHistory.V_soma;

    tt = 2001:10000;

    % Remove the beginning for better processing
    z1 = double(z(:, tt)');
    % Data without a dc component
    z2 = bsxfun(@minus, z1, mean(z1,1));

     
    % Here we first remove the non-oscllating cells.
    z3 = hilbert(z2);
    logamp = log(mean(abs(z3)));
    if plot_all
        hist(logamp, 500)
        box off
    end

    % Two peaks at the end are from non-oscillating cells. We remove them from our data...

     
    igood = find(logamp>=-5.4);
    ibad = find(logamp<-5.4);

    z4 = z2(:,igood);
    if plot_all
        % non-oscillating cells
        plot(z2(:,ibad))

    end


    % Now we compute the normalized data and (partial) correlation matrix.

     
    z4 = bsxfun(@times, z4, 1./std(z4)); % Normalize the amplitude as well
    z4m = mean(z4, 2);

    % zn = (z4'*z4)/size(z4, 1); % either correlation
    zn = partialcorr(z4, z4m); % partial correlation w.r.t the mean oscillation
    zn(zn<0) = 0; % Ignore negative correlation...
     

    if plot_all
        figure
        imagesc(zn)
        xlabel('Cell')
        ylabel('Cell')
        c = colorbar();
        c.Label.String = 'Correlation';
    end


    % Cluster stability is a non-parametric method to determine the best number of clusters. First, we make  data sets from randomly drawn samples (a fraction of 80% below), obtain cluster labels with a given number of clusters for each resampled data set, and compare how well they match with each other.

     
    
    [zn2, idxz] = make_rand_samples2(zn, Nshuffle, 0.8);
     
    v = spctcldata(zn);
    vs = cellfun(@spctcldata, zn2, 'UniformOutput', false);

    % Now, we evaluate the distance of clustering results between each pair of random sample, with changing the number of clusters from 2 to 120.

     
    ngrid = max_no_clusters - 1;
    met = zeros(ngrid, Nshuffle*(Nshuffle-1)/2);
    rmet = zeros(ngrid, Nshuffle*(Nshuffle-1)/2);
    fprintf('i, Ncl\n');
    for i=1:ngrid
        Ncl = i+1;
        
        % Perform clustering on the random samples
        clidx = arrayfun(@(k) spctcl(vs{k}, Ncl, 1), 1:Nshuffle, 'UniformOutput', false);
        
        % Evaluate distance between each pair
        [pd, rpd] = cij2(clidx, idxz);
        met(i,:) = pd;
        rmet(i,:) = rpd;
        
        fprintf('%d, %d\n', i, Ncl);
    end

    % met and rmat are the cosine distances 
     % for the actual and randomized (control) data.
      % We correct for the data and label size effect by using rmet and compute the stability . 
       % It turned out that maximizes when .

    figure 
        tcl = 2:max_no_clusters;
        % nothing
        mmet = 1-mean(met'./rmet');
        sdmet = std(met'./rmet');
        plot(tcl, mmet+sdmet,':k',...
             tcl, mmet-sdmet,':k',...
             tcl, mmet,'k')
        box off
        axis tight
        ylim([0.45 1])
        xlabel('Clusters')
        ylabel('Stability')
        title(num2str(nsims))

        output.tcl = tcl;
        output.mmet{nsims} = mmet;
        output.sdmet{nsims} = sdmet;

    if plot_all


        % The correlation matrix looks nice when ncl = ?.

        % figure
        % ncl = 12;
        % idx = spctcl(v, ncl, 20); %Recluster the data
        % [~, sidx] = sort(idx);
        % imagesc(zn(sidx, sidx))
        % xlabel('Cell (sorted)')
        % ylabel('Cell (sorted)')


        % Examining the membrane potentials for each cluster shows that clustering is mostly driven by locking at different phases.

         
        %plot inline -s 800,600


        figure
        %plot inline -s 1000,800
        subplot(121)
        imagesc(z4(:, sidx)')
        xlabel('Time (ms)')
        ylabel('Cell (sorted)')
        caxis([-2 2])
        title('Membrane potential (normalized)')
        subplot(122)
        imagesc(zn(sidx, sidx))
        xlabel('Cell (sorted)')
        ylabel('Cell (sorted)')
        title('Correlation matrix' )


        figure 
        hold on
        plot(tcl, cell2mat(output.mmet'))
        legend

        box off
        axis tight
        ylim([0.45 1])
        xlabel('Clusters')
        ylabel('Stability')
        title(num2str(nsims))

        output.tcl = tcl;
        output.mmet{nsims} = mmet;
        output.sdmet{nsims} = sdmet;





    end
end


bla = cell2mat(output.mmet')
 plot(tcl, bla([1 3 5],:))

        box off
        axis tight
        ylim([0.45 1])
        xlabel('Clusters')
        ylabel('Stability')
legend({'brick', 'bridges', 'clusters'})