% Figure 7 of Loyola et al. 2023

% What to do
prepare = 0;
simulate = 0;
plots = 1;


try
    cmap = gen_divergent_colormap;
    cmap = flipud(cmap);
catch
    disp('missing user colormaps')
end


if prepare
    npert = 16;
        
        %% Make neurons

    p1 = [1.1]; 		% Calcium v - conductance range
    p2 = [4.5];  		% g_CaH (does not impact PRC without CaH.)
    

    [p{1} p{2} ] = ndgrid(p1,p2);

    nneurons = prod(size(p{1}));
            
        for n = 1:nneurons
            neuron{n} = createDefaultNeurons(1, 'celltypes',  'clones');
            neuron{n}.g_CaL     = p{1}(n);
            neuron{n}.g_CaH      = p{2}(n);
        end


        
    %% template for perturbation
        pert.mask  	  {1} = 1;
    	pert.amplitude{1} = 1;
    	pert.duration {1} = 1;
    	pert.triggers{1} = [100];   
    %%
        
    end

    if simulate
    %% perturbation AMPA SUB
    	pert.type	  {1} = 'ampa_dend';
    	
        for n = 1:nneurons 
            neuron{n}.gbar_ampa_dend = .1;
            prc_ampa_sub{n} = IO_PRC(neuron{n},pert, npert);
        end
        
    %% perturbation AMPA SUPRA
        if 0
        	pert.type	  {1} = 'ampa_dend';

            for n = 1:nneurons 
                neuron{n}.gbar_ampa_dend = .5;
                prc_ampa_supra{n} = IO_PRC(neuron{n},pert, npert);
            end
        end
        
    %% GABA SUPRATHRESHOLD REBOUND
        if 0
        	pert.duration {1} = 4;
        	pert.type	  {1} = 'gaba_dend';

            for n = 1:nneurons
                neuron{n}.gbar_gaba_dend = 2;
                prc_gaba_supra{n} = IO_PRC(neuron{n},pert, npert);
            end
        end
        
    %% GABA SUBTHRESHOLD
        pert.duration {1} = 1;
    	pert.type	  {1} = 'gaba_dend';
    	
        for n = 1:nneurons
            neuron{n}.gbar_gaba_dend = 1;
            prc_gaba_sub{n} = IO_PRC(neuron{n},pert, npert);
        end

    if 0
    %% perturbation GABA DENDRITE PSPACE
    	pert.duration {1} = 2;
    	pert.type	  {1} = 'gaba_dend';
    	
        gabavals = 10;
         
        for n = 1:gabavals
            neuron{1}.gbar_gaba_dend =(n-1)*.4;
            prc_gaba_pspace{n} = IO_PRC(neuron{1},pert, npert);
        end
    end


    %% perturbation GABA vs AMPA
    pert.mask  	  {1} = 1;
    pert.amplitude{1} = 1;
    pert.duration {1} = 1;
    pert.type	  {1} = 'ampa_dend';
    pert.triggers{1} = [60];  

    pert.mask  	  {2} = 1;
    pert.amplitude{2} = 1;
    pert.duration {2} = 3;
    pert.type	  {2} = 'gaba_dend';
    pert.triggers{2} = [0];   


    for n = 1:nneurons
        neuron{n}.gbar_gaba_dend = 1;
        prc_comb_sub{n} = IO_PRC(neuron{n},pert, npert);
    end

end

%% PLOTS

%%
if plots



nneurons = 1;
    % cmap = cbrewer('div', 'RdYlBu', 64);
    colorlim = [-70 -40];


    for n = 1:nneurons
        f1 = figure;
        f1.Color = [1 1 1];            
        subplot(2,1,1)

        pph = prc_gaba_sub{n}.pertphases; %: [229 237 245 253 261 269 277 285 292 300 308 316 324 332 340 348]
        pph = -[pph(1) pph];

        stst = prc_gaba_sub{1}.steady_state.networkHistory.V_soma;
        VS = prc_gaba_sub{n}.VS;
        


        stst = repmat(stst, 17,1);
        whole_traces = [stst VS];
        whole_traces_cropped = whole_traces(:,300:end-300);
        whole_traces_shifted = cell2mat(arrayfun(@(row) circshift(whole_traces(row, :), [0 pph(row)-pph(1) ]) , [1:17], 'uniformoutput', 0)');
        whole_traces_shifted_cropped = whole_traces_shifted(:,300:end-300);
    

        hilb_vs = hilbert_of_membranepotential(whole_traces_shifted_cropped(2:end,:))


        sp(1) = subplot(5,1,1)
        % plot(whole_traces_shifted_cropped(1,:),'k-')
        % hold on
        % plot(whole_traces_shifted_cropped(7,:),'k--')
        plot(whole_traces_cropped(1,:),'k-')
        hold on
        plot(whole_traces_cropped(9,:),'k--')
        axis tight
        axis off
        ylim([-80,0])

        title({['gaba sub']; ['neuron' num2str(n)] ; ['CaL: ' num2str(p{1}(n)) ]; ['CaH: ' num2str(p{2}(n))]})
        
        sp(2) = subplot(5,1,2)
        
        imagesc(whole_traces_shifted_cropped(2:end,:),colorlim)
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])

        colormap(sp(2), parula)
        axis tight
        axis off
        

        sp(3) = subplot(5,1,3)

        plot(whole_traces_shifted_cropped(2:end-1,:)','k')
        hold on
        plot(mean(whole_traces_shifted_cropped(2:end-1,:)), 'r', 'linewidth',2)
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end-1,:))
        axis tight
        axis off

        
        % colormap(sp(3), parula)

        


        sp(4) = subplot(5,1,4)

        imagesc(hilb_vs.hilbert)
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])
        colormap(sp(4), cmap)
        axis off

        

        sp(5) = subplot(5,1,5)
        
        plot_mean_and_std(hilb_vs.hilbert(1:end-1,:))
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end,:))

        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])


    end

    %%


    
    for n = 1:nneurons
        f1 = figure;
        f1.Color = [1 1 1];            

        pph = prc_ampa_sub{n}.pertphases; %: [229 237 245 253 261 269 277 285 292 300 308 316 324 332 340 348]
        pph = -[pph(1) pph];

        stst = prc_ampa_sub{1}.steady_state.networkHistory.V_soma;
        VS = prc_ampa_sub{n}.VS;
        stst = repmat(stst, 17,1);
        whole_traces = [stst VS];
        whole_traces_cropped = whole_traces(:,300:end-300);
        whole_traces_shifted = cell2mat(arrayfun(@(row) circshift(whole_traces(row, :), [0 pph(row)-pph(1) ]) , [1:17], 'uniformoutput', 0)');
        whole_traces_shifted_cropped = whole_traces_shifted(:,300:end-300);
    

        hilb_vs = hilbert_of_membranepotential(whole_traces_shifted_cropped(2:end,:))


        sp(1) = subplot(5,1,1)
        % plot(whole_traces_shifted_cropped(1,:),'k-')
        % hold on
        % plot(whole_traces_shifted_cropped(7,:),'k--')
        plot(whole_traces_cropped(1,:),'k-')
        hold on
        plot(whole_traces_cropped(9,:),'k--')
        axis tight
        axis off
        ylim([-80,0])

        title({['ampa sub']; ['neuron' num2str(n)] ; ['CaL: ' num2str(p{1}(n)) ]; ['CaH: ' num2str(p{2}(n))]})
        
        sp(2) = subplot(5,1,2)
        
        imagesc(whole_traces_shifted_cropped(2:end,:),colorlim)
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])

        colormap(sp(2), parula)
        axis tight
        axis off
        

        sp(3) = subplot(5,1,3)
  plot(whole_traces_shifted_cropped(2:end-1,:)','k')
        hold on
        plot(mean(whole_traces_shifted_cropped(2:end-1,:)), 'r', 'linewidth',2)
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end-1,:))
        axis tight
        axis off

        


        sp(4) = subplot(5,1,4)

        imagesc(hilb_vs.hilbert)
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])
        colormap(sp(4), cmap)
        axis off

        

        sp(5) = subplot(5,1,5)
        
        plot_mean_and_std(hilb_vs.hilbert(1:end-1,:))
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end,:))

        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])


    




    end








    
    for n = 1:nneurons
        f1 = figure;
        f1.Color = [1 1 1];            
        subplot(2,1,1)

        pph = prc_comb_sub{n}.pertphases; %: [229 237 245 253 261 269 277 285 292 300 308 316 324 332 340 348]
        pph = -[pph(1) pph];

        stst = prc_comb_sub{1}.steady_state.networkHistory.V_soma;
        VS = prc_comb_sub{n}.VS;
        stst = repmat(stst, 17,1);
        

        whole_traces = [stst VS];
        whole_traces_cropped = whole_traces(:,300:end-300);
        whole_traces_shifted = cell2mat(arrayfun(@(row) circshift(whole_traces(row, :), [0 pph(row)-pph(1) ]) , [1:17], 'uniformoutput', 0)');
        whole_traces_shifted_cropped = whole_traces_shifted(:,300:end-300);
    

        hilb_vs = hilbert_of_membranepotential(whole_traces_shifted_cropped)


        sp(1) = subplot(5,1,1)
        % plot(whole_traces_shifted_cropped(1,:),'k-')
        % hold on
        % plot(whole_traces_shifted_cropped(7,:),'k--')
        plot(whole_traces_cropped(1,:),'k-')
        hold on
        plot(whole_traces_cropped(9,:),'k--')
        
        axis tight
        axis off
        ylim([-80,0])

        title({['combined']; ['neuron' num2str(n)] ; ['CaL: ' num2str(p{1}(n)) ]; ['CaH: ' num2str(p{2}(n))]})
        
        sp(2) = subplot(5,1,2)
        
        imagesc(whole_traces_shifted_cropped(2:end,:),colorlim)
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])

        colormap(sp(2), parula)
        axis tight
        axis off
        

        sp(3) = subplot(5,1,3)
        
        plot(whole_traces_shifted_cropped(2:end-1,:)','k')
        hold on
        plot(mean(whole_traces_shifted_cropped(2:end-1,:)), 'r', 'linewidth',2)
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end-1,:))
        axis tight
        axis off

        
        


        sp(4) = subplot(5,1,4)

        imagesc(hilb_vs.hilbert(2:end,:))
        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])
        colormap(sp(4), cmap)
        axis off

        

        sp(5) = subplot(5,1,5)
        
        plot_mean_and_std(hilb_vs.hilbert(2:end-1,:))
        % plot_mean_and_std(whole_traces_shifted_cropped(2:end,:))

        % colororder = flipud(cbrewer('seq', 'Blues', 18));
        % colororder(1:2,:)=[];
        % ab.ColorOrder = colororder;
        xlim([0 2000])


    
    end
    
    


end