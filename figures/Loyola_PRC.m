% Figure 7 of Loyola et al. 2023

% What to do

simulate = 1;
plots = 1;



if simulate
    npert = 16;
        
        %% Make neurons

    p1 = [1.1]; 		% Calcium v - conductance range
    p2 = [4.5];  		% g_CaH (needs CaH differences to impact PRC .)
    

    [p{1} p{2} ] = ndgrid(p1,p2);

    nneurons = prod(size(p{1}));
            
        for n = 1:nneurons
            neuron{n} = createDefaultNeurons(1, 'celltypes',  'clones');
            neuron{n}.g_CaL      = p{1}(n);
            neuron{n}.g_CaH      = p{2}(n);
            neuron{n}.gbar_gaba_dend = .1;
            neuron{n}.gbar_gaba_soma = .05;
            neuron{n}.gbar_ampa_dend = .12;
        end


        
    %% template for perturbation
        
    %%
        
    end

    if simulate
    %% perturbation AMPA SUB
        if 1
    	    pert.type	  {1} = 'ampa_dend';
            pert.mask  	  {1} = 1;
    	    pert.duration {1} = 1;
    	    pert.triggers {1} = 0;   
    	    
            for n = 1:nneurons 
                prc_ampa_sub{n} = IO_PRC(neuron{n},pert, npert);
            end
        
            if plots
                plotstuff(prc_ampa_sub)
                title('ampa sub')
            end


            
        end
        
    %% perturbation AMPA SUPRA
        if 0
        	pert.type	  {1} = 'ampa_dend';

            for n = 1:nneurons 
                prc_ampa_supra{n} = IO_PRC(neuron{n},pert, npert);
            end
        end
        
    %% GABA DENDRITE
        if 0
        	pert.duration {1} = 60;
        	pert.type	  {1} = 'gaba_dend';
            pert.trigger  {1} = 0;

            for n = 1:nneurons
                prc_gaba_dend{n} = IO_PRC(neuron{n},pert, npert);
            end

            if plots
                plotstuff(prc_gaba_dend)
                title('gaba dend')
            end
        end
        
    %% GABA SOMA and DEND
        if 1
            pert = [];
        pert.mask  	  {1} = 1;
        pert.amplitude{1} = 1;
        pert.duration {1} = 60;
        pert.type	  {1} = 'gaba_dend';
        pert.triggers {1} = 0;   

        pert.mask  	  {2} = 1;
        pert.amplitude{2} = 1;
        pert.duration {2} = 60;
        pert.type	  {2} = 'gaba_soma';
        pert.triggers {2} = 0;   
    	    
        
            for n = 1:nneurons
                prc_gaba_soma{n} = IO_PRC(neuron{n},pert, npert);
            end
            if plots
                plotstuff(prc_gaba_soma)
                title('gaba dend and soma')
            end
        end

    %% perturbation GABA DENDRITE PSPACE
    if 0
    
    	pert.duration {1} = 2;
    	pert.type	  {1} = 'gaba_dend';
    	
        gabavals = 10;
         
        for n = 1:gabavals
            neuron{1}.gbar_gaba_dend =(n-1)*.4;
            prc_gaba_pspace{n} = IO_PRC(neuron{1},pert, npert);
        end
%%
        if plots
                plotstuff(prc_gaba_soma)
                title('gaba dend and soma')
        end
    end

    
    %% perturbation GABA vs AMPA
    if 1
        pert.mask  	  {1} = 1;
        pert.amplitude{1} = 1;
        pert.duration {1} = 60;
        pert.type	  {1} = 'gaba_dend';
        pert.triggers {1} = 0;   

        pert.mask  	  {2} = 1;
        pert.amplitude{2} = 1;
        pert.duration {2} = 60;
        pert.type	  {2} = 'gaba_soma';
        pert.triggers {2} = 0;   

        pert.mask  	  {3} = 1;
        pert.amplitude{3} = 1;
        pert.duration {3} = 1;
        pert.type	  {3} = 'ampa_dend';
        pert.triggers {3} = 75; 
    
        
    
        neuron{1}.gbar_gaba_dend = .1;
        neuron{1}.gbar_gaba_soma = .1;
        neuron{1}.gbar_ampa_soma = .1;
        for n = 1:nneurons
            prc_comb_sub{n} = IO_PRC(neuron{n},pert, npert);
        end

        if plots
                plotstuff(prc_comb_sub)
                title('comb')
        end

    end

end


%% PLOTING FUNCTION



    


function f = plotstuff(data)
    
    colorlim = [-70 -40];

    try
    cmap = gen_divergent_colormap;
    cmap = flipud(cmap);
    catch
        disp('missing user colormaps')
    end

    

    nneurons = length(data);
    
    for n = 1:nneurons
        f(n) = figure;
        f(n).Color = [1 1 1];            
        subplot(2,1,1)

        pph = data{n}.pertphases; %: [229 237 245 253 261 269 277 285 292 300 308 316 324 332 340 348]
        pph = -[pph(1) pph];

        stst = data{1}.steady_state.networkHistory.V_soma;
        VS = data{n}.VS;
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
        % plot(whole_traces_cropped(1,:),'k-')
        % plot(whole_traces_cropped(9,:),'k--')


        X = whole_traces_shifted_cropped;
        sz = size(X);
        ampx = (max(X(:)) - min(X(:)));
        traces = X/ampx;
        plot([traces + repmat([sz(1):-1:1]',1, sz(2))]', 'k')
        axis tight
        axis off

        title({['neuron' num2str(n)] ; ['CaL: ' num2str(data{n}.neuron.g_CaL) ]; ['CaH: ' num2str(data{n}.neuron.g_CaH)]})
        
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
    