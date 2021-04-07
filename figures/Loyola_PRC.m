
%% Make neurons

npert = 16;

		p1 = [.8 1.1]; 		% CalciumL - conductance range
		p2 = [.12 1.2];  		% g_h
		
		[p{1} p{2} ] = ndgrid(p1,p2);

nneurons = prod(size(p{1}));
        
    for n = 1:nneurons
        neuron{n} = createDefaultNeurons(1, 'celltypes',  'clones');
        neuron{n}.g_CaL    = p{1}(n);
        neuron{n}.g_h 	    = p{2}(n);
    end
       
 %% make perturbation
    
    pert.mask  	  {1} = 1;
	pert.amplitude{1} = 1;
	pert.duration {1} = 1;
	pert.triggers{1} = [100];   

%% perturbation AMPA SUB
	pert.type	  {1} = 'ampa_dend';
	
    for n = 1:nneurons 
        neuron{n}.gbar_ampa_dend = .1;
        prc_ampa_sub{n} = IO_PRC(neuron{n},pert, npert);
    end
    
%% perturbation AMPA SUPRA
	pert.type	  {1} = 'ampa_dend';

    for n = 1:nneurons 
        neuron{n}.gbar_ampa_dend = .5;
        prc_ampa_supra{n} = IO_PRC(neuron{n},pert, npert);
    end

    
%% GABA SUPRATHRESHOLD
	pert.duration {1} = 4;
	pert.type	  {1} = 'gaba_dend';

    for n = 1:nneurons
        neuron{n}.gbar_gaba_dend = 2;
        prc_gaba_supra{n} = IO_PRC(neuron{n},pert, npert);
    end
    
%% GABA SUBTHRESHOLD
	pert.type	  {1} = 'gaba_dend';
	
    for n = 1:nneurons
        neuron{n}.gbar_gaba_dend = 1;
        prc_gaba_sub{n} = IO_PRC(neuron{n},pert, npert);
    end

 
%% perturbation GABA DENDRITE PSPACE
	pert.duration {1} = 4;
	pert.type	  {1} = 'gaba_dend';
	
    gabavals = 10;
    pertphases = 16;
     
    for n = 1:gabavals
        neuron{1}.gbar_gaba_dend =(n-1)*.4;
        prc_gaba_pspace_d4{n} = IO_PRC(neuron{1},pert, pertphases);
    end
    
 
 %%
    
    f1 = figure;
    f1.Color = [1 1 1];
    for n = 1:nneurons
      
        ax(n) = subplot(2,2,n)
        ttt = [1:500]- prc_gaba_sub{2}.peaktimes(1);
        waterfall(prc_gaba_sub{n}.VS)
        title({['neuron' num2str(n)] ; num2str(p{1}(n)) ; num2str(p{2}(n))})
        
    end
    
    %%
    f2 = figure;
    f1.Color = [1 1 1];
    for n = 1:nneurons
        plot(prc_gaba_sub{n}.PRC); hold on;
                set(gca,'Colormap', cbrewer('qual', 'Set2', 11))
    end
    
    
    %%
    
    

    f3 = figure;
    f3.Color = [1 1 1];
    for n = 1:nneurons
      
        ax3(n) = subplot(2,2,n)
        ttt = [1:500]- prc_gaba_sub{2}.peaktimes(1);
        waterfall(prc_gaba_sub{n}.VS)
        title({['neuron' num2str(n)] ; num2str(p{1}(n)) ; num2str(p{2}(n))})
        
    end
    
    %%
    f4 = figure;
    f4.Color = [1 1 1];
    for n = 1:nneurons
        plot(prc_gaba_sub{n}.PRC); hold on;
        set(gca,'Colormap', cbrewer('qual', 'Set2', 11))
    end
    

    
 %%
 
  %%
    f5 = figure;
    f5.Color = [1 1 1];
    for n = 1:nneurons
      
        ax5(n) = subplot(2,2,n)
        ttt = [1:500]- prc_ampa_sub{2}.peaktimes(1);
        waterfall(prc_ampa_sub{n}.VS)
        title({['neuron' num2str(n)] ; num2str(p{1}(n)) ; num2str(p{2}(n))})
        
    end
    %%
    
    f6 = figure; clf;
    f6.Color = [1 1 1];
    for n = 1:gabavals
        plot(prc_ampa_sub{n}.PRC); hold on;
    end
    
 
 %% perturbation GABA vs AMPA
	pert.mask  	  {1} = 1;
	pert.amplitude{1} = 1;
	pert.duration {1} = 1;
	pert.type	  {1} = 'ampa_dend';
	pert.triggers{1} = [0];   
   
    pert.mask  	  {2} = 1;
	pert.amplitude{2} = 1;
	pert.duration {2} = 4;
	pert.type	  {2} = 'gaba_dend';
	pert.triggers{2} = [30];   
   
    
    for n = 1:nneurons
        neuron{n}.gbar_gaba_dend = 1;
        prc_comb_sub{n} = IO_PRC(neuron{n},pert, npert);
    end
    
    f7 = figure;
    f7.Color = [1 1 1];
    for n = 1:nneurons
      
        ax7(n) = subplot(2,2,n)
        ttt = [1:500]- prc_comb_sub{2}.peaktimes(1);
        waterfall(prc_comb_sub{n}.VS)
        title({['neuron' num2str(n)] ; num2str(p{1}(n)) ; num2str(p{2}(n))})
        
    end
    
    %%
    f8 = figure;
    f8.Color = [1 1 1];
    for n = 1:nneurons
        plot(prc_gaba_sub{n}.PRC); hold on;
    end
    


%%

        %%
    f9 = figure;
    f9.Color = [1 1 1];
    for n = 1:gabavals
        plot(prc_gaba_pspace{n}.PRC); hold on;
    end
    
    
    %%
    
        f10 = figure;
    f10.Color = [1 1 1];
    for n = 1:gabavals
      
        plot(prc_gaba_pspace{n}.VS'); hold on

        
    end
    