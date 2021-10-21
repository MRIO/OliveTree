
if 1
npert = 16;
    
    %% Make neurons

p1 = [.9 1.2]; 		% CalciumL - conductance range
p2 = [1.5 3.5];  		% g_CaK (does not impact PRC without CaH.)
% p2 = [1.2 2.4];  		% g_CaK (does not impact PRC without CaH.)

[p{1} p{2} ] = ndgrid(p1,p2);

nneurons = prod(size(p{1}));
        
    for n = 1:nneurons
        neuron{n} = createDefaultNeurons(1, 'celltypes',  'clones');
        neuron{n}.g_CaL     = p{1}(n);
        neuron{n}.g_CaH 	= p{2}(n);
%         neuron{n}.g_h 	= p{2}(n);
    end

    
    
    
%% template for perturbation
    pert.mask  	  {1} = 1;
	pert.amplitude{1} = 1;
	pert.duration {1} = 1;
	pert.triggers{1} = [100];   
%%
    
end

if 0
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
    pert.duration {1} = 1;
	pert.type	  {1} = 'gaba_dend';
	
    for n = 1:nneurons
        neuron{n}.gbar_gaba_dend = 1;
        prc_gaba_sub{n} = IO_PRC(neuron{n},pert, npert);
    end


%% perturbation GABA DENDRITE PSPACE
	pert.duration {1} = 2;
	pert.type	  {1} = 'gaba_dend';
	
    gabavals = 10;
     
    for n = 1:gabavals
        neuron{1}.gbar_gaba_dend =(n-1)*.4;
        prc_gaba_pspace{n} = IO_PRC(neuron{1},pert, npert);
    end
%%

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

end

%% PLOTS

%%
if 1


    f1 = figure;
    f1.Color = [1 1 1];
    for n = 1:nneurons
        ax(n) = subplot(2,2,n)
        ttt = [1:500]- prc_gaba_sub{n}.peaktimes(1);
        waterfall(prc_gaba_sub{n}.VS)
        title({['gaba sub']; ['neuron' num2str(n)] ; ['CaL: ' num2str(p{1}(n)) ]; ['Ih: ' num2str(p{2}(n))]})
    end

    %%
    nneurons = 4
    f2 = figure;
    f2.Color = [1 1 1];
    for n = 1:nneurons
    plot(prc_gaba_sub{n}.PRC(1,:), prc_gaba_sub{n}.PRC(2,:)); hold on;
    set(gca,'Colormap', cbrewer('qual', 'Set2', 11))
    title({['gaba sub PRC']})
    end
    axis tight
    xlabel('rad')

    legend({num2str([1:4]')})

    %%
    f3 = figure;
    f3.Color = [1 1 1];
    for n = 1:nneurons

        ax3(n) = subplot(2,2,n)
        ttt = [1:500]- prc_gaba_supra{n}.peaktimes(1);
        waterfall(prc_gaba_supra{n}.VS)
        title({['gaba supra']; num2str(p{1}(n)) ; num2str(p{2}(n))})

    end

    %%

    f4 = figure;
    f4.Color = [1 1 1];
    for n = 1:nneurons
    plot(prc_gaba_supra{n}.PRC(1,:), prc_gaba_supra{n}.PRC(2,:)); hold on;
    set(gca,'Colormap', cbrewer('qual', 'Set2', 11))
    title('PRC gaba supra')
    end
    legend({num2str([1:4]')})
    axis tight

    %%
    f5 = figure;
    f5.Color = [1 1 1];
    for n = 1:nneurons

        ax5(n) = subplot(2,2,n)
        ttt = [1:500]- prc_ampa_sub{n}.peaktimes(1);
        waterfall(prc_ampa_sub{n}.VS)
        title({['ampa sub' ] ; num2str(p{1}(n)) ; num2str(p{2}(n))})
    end

    %%
    f6 = figure; clf;
    f6.Color = [1 1 1];
    for n = 1:nneurons
        plot(prc_ampa_sub{n}.PRC(1,:), prc_ampa_sub{n}.PRC(2,:)); hold on;
    end
    title('PRC ampa sub')
    legend({num2str([1:4]')})
    axis tight

    %%
    f7 = figure;
    f7.Color = [1 1 1];
    for n = 1:nneurons

        ax7(n) = subplot(2,2,n)
        ttt = [1:500]- prc_ampa_supra{n}.peaktimes(1);
        waterfall(prc_ampa_supra{n}.VS)
        title({['ampa supra' ] ; num2str(p{1}(n)) ; num2str(p{2}(n))})

    end
    %%
    f8 = figure; clf;
    f8.Color = [1 1 1];
    for n = 1:nneurons
        plot(prc_ampa_supra{n}.PRC(1,:), prc_ampa_supra{n}.PRC(2,:)); hold on;
    end
    title('PRC ampa supra')
    legend({num2str([1:4]')})
    axis tight

    %%

    f9 = figure;
    f9.Color = [1 1 1];
    for n = 1:nneurons

        ax9(n) = subplot(2,2,n)
        ttt = [1:500]- prc_comb_sub{n}.peaktimes(1);
        waterfall(prc_comb_sub{n}.VS)
        title({['PRC combined'] ; num2str(p{1}(n)) ; num2str(p{2}(n))})

    end

    %%

    %%
    f10 = figure;
    f10.Color = [1 1 1];
    for n = 1:gabavals
    plot(prc_gaba_pspace{n}.PRC(1,:),prc_gaba_pspace{n}.PRC(2,:)); hold on;
    end
    legend(num2str(([1:gabavals]'-1)*.4))
    title('gaba prc pspace (stim=2ms)')
    axis tight
    %%

    f11 = figure;
    f11.Color = [1 1 1];
    for n = 1:gabavals
        subplot(1,10,n)
        imagesc(prc_gaba_pspace{n}.VS', [-70 -40]); hold on
    end
    title('gaba prc pspace (stim=2ms)')

    %%
    f12 = figure;
    f12.Color = [1 1 1];
    for n = 1:gabavals
        plot(prc_gaba_pspace{n}.PRC(1,:),prc_gaba_pspace{n}.PRC(2,:)); hold on;
    end
    title('gaba prc pspace (stim=2ms)')
    legend(num2str(([1:gabavals]'-1)*.4))
    axis tight


%     %%
%     f13 = figure;
%     f13.Color = [1 1 1];
%     for n = 1:gabavals
%     subplot(1,10,n)
%     imagesc(prc_gaba_pspace_d4{n}.VS', [-70 -40]); hold on
%     end
%     title('gaba prc pspace (stim=4ms)')

end
