function [SA DA MA] =  plot_current_contributions(sim)
   % obviously only works if we saved all currents to the simulation results

plot_individual_neurons = 0;
plot_individual_neurons_area = 1;

% sim = st_st;

switch sim.cellFunction
   case 'devel'
      %                    1      2         3         4           5               6               7           8             9           10        11       12         13      14         15        16    17      18        19     20          21       
      currents =     {'V_soma','V_dend','V_axon', 'I_CaL',       'I_ds',         'I_as',         'I_Na_s',  'I_ls',      'I_Kdr_s',  'I_K_s',   'I_h_s' , 'I_CaH',   'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'}; 
      currentnames = {'V_soma','V_dend','V_axon', 'Calcium low', 'leak to dend', 'leak to axon', 'Na soma', 'soma leak', 'Kr soma ', 'Ks soma', 'Ih_soma','Ca High', 'dend-soma leak', 'dend leak', 'Ca act K (SK)', 'Cx36', 'I h', 'K axon', 'axon-soma leak', 'axon leak', 'Na axon'};

      somatic_currents = [4:11];
      dendritic_currents = [12:17];
      axon_currents      = [18:21]
   case 'vanilla'
%                         1      2         3         4           5               6               7           8             9           10        11       12         13      14         15        16    17      18        19     20          21       
      currents =     {'V_soma','V_dend','V_axon', 'I_CaL',       'I_ds',         'I_as',         'I_Na_s',  'I_ls',      'I_Kdr_s',  'I_K_s',  'I_CaH',   'I_sd', 'I_ld', 'I_K_Ca', 'I_cx36', 'I_h', 'I_K_a', 'I_sa', 'I_la', 'I_Na_a'}; 
      currentnames = {'V_soma','V_dend','V_axon', 'Calcium low', 'leak to dend', 'leak to axon', 'Na soma', 'soma leak', 'Kr soma ', 'Ks soma', 'Ca High', 'dend-soma leak', 'dend leak', 'Ca act K (SK)', 'Cx36', 'I h', 'K axon', 'axon-soma leak', 'axon leak', 'Na axon'};

      somatic_currents = [4:10];
      dendritic_currents = [11:16];
      axon_currents      = [17:20]
end
try; set(0,'defaultaxescolororder', linspecer(8)); catch set(0,'defaultaxescolororder', jet(8)); end

% 4:10 somatic

Plist = sim.Plist;
neurons = [1:size(Plist,1)];
neurons = neurons(end-5:end);


 for neuron = neurons
   clf

      f = 0;
      for i = somatic_currents
         f = f+1;
         eval([ 'SA{neuron}(' num2str(f) ', :) = sim.networkHistory.' currents{i} '(' num2str(neuron) ',:);']   )
      end
      
      f = 0;
      for i = dendritic_currents
         f = f+1;
         eval([ 'DA{neuron}(' num2str(f) ', :) = sim.networkHistory.' currents{i} '(' num2str(neuron) ',:);']   )
      end
      
      f = 0;
      for i = axon_currents
         f = f+1;
         eval([ 'AA{neuron}(' num2str(f) ', :) = sim.networkHistory.' currents{i} '(' num2str(neuron) ',:);']   )
      end
      
end


if plot_individual_neurons
   figure
   for neuron = neurons
   clf

   
   	a(1) = subplot(3,1,1); plot(real(SA{neuron})'), legend(currentnames{somatic_currents}), hold on, plot(zscore(sim.networkHistory.V_soma(neuron,:)), 'linewidth', 2, 'color', 'r')
   	title({ ['neuron: ' num2str(neuron)] ; ['CaT, Kca, Ih, Vh' ] ; num2str(Plist(neuron,:))})
      ylim([-50 50])

   	f = 0;
   

   	a(2) = subplot(3,1,2); plot(real(DA{neuron})'), legend(currentnames{dendritic_currents}), hold on, plot(zscore(sim.networkHistory.V_dend(neuron,:)), 'linewidth', 2,'color', 'r')
      ylim([-50 50])

   	f = 0;
   	
   	a(3) = subplot(3,1,3); plot(real(AA{neuron})'), legend(currentnames{axon_currents}), hold on, plot(zscore(sim.networkHistory.V_axon(neuron,:)), 'linewidth', 2, 'color','r')
   	colormap(linspecer(11))
   	xlabel('ms')
      ylabel('Current - uA/cm^2')

   	linkaxes(a,'x')
      ylim([-50 50])

   	pause

   end

end


if plot_individual_neurons_area
figure
   for neuron = neurons
   clf

      f = 0;
   a(2) = subplot(3,1,1); 
      bla1 = subplus(real(SA{neuron}));
      bla2 = -subplus(real(-SA{neuron}));
      area(bla1'),hold on
      area(bla2'),
       legend(currentnames{somatic_currents}), hold on, plot(zscore(sim.networkHistory.V_soma(neuron,:)), 'linewidth', 2, 'color', 'r')
      title({ ['neuron: ' num2str(neuron)] ; ['CaT, Kca, Ih, Vh' ] ; num2str(Plist(neuron,:))})

      
      a(2) = subplot(3,1,2); 
      bla1 = subplus(real(DA{neuron}));
      bla2 = -subplus(real(-DA{neuron}));
      area(bla1'),hold on
      area(bla2'),
      legend(currentnames{dendritic_currents}),
       hold on, plot(zscore(sim.networkHistory.V_dend(neuron,:)), 'linewidth', 2,'color', 'r')

      
      a(3) = subplot(3,1,3); 
      bla1 = subplus(real(AA{neuron}));
      bla2 = -subplus(real(-AA{neuron}));
      area(bla1'),hold on
      area(bla2'),
      legend(currentnames{axon_currents}),  plot(zscore(sim.networkHistory.V_axon(neuron,:)), 'linewidth', 2, 'color','r')

      colormap(linspecer(10))
      xlabel('ms')

      linkaxes(a,'x')

      pause

   end

end

