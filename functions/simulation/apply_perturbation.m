% apply_perturbation.m
% here we calculate the conductances for gabas and ampas as a function of stimuli
% function [g_AMPA g_GABA_dend g_GABA_soma Ca2_soma] = apply_perturbation(pert,t,clock_freq, g_AMPA, g_GABA_dend, g_GABA_soma, Ca2_soma)
function [g_AMPA g_AMPA_dend g_GABA_dend g_GABA_soma Ca2_soma currnoise] = apply_perturbation(pert,t,clock_freq, g_AMPA, g_AMPA_dend, g_GABA_dend, g_GABA_soma, Ca2_soma)

% # TODO introduce ampa in the dendrite

persistent AMPAGAUGE AMPAGAUGE_DEND GABAGAUGE_DEND GABAGAUGE_SOMA dAMPAS_dt dAMPAS_dend_dt dg_GABA_dend_dt dg_GABA_soma_dt dCa2_soma_dt INIT CURRNOISE RNGSTATE

Ca_0  = 2e-3;
noneurons = numel(g_AMPA);

if isempty(INIT)
    
    AMPAGAUGE                         = zeros(noneurons,1);
    AMPAGAUGE_DEND                    = zeros(noneurons,1);
    GABAGAUGE_DEND                    = zeros(noneurons,1);
    GABAGAUGE_SOMA                    = zeros(noneurons,1);
    CURRNOISE                         = zeros(noneurons,1);
    INIT = 1;                 
end                 
                  
                  
% pert.mask                 
% pert.amplitude
% pert.triggers
% pert.pulse_duration
% pert.type

nm = 0;
for P = 1:length(pert.type)

  if isempty(pert.type{P})
    continue
  end

  try
  switch pert.type{P}
    case 'ampa'

      % dirac triggers to release neurotransmitter according to mask.    
      % AMPA

      
        if ismember(t, pert.triggers{P}*clock_freq);
          AMPAGAUGE = max(AMPAGAUGE,  pert.mask{P}.*pert.duration{P}*clock_freq);
        end 
    
    case 'ampa_dend'

      % dirac triggers to release neurotransmitter according to mask.    
      % AMPA

      
        if ismember(t, pert.triggers{P}*clock_freq);
          AMPAGAUGE_DEND = max(AMPAGAUGE_DEND,  pert.mask{P}.*pert.duration{P}*clock_freq);
        end 
       
    

    case 'gaba_dend'
      % dirac triggers to release neurotransmitter according to mask.
      
        if ismember(t, pert.triggers{P}*clock_freq);
          GABAGAUGE_DEND = max(GABAGAUGE_DEND,  pert.mask{P}.*pert.duration{P}*clock_freq);
        end

    case 'gaba_soma'
    
        if ismember(t, pert.triggers{P}*clock_freq);
          GABAGAUGE_SOMA = max(GABAGAUGE_SOMA,  pert.mask{P}.*pert.duration{P}*clock_freq);
          
        end
    
    case 'ampa_noise'
 
      AMPAGAUGE = max(AMPAGAUGE,...
         pert.mask{P} .* pert.duration{P} .* (rand(size(pert.mask{P})) < pert.triggers{P}) * clock_freq);

    case 'ou_noise'
      nm = nm +1;

      if t >= pert.triggers{P}*clock_freq & t < pert.triggers{P}*clock_freq + pert.duration{P}*clock_freq;
        % disp('yeah')

        th = pert.param{P}(1);
        mu = pert.param{P}(2);
        sig = pert.param{P}(3);
        mix = pert.param{P}(4);

          CURRNOISE(:,nm) = CURRNOISE(:,nm) + th*(mu - CURRNOISE(:,nm)) * (1/clock_freq) .* pert.mask{P} + ...
                      + mix * sig * sqrt(1/clock_freq) * randn .* pert.mask{P}+...
                      + (1-mix) * sig * sqrt(1/clock_freq) * randn(size(pert.mask{P})) .* pert.mask{P};
  
      else
          CURRNOISE(:,nm) = zeros(noneurons,1);
      end

     case 'ou_noise_pos'
      nm = nm +1;

      if t >= pert.triggers{P}*clock_freq & t < pert.triggers{P}*clock_freq + pert.duration{P}*clock_freq;
        % disp('yeah')

        th = pert.param{P}(1);
        mu = pert.param{P}(2);
        sig = pert.param{P}(3);
        mix = pert.param{P}(4);

          CURRNOISE(:,nm) = CURRNOISE(:,nm) + th*(mu - CURRNOISE(:,nm)) * (1/clock_freq) .* pert.mask{P} + ...
                      + mix * sig * sqrt(1/clock_freq) * abs(randn) .* pert.mask{P}+...
                      + (1-mix) * sig * sqrt(1/clock_freq) * abs(randn(size(pert.mask{P}))) .* pert.mask{P};
  
      else
          CURRNOISE(:,nm) = zeros(noneurons,1);
      end


    case 'ou_noise_neg'
      nm = nm +1;

      if t >= pert.triggers{P}*clock_freq & t < pert.triggers{P}*clock_freq + pert.duration{P}*clock_freq;
        % disp('yeah')

        th = pert.param{P}(1);
        mu = pert.param{P}(2);
        sig = pert.param{P}(3);
        mix = pert.param{P}(4);

          CURRNOISE(:,nm) = CURRNOISE(:,nm) + th*(mu - CURRNOISE(:,nm)) * (1/clock_freq) .* pert.mask{P} + ...
                      + mix * sig * sqrt(1/clock_freq) * (-1)* abs(-randn) .* pert.mask{P} ...
                      + (1-mix) * sig * sqrt(1/clock_freq) * (-1)* abs(-randn(size(pert.mask{P}))) .* pert.mask{P};
  
      else
          CURRNOISE(:,nm) = zeros(noneurons,1);
      end



    otherwise

      disp('could not find perturbation type!')
      pert.type{P}

  end

    catch E
        disp('apply perturbation problem -- check mask size and perturbation triggers for appropriate size and orientation (should be a column vector).')
        keyboard
      
  end



end

% TODO: MAN AT WORK
% AMPAGAUGE = max(AMPAGAUGE,  pert.mask{P}.*pert.duration{P}*clock_freq);



   % compute ampa state
   if sum(g_AMPA) | sum(AMPAGAUGE)
       [dAMPAS_dt] = arrayfun(@computeAMPA_state, AMPAGAUGE > 1, g_AMPA);
       g_AMPA = g_AMPA + dAMPAS_dt*(1e-3)/clock_freq;
   end

 % compute ampa state DEND
   if sum(g_AMPA_dend) | sum(AMPAGAUGE_DEND)
       [dAMPAS_dend_dt] = arrayfun(@computeAMPA_state, AMPAGAUGE_DEND > 1, g_AMPA_dend);
      

       g_AMPA_dend = g_AMPA_dend + dAMPAS_dend_dt*(1e-3)/clock_freq;
   end


    % compute gaba state dendrite
    if sum(g_GABA_dend) | sum(GABAGAUGE_DEND)
       [dg_GABA_dend_dt] = arrayfun(@computefastGABA_state, GABAGAUGE_DEND>1 , g_GABA_dend);
       g_GABA_dend = g_GABA_dend + dg_GABA_dend_dt*1e-3/clock_freq; % *1e-3 : converting derivatives to ms (see computeGABA.m)
       
    end

    % compute gaba state soma
    if sum(g_GABA_soma)>eps | sum(GABAGAUGE_SOMA) 
       [dg_GABA_soma_dt, dCa2_soma_dt] = arrayfun(@computeslowGABA_state, GABAGAUGE_SOMA>1 , g_GABA_soma, Ca2_soma);
       g_GABA_soma = g_GABA_soma + dg_GABA_soma_dt*1e-3/clock_freq; % *1e-3 : converting derivatives to ms (see computeGABA.m)
       Ca2_soma   = Ca2_soma + dCa2_soma_dt*1e-3/clock_freq; % *1e-3 : converting derivatives to ms
    end

    
% if sum(AMPAGAUGE+GABAGAUGE_SOMA+GABAGAUGE_DEND)

%   imagesc([AMPAGAUGE GABAGAUGE_SOMA GABAGAUGE_DEND])
%   colorbar
%   sum(sum(AMPAGAUGE+GABAGAUGE_SOMA+GABAGAUGE_DEND))
%   pause(.1)

% end

%                    __      __                                     
%   __  ______  ____/ /___ _/ /____     ____ _____ ___  ______ ____ 
%  / / / / __ \/ __  / __ `/ __/ _ \   / __ `/ __ `/ / / / __ `/ _ \
% / /_/ / /_/ / /_/ / /_/ / /_/  __/  / /_/ / /_/ / /_/ / /_/ /  __/
% \__,_/ .___/\__,_/\__,_/\__/\___/   \__, /\__,_/\__,_/\__, /\___/ 
%     /_/                            /____/            /____/       

  AMPAGAUGE      = subplus(AMPAGAUGE      -1);
  AMPAGAUGE_DEND = subplus(AMPAGAUGE_DEND -1);
  GABAGAUGE_DEND = subplus(GABAGAUGE_DEND -1);
  GABAGAUGE_SOMA = subplus(GABAGAUGE_SOMA -1);
  currnoise = sum(CURRNOISE,2);

