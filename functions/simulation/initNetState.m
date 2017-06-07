function state = initNetState(noNeurons, use_gpu, tempState)
% create clones (with or without random mutation)

% if we receive a struct, load that one and return
if isstruct(tempState)
        disp('attempting to use tempstate!')
        for eachfield = fieldnames(tempState)'
                eval(['state.' eachfield{1} ' = tempState.' eachfield{1} ';']);
        end

        return

end


% compatibility
        state.ou_noise = zeros(noNeurons,1);
        state.curr_noise_pert = zeros(noNeurons,1);

try
    % gpuDeviceCount needs parallel toolbox
    if gpuDeviceCount == 0; use_gpu = 0; disp('no gpu found in this machine'); end
catch
    use_gpu = 0; disp('no parallel computing toolbox present');
end

% else, create neurons from scratch
randomize = 0;
if ~randomize

        %% Initial somatic parameters
        state.V_soma = -60*ones(noNeurons,1); % state.V_soma(1,1) = -100;
        state.Sodium_h = 0.3596066*ones(noNeurons,1);

        % Potassium
        state.Potassium_n   = 0.2369847*ones(noNeurons,1);
        state.Potassium_x_s = 0.1*ones(noNeurons,1);

        % Low-threshold calcium
        state.Calcium_k = 0.7423159*ones(noNeurons,1);
        state.Calcium_l = 0.0321349*ones(noNeurons,1);

        %% Initial dendritic parameters
        state.V_dend = -60*ones(noNeurons,1);

        % High-threshold calcium
        state.Calcium_r = 0.0112788*ones(noNeurons,1);

        % Calcium-dependent potassium
        state.Potassium_s = 0.0049291*ones(noNeurons,1);

        % H current
        state.Hcurrent_q = 0.0337836*ones(noNeurons,1);
        state.Hcurrent_q_s = 0.0337836*ones(noNeurons,1);

        % Calcium concentration
        state.Ca2Plus = 3.7152*ones(noNeurons,1);

        % Currents
        state.I_CaH = 0.5*ones(noNeurons,1);
        state.I_cx36 = zeros(noNeurons,1);
        state.curr_noise = zeros(noNeurons,1); % alloc
        state.current = zeros(noNeurons,1);
        state.vclamp     = zeros(noNeurons,1);



        %% Initial axonal parameters
        state.V_axon = -60*ones(noNeurons,1);

        % Sodium
        state.Sodium_m_a    = 0.003596066*ones(noNeurons,1);
        state.Sodium_h_a    = 0.9*ones(noNeurons,1);

        % Potassium
        state.Potassium_x_a  = 0.2369847*ones(noNeurons,1);

        state.Ca2_soma = 50e-9*ones(noNeurons,1);
        state.g_gaba_dend = 0*ones(noNeurons,1);
        state.g_gaba_soma = 0*ones(noNeurons,1);
        state.g_ampa_soma = 0*ones(noNeurons,1);
        state.g_ampa_dend = 0*ones(noNeurons,1);


        state.backgroundnoise = 0*ones(noNeurons,1);

        state.curr_noise_pert = zeros(noNeurons,1);
        state.ou_noise = zeros(noNeurons,1);
else
    
    
    
    % #TODO: randomize cell initial states


end


        
for eachfield = fieldnames(state)'
        if use_gpu
            eval(['state.' eachfield{1} ' = gpuArray(state.' eachfield{1} ');';]);
        else
            eval(['state.' eachfield{1} ' = state.' eachfield{1} ';']);
        end
end



