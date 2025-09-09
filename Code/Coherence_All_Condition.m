%% CaCoh
addpath('Processes_EEG\Long')

Files = dir('.\Processes_EEG\Long\EEG*');
channels = {'T7', 'T8'}; % Only T7 and T8 electrodes
reference_channel = 'env'; % Reference channel for coherence analysis

% Initialize an empty cell array to store data for the table
table_data = {};

for i = 1:length(Files)
    data = load(Files(i).name);
    name = fieldnames(data);
    data = data.(name{1});
    fs = data.fsample;

    for ch = channels % Loop over T7 and T8
        current_channel = ch{1}; % Extract channel name
        ch_number = find(strcmp(channels, current_channel)); % Get the index of the current channel

        % Preprocess EEG data: Select only the current pair of channels
        cfg = [];
        cfg.channel = {current_channel, reference_channel}; % Select pair: current channel and env
        cfg.demean = 'yes';
        cfg.baselinewindow = [-1 0];
        eegdata = ft_preprocessing(cfg, data); % Preprocess only the selected pair

        for k = 1:6 % Loop over stimulation types
            cfg = [];
            cfg.method = 'wavelet';
            cfg.output = 'powandcsd';
            cfg.pad = 'nextpow2';
            cfg.foi = 0.5:0.5:10; % Frequency range
            cfg.toi = 0.6:0.01:4; % Time of interest
            cfg.width = 4;
            cfg.keeptrials = 'no';
            cfg.trials = find(eegdata.trialinfo(:, 1) == k); % Select trials for the current stimulation type
            Rdata{i} = ft_freqanalysis(cfg, eegdata);

            % Coherence analysis for current channel vs reference
            cfg = [];
            cfg.method = 'coh';
            Rstat{i} = ft_connectivityanalysis(cfg, Rdata{i});

            % Extract coherence values
            if ~isempty(Rstat{i}.cohspctrm)
                coh_values = squeeze(mean(Rstat{i}.cohspctrm, 3, 'omitnan')); % Average across trials
                freqs = Rstat{i}.freq;

                % Add data to the table_data cell array
                for f = 1:length(freqs)
                    table_data = [table_data; {Files(i).name, current_channel, ch_number, k, freqs(f), coh_values(f)}];
                end
            end
        end
    end
end

% Convert the cell array to a table
table_headers = {'FileName', 'Electrode', 'ChannelNumber', 'ConditionNumber', 'Frequency', 'CaCoh'};
output_table = cell2table(table_data, 'VariableNames', table_headers);

% Save the table to a file
writetable(output_table, 'coherence_table.csv', 'Delimiter',';');