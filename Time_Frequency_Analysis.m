addpath('Processes_EEG\Long')

Files = dir('.\Processes_EEG\Long\EEG*');
channels = { 'T7', 'T8'}; % EEG channels of interest


% Pre-allocate containers as a struct-of-fields by channel
res_evoked     = struct();  % res_evoked.T7(i).pow{k}
res_evoked_bas = struct();  % res_evoked_bas.T7(i).pow{k}
for c = 1:numel(channels)
    res_evoked.(channels{c})     = struct([]);
    res_evoked_bas.(channels{c}) = struct([]);
end

for i = 1:length(Files)
    S = load(Files(i).name);
    fn = fieldnames(S);
    data = S.(fn{1});
    fs   = data.fsample;

    % (Optional) normalize labels once
    % data.label = upper(strtrim(data.label));  % uncomment if helpful

    for t = 1:numel(channels)
        current_channel = channels{t};  % 'T7', then 'T8'

        % --- Preprocess: select just this channel ---
        cfg = [];
        cfg.channel        = {current_channel}; % select only this channel
        cfg.demean         = 'yes';
        cfg.baselinewindow = [-1 0];
        eegdata = ft_preprocessing(cfg, data);

        % Safety check: ensure we actually got the intended channel
        if ~any(strcmpi(current_channel, eegdata.label))
            warning('Subject %d: requested %s but selected labels: {%s}', ...
                i, current_channel, strjoin(eegdata.label, ', '));
            continue  % skip gracefully
        end

        for k = 1:6
            % ---- Time-frequency analysis ----
            cfg = [];
            cfg.method     = 'wavelet';
            cfg.output     = 'pow';
            cfg.pad        = 'nextpow2';
            cfg.foi        = 1:0.5:45;
            cfg.toi        = -1:0.01:4.5;
            cfg.width      = 7;
            cfg.trials     = find(eegdata.trialinfo(:,1) == k);
            cfg.keeptrials = 'no';
            tf = ft_freqanalysis(cfg, eegdata);

            % ---- Baseline ----
            cfg = [];
            cfg.baseline     = [-0.5 0];
            cfg.baselinetype = 'relative';
            tfb = ft_freqbaseline(cfg, tf);

            % ---- Store by channel name (no eval) ----
            res_evoked.(current_channel)(i).pow{k}     = tf;
            res_evoked_bas.(current_channel)(i).pow{k} = tfb;
        end
    end
end

% If you still want variables like res_evoked_T7, res_evoked_T8:
res_evoked_T7     = res_evoked.T7;
res_evoked_T8     = res_evoked.T8;
res_evoked_bas_T7 = res_evoked_bas.T7;
res_evoked_bas_T8 = res_evoked_bas.T8;


% Power spectrum averaged

channels = {'T7', 'T8'}; % EEG channels of interest

for ch = channels
    current_channel = ch{1}; % e.g. 't7', 't8'

    % Dynamically resolve the whole variable (once)
    res_evo     = eval(['res_evoked_' current_channel]);
    res_evo_bas = eval(['res_evoked_bas_' current_channel]);

    % Initialize containers
    gdavg_evo_channel      = cell(1,6);
    gdavg_evo_bas_channel  = cell(1,6);

    for k = 1:6
        cfg = [];
        cfg.keepindividual = 'no';

        % Extract inputs explicitly (like your working TestGA)
        
        in1 = res_evo(1).pow{k};
        in2 = res_evo(2).pow{k};
        in3 = res_evo(3).pow{k};

        gdavg_evo_channel{k} = ft_freqgrandaverage(cfg, in1, in2);

        b1 = res_evo_bas(1).pow{k};
        b2 = res_evo_bas(2).pow{k};
        b3 = res_evo_bas(3).pow{k};

        gdavg_evo_bas_channel{k} = ft_freqgrandaverage(cfg, b1, b2);
    end

    % Save results under channel-specific variable names
    assignin('base', ['gdavg_evo_' current_channel], gdavg_evo_channel);
    assignin('base', ['gdavg_evo_bas_' current_channel], gdavg_evo_bas_channel);
end


%% plot


% Loop over each channel
for ch = channels
    current_channel = ch{1}; % Extract channel name

    % Dynamically construct the variable name for the baseline-corrected grand average
    gdavg_evo_bas_var = ['gdavg_evo_bas_' current_channel];

       % Extract power spectrum for all conditions (1 to 6)
    for k = 1:6
        % Dynamically access the baseline-corrected grand average
        Power_corrected_channel(k, :, :) = eval([gdavg_evo_bas_var '{k}.powspctrm']);
    end

    % Save Power_corrected for this channel
    eval(['Power_corrected_' current_channel ' = Power_corrected_channel;']);
end

mean_power = (Power_corrected_T7 + Power_corrected_T8)/2;


%% Plot Baboon Data

Power_Baboon= zscore(squeeze(mean(mean_power(1,:,51:501), [1 3])));
Power_Human= zscore(squeeze(mean(mean_power(3,:,51:501), [1 3])));
Power_Noise= zscore(squeeze(mean(mean_power(5,:,51:501), [1 3])));


Power_Baboon_2= zscore(squeeze(mean(mean_power(2,:,51:501), [1 3])));
Power_Human_2= zscore(squeeze(mean(mean_power(4,:,51:501), [1 3])));
Power_Noise_2= zscore(squeeze(mean(mean_power(6,:,51:501), [1 3])));


time = gdavg_evo_T7{1, 1}.time;
freq = gdavg_evo_T7{1, 1}.freq;

% plot baboon

figure


sfh1= subplot(1,2,1);
imagesc(time, freq, squeeze(zscore(mean(mean_power(1,:,:), 1)))); axis xy; colormap(parula);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title('Baboon_1')
caxis([-2 3])

sfh2=subplot(1,2,2);
imagesc(time, freq, squeeze(zscore(mean(mean_power(2,:,:), 1)))); axis xy; colormap(parula);

xlabel('Time (s)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title ('Baboon_2')
caxis([-2 3])


c = colorbar('southoutside');
c.Label.String = 'relative power (zscore)';


figure

hold on
plot(Power_Baboon,freq, 'LineWidth', 3, 'Color', [225/255 147/255 14/255])
plot(Power_Baboon_2,freq, 'LineWidth', 3, 'Color', [42/255 144/255 173/255])

title ('Power spectrum')
ylim([1 45])
axis square
xlabel('Power (a.u.)')
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');

% plot human

figure


sfh1= subplot(1,2,1);
imagesc(time, freq, squeeze(zscore(mean(mean_power(3,:,:), 1)))); axis xy; colormap(parula);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title('Human_1')
caxis([-2 3])

sfh2=subplot(1,2,2);
imagesc(time, freq, squeeze(zscore(mean(mean_power(4,:,:), 1)))); axis xy; colormap(parula);

xlabel('Time (s)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title ('Human_2')
caxis([-2 3])


c = colorbar('southoutside');
c.Label.String = 'relative power (zscore)';

figure

hold on
plot(Power_Human,freq, 'LineWidth', 3, 'Color', [225/255 147/255 14/255])
plot(Power_Human_2,freq, 'LineWidth', 3, 'Color', [42/255 144/255 173/255])

title ('Power spectrum')
ylim([1 45])
axis square
xlabel('Power (a.u.)')
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');

 

% plot noise

figure


sfh1= subplot(1,2,1);
imagesc(time, freq, squeeze(zscore(mean(mean_power(5,:,:), 1)))); axis xy; colormap(parula);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title('Noise_1')
caxis([-2 3])

sfh2=subplot(1,2,2);
imagesc(time, freq, squeeze(zscore(mean(mean_power(6,:,:), 1)))); axis xy; colormap(parula);

xlabel('Time (s)')
xlim([-0.5 4.5])
ylim([1 45])
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');
axis square
title ('Noise_2')
caxis([-2 3])

figure

hold on
plot(Power_Noise,freq, 'LineWidth', 3, 'Color', [225/255 147/255 14/255])
plot(Power_Noise_2,freq, 'LineWidth', 3, 'Color', [42/255 144/255 173/255])

title ('Power spectrum')
ylim([1 45])
axis square
xlabel('Power (a.u.)')
set(gca,'FontName','Arial','FontSize',12);
set(gcf, 'Color', 'w');








