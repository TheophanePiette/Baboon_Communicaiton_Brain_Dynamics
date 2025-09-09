% Define base folders for all subjects
baseFolders = {
    'C:\Users\piette\Desktop\BABOON_EEG_STUDY\Feline\Rhythm_Processing\Pac_Standard_Aout_2025',
    'C:\Users\piette\Desktop\BABOON_EEG_STUDY\Formule\Rhythm_Processing\Pac_Standard_Aout_2025',
    'C:\Users\piette\Desktop\BABOON_EEG_STUDY\Ozone\Rhythm_Processing\Pac_Standard_Aout_2025'
};

% Condition to analyze (ONLY Baboon_1)
conditions = 'Human_2';

% Electrode to analyze
electrode = 'T7'; 

% Initialize PAC matrix storage for T7
allPAC_matrices = [];

% Define frequency axes
phaseFreq = 2:0.5:10; % 17 phase frequencies
ampFreq = 20:1:45;  % 61 amplitude frequencies

% Loop through subjects
for subjIdx = 1:length(baseFolders)
    baseFolder = baseFolders{subjIdx};
    subjectData = [];  % Store PAC matrices for this subject (T7 only)
    
    % Loop through the selected condition (only Baboon_1)
    condFolder = fullfile(baseFolder, conditions);
    elecFolder = fullfile(condFolder, electrode);
    
    % Check if folder exists
    if ~isfolder(elecFolder)
        disp(['Folder does not exist: ', elecFolder]);
        continue;
    end
    
    % Load all Z_score.mat files in this folder
    files = dir(fullfile(elecFolder, 'Z_score.mat'));
    for fileIdx = 1:length(files)
        data = load(fullfile(files(fileIdx).folder, files(fileIdx).name));
        
        % Ensure Z_Score exists and has correct dimensions
        if isfield(data, 'Z_Score') && isequal(size(data.Z_Score), [17, 26])
            subjectData = cat(3, subjectData, data.Z_Score);
        else
            disp(['Invalid or missing Z_Score in file: ', files(fileIdx).name]);
        end
    end
    
    % Compute mean PAC matrix for this subject
    if ~isempty(subjectData)
        meanPAC_subject = mean(subjectData, 3);
        allPAC_matrices = cat(3, allPAC_matrices, meanPAC_subject); 
    end
end

% Compute the grand average PAC matrix for T7 across all subjects (Baboon_1 only)
grandMeanPAC_T7 = mean(allPAC_matrices, 3); 

% ----- Plot Heatmap -----
figure;
imagesc(phaseFreq,ampFreq, grandMeanPAC_T7');
set(gca, 'YDir', 'normal');
ylabel('Amplitude Frequency (Hz)');
xlabel('Phase Frequency (Hz)');
title('Grand Mean Z-Score PAC (T7, Baboon_1)');
colorbar;
caxis([-2.5 2.5]); % Set colorbar limits

% ----- Interpolate for Smoother Heatmap -----
[x, y] = meshgrid(1:size(grandMeanPAC_T7, 2), 1:size(grandMeanPAC_T7, 1)); 
[xq, yq] = meshgrid(linspace(1, size(grandMeanPAC_T7, 2), 200), linspace(1, size(grandMeanPAC_T7, 1), 200)); 

% Interpolate data using cubic interpolation
smoothedPAC = interp2(x, y, grandMeanPAC_T7, xq, yq, 'cubic'); 

% Plot smoothed heatmap
figure;
imagesc(phaseFreq,ampFreq, smoothedPAC');
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Amplitude Frequency (Hz)');
ylabel('Phase Frequency (Hz)');
title('Smoothed PAC Heatmap (T7, Baboon_1)');
colormap(parula); 
caxis([-2 3]); 


% ----- Cluster Analysis (1.96 / 1.75 threshold) -----

% Define thresholds
seedThreshold = 1.96;

% Step 1: create a binary mask for Z >= 1.75 (extension threshold)
extensionMask = grandMeanPAC_T7 >= seedThreshold;

% Step 2: label all connected regions using 8-connectivity
CC = bwconncomp(extensionMask, 8); % returns connected components

% Initialize cluster storage
validClusters = struct('PhaseFreqs', {}, 'AmpFreqs', {}, 'ZScores', {}, 'MaxZ', {}, 'Size', {});

clusterCount = 0;

% Step 3: Check each cluster for seed points (Z >= 1.96)
for i = 1:CC.NumObjects
    idx = CC.PixelIdxList{i}; % linear indices of the cluster
    zvals = grandMeanPAC_T7(idx);
    
    if any(zvals >= seedThreshold)
        % Convert linear indices to subscript (row/col) to get freq values
        [rows, cols] = ind2sub(size(grandMeanPAC_T7), idx);
        
        clusterCount = clusterCount + 1;
        validClusters(clusterCount).PhaseFreqs = phaseFreq(rows);
        validClusters(clusterCount).AmpFreqs = ampFreq(cols);
        validClusters(clusterCount).ZScores = zvals;
        validClusters(clusterCount).MaxZ = max(zvals);
        validClusters(clusterCount).Size = length(zvals);
    end
end

% ----- Display detected clusters -----
fprintf('\n=== Detected PAC Clusters (Seed ≥ %.2f, Extend ≥ %.2f) ===\n', seedThreshold, extendThreshold);
for i = 1:length(validClusters)
    cl = validClusters(i);
    fprintf('Cluster %d: Size = %d, Max Z = %.2f\n', i, cl.Size, cl.MaxZ);
    for j = 1:length(cl.ZScores)
        fprintf('  Phase: %.2f Hz | Amplitude: %.2f Hz | Z = %.2f\n', ...
            cl.PhaseFreqs(j), cl.AmpFreqs(j), cl.ZScores(j));
    end
end

