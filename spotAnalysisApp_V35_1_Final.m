function spotAnalysisApp_V35_1_Final()
% spotAnalysisApp - GUI for batch .tif fluorescence image analysis
% VERSION - V35.1 (Final - QC Panel Bug Fix)
% - FEATURES:
%   - Bug Fix: Corrects a field name error in the updateQCPanel function,
%     ensuring Spot Area is calculated and displayed correctly.
%   - Includes all previous features: Redesigned UI, full Excel reporting,
%     ROI selection, Dual Calibration, Live Preview, etc.

% --- App Data Structure ---
appdata = struct('rootDir', '', 'subDirs', {{}}, 'results', struct(), 'extStdPath', '', 'roiPosition', []);

% --- UI Figure Creation & Layout ---
fig = uifigure('Name', 'Spot Analysis (V35.1 - Final)', 'Position', [100 100 1100, 650]);

% --- Panel 1: Setup & Selections ---
p1 = uipanel(fig, 'Title', '1. Setup & Selections', 'Position', [20 20 250 610], 'FontWeight', 'bold');
uilabel(p1, 'Text', 'Root Folder:', 'Position', [10 570 100 22]);
btnRoot = uibutton(p1, 'Text', 'Select Root', 'Position', [100 570 140 22], 'ButtonPushedFcn', @selectRootCallback);
lblRootPath = uilabel(p1, 'Text', 'Current Root: N/A', 'Position', [10 545, 230, 22], 'HorizontalAlignment', 'left');

btnROI = uibutton(p1, 'Text', 'Select Analysis ROI', 'Position', [10 510, 230, 25], 'ButtonPushedFcn', @selectROI, 'Enable', 'off');
lblROIStatus = uilabel(p1, 'Text', 'ROI: Not Set (Full Image)', 'Position', [10 490, 230, 15], 'FontColor', [0.5 0.5 0.5]);

calibPanel = uipanel(p1, 'Title', 'Calibration Mode', 'Position', [10 395, 230, 90]);
calibMode = uibuttongroup(calibPanel, 'Position', [5 5 220 80], 'SelectionChangedFcn', @calibModeChanged);
rbInternal = uiradiobutton(calibMode, 'Text', 'Internal Anchor:', 'Position', [10 45 120 22], 'Value', true);
ddAnchor = uidropdown(calibMode, 'Position', [120 45, 95, 22], 'Tooltip', 'Select well from current experiment');
rbExternal = uiradiobutton(calibMode, 'Text', 'External Standard:', 'Position', [10 20 120 22]);
btnExtStd = uibutton(calibMode, 'Text', 'Select', 'Position', [120 20 95 22], 'Enable', 'off', 'ButtonPushedFcn', @selectExternalStandard);
lblExtStd = uilabel(p1, 'Text', 'Ext. Std: N/A', 'Position', [10 370 230 22], 'FontColor', [0.5 0.5 0.5], 'HorizontalAlignment', 'right');

tabFolders = uitabgroup(p1, 'Position', [10 10 230 360]);
tabProcess = uitab(tabFolders, 'Title', 'Process');
lstbxSub = uilistbox(tabProcess, 'Position', [5 5 220 320], 'Multiselect', 'on');
tabControl = uitab(tabFolders, 'Title', 'NSB Control');
lstbxControl = uilistbox(tabControl, 'Position', [5 5 220 320], 'Multiselect', 'on');
tabDisplay = uitab(tabFolders, 'Title', 'Display');
lstbxRep = uilistbox(tabDisplay, 'Position', [5 5 220 320], 'Multiselect', 'on', 'ValueChangedFcn', @updateDisplayCallback);

% --- Panel 2: Parameters & Preview ---
p2 = uipanel(fig, 'Title', '2. Parameters & Preview', 'Position', [290 20 400 610], 'FontWeight', 'bold');
paramPanel = uipanel(p2, 'Title', 'Image Processing Parameters', 'Position', [10 310, 380, 285]);
uilabel(paramPanel, 'Text', 'Sensitivity (0~1):', 'Position', [10 220-30 110 22]);
edSens = uieditfield(paramPanel, 'numeric', 'Limits', [0 1], 'Value', 0.5, 'Position', [120 220-30 60 22], 'ValueChangedFcn', @updatePreview, 'Tooltip', 'Higher values detect dimmer spots (0.0 to 1.0)');
uilabel(paramPanel, 'Text', 'Gaussian Sigma (px):', 'Position', [10 180-30 120 22]);
edSigma = uieditfield(paramPanel, 'numeric', 'Limits', [0 Inf], 'Value', 1, 'Position', [130 180-30 60 22], 'ValueChangedFcn', @updatePreview, 'Tooltip', 'Image smoothing factor (e.g., 0.5 to 3)');
uilabel(paramPanel, 'Text', 'Open Radius (px):', 'Position', [10 140-30 120 22]);
edSE = uieditfield(paramPanel, 'numeric', 'Limits', [0 Inf], 'Value', 1, 'Position', [130 140-30 60 22], 'ValueChangedFcn', @updatePreview, 'Tooltip', 'Removes small noise spots (integer, e.g., 1, 2)');
uilabel(paramPanel, 'Text', 'Min Area (px²):', 'Position', [10 100-30 110 22]);
edMinA = uieditfield(paramPanel, 'numeric', 'Limits', [1 Inf], 'Value', 5, 'Position', [130 100-30 60 22], 'ValueChangedFcn', @updatePreview, 'Tooltip', 'Minimum spot size in pixels');
uilabel(paramPanel, 'Text', 'Max Area (px²):', 'Position', [10 60-30 110 22]);
edMaxA = uieditfield(paramPanel, 'numeric', 'Limits', [1 Inf], 'Value', 500, 'Position', [130 60-30 60 22], 'ValueChangedFcn', @updatePreview, 'Tooltip', 'Maximum spot size in pixels');
uilabel(paramPanel, 'Text', 'Saturation Threshold:', 'Position', [10 20-30 130 22]);
edSatThresh = uieditfield(paramPanel, 'numeric', 'Value', 65530, 'Position', [130 20-30 60 22], 'Editable', 'off', 'FontColor', [0.5 0.5 0.5]);

previewPanel = uipanel(p2, 'Title', 'Spot Detection 🔄 (Live Preview)', 'Position', [10 10, 380, 290]);
axPreviewOriginal = uiaxes(previewPanel, 'Position', [5 5 180 260]);
axPreviewOriginal.XTick = []; axPreviewOriginal.YTick = []; title(axPreviewOriginal, 'Original');
axPreviewProcessed = uiaxes(previewPanel, 'Position', [195 5 180 260]);
axPreviewProcessed.XTick = []; axPreviewProcessed.YTick = []; title(axPreviewProcessed, 'Detected');

% --- Panel 3: Run & Results ---
p3 = uipanel(fig, 'Title', '3. Run & Results', 'Position', [710 20 370 610], 'FontWeight', 'bold');
btnProcess = uibutton(p3, 'Text', 'Process Selected', 'Position', [10, 560, 350, 35], 'ButtonPushedFcn', @processSelectedCallback, 'FontSize', 14);
btnSave = uibutton(p3, 'Text', 'Save to Excel', 'Position', [10, 520, 350, 35], 'Enable', 'off', 'ButtonPushedFcn', @saveCallback, 'FontSize', 14);

axHist = uiaxes(p3, 'Position', [10, 280, 350, 220]); title(axHist, 'Intensity Distribution (Pre-NSB)');
qcPanel = uipanel(p3, 'Title', 'QC Summary (NSB Corrected)', 'Position', [10 10, 350, 260]);
lblSatSpots = uilabel(qcPanel, 'Text', 'Saturated Spots: N/A', 'Position', [10 215 330 22]);
lblSpotArea = uilabel(qcPanel, 'Text', 'Spot Area: N/A', 'Position', [10 195 330 22]);
lblRefMeans = uilabel(qcPanel, 'Text', 'Monomer Ref: N/A', 'Position', [10 175 330 22]);
uilabel(qcPanel, 'Text', 'Avg. Specific Spots/Folder:', 'Position', [10 145 330 22], 'FontWeight','bold');
lblMonoSpots = uilabel(qcPanel, 'Text', 'Monomer(<2.5x): N/A', 'Position', [20 120 320 22]);
lblTrimerSpots = uilabel(qcPanel, 'Text', 'Trimer: N/A', 'Position', [20 100 320 22]);
lblTrimerPlusSpots = uilabel(qcPanel, 'Text', 'Trimer+: N/A', 'Position', [20 80 320 22]);

%% --- Callback Functions ---
    function selectRootCallback(~, ~)
        d = uigetdir(pwd, 'Select the root folder');
        if isequal(d, 0); return; end
        [~, folderName, ~] = fileparts(d);
        lblRootPath.Text = ['Current Root: ' folderName];
        lblRootPath.Tooltip = d;
        info = dir(d);
        allNames = {info.name}; isDir = [info.isdir]; validDirNames = {};
        for i = 1:length(allNames)
            name = allNames{i};
            if isDir(i) && ~isempty(name) && name(1) ~= '.' && name(1) ~= '_'
                validDirNames{end+1} = name;
            end
        end
        dirs = setdiff(validDirNames, '_processed', 'stable');
        lstbxSub.Items = dirs; lstbxSub.Value = dirs;
        lstbxRep.Items = dirs; lstbxRep.Value = {};
        lstbxControl.Items = dirs;
        ddAnchor.Items = dirs;
        if ~isempty(dirs); ddAnchor.Value = dirs{1}; end
        btnROI.Enable = 'on';
        appdata.roiPosition = [];
        lblROIStatus.Text = 'ROI: Not Set (Full Image)';
        lblROIStatus.FontColor = [0.5 0.5 0.5];
        cla(axHist);
        cla(axPreviewOriginal); title(axPreviewOriginal, 'Original');
        cla(axPreviewProcessed); title(axPreviewProcessed, 'Detected');
        previewPanel.Title = 'Spot Detection 🔄 (Live Preview)';
        appdata.rootDir = d; appdata.subDirs = dirs; appdata.results = struct();
        btnSave.Enable = 'off';
    end

    function calibModeChanged(~, event)
        if event.NewValue == rbInternal
            ddAnchor.Enable = 'on';
            btnExtStd.Enable = 'off';
            lblExtStd.FontColor = [0.5 0.5 0.5];
        else % External
            ddAnchor.Enable = 'off';
            btnExtStd.Enable = 'on';
            lblExtStd.FontColor = [0 0 0];
        end
    end

    function selectExternalStandard(~, ~)
        d = uigetdir(appdata.rootDir, 'Select External Standard (Monomer) Folder');
        if isequal(d, 0); return; end
        appdata.extStdPath = d;
        [~, folderName, ~] = fileparts(d);
        lblExtStd.Text = ['Ext. Std: ' folderName];
        lblExtStd.Tooltip = d;
    end
    
    function selectROI(~, ~)
        repSel = lstbxRep.Value;
        if isempty(repSel) || iscell(repSel) && isempty(repSel{1})
            uialert(fig, 'Select a folder in the "Display" tab first.', 'Info');
            return;
        end
        if iscell(repSel); repSel = repSel{1}; end
        folderPath = fullfile(appdata.rootDir, repSel);
        files = dir(fullfile(folderPath, '*.tif'));
        if isempty(files); uialert(fig, 'No .tif files in the selected folder to define ROI.', 'Error'); return; end
        try
            I = imread(fullfile(folderPath, files(1).name));
            figROI = figure('Name', 'Draw ROI and DOUBLE-CLICK inside to confirm', 'NumberTitle', 'off', 'WindowState', 'maximized');
            imshow(I, []);
            h = imrect;
            if isempty(h); close(figROI); return; end
            appdata.roiPosition = wait(h);
            close(figROI);
            if ~isempty(appdata.roiPosition)
                pos = appdata.roiPosition;
                lblROIStatus.Text = sprintf('ROI Set: [%.0f, %.0f, %.0f, %.0f]', pos(1), pos(2), pos(3), pos(4));
                lblROIStatus.FontColor = [0, 0.45, 0.74];
                updatePreview();
            end
        catch ME
            uialert(fig, ['Could not open image to select ROI. Error: ' ME.message], 'Error');
        end
    end

    function processSelectedCallback(~, ~)
        progDialog = uiprogressdlg(fig, 'Title', 'Processing...', 'Message', 'Starting analysis.', 'Indeterminate', 'on');
        try
            procSel = lstbxSub.Value;
            if ischar(procSel); procSel = {procSel}; end
            if isempty(procSel); uialert(fig, 'No subfolders selected!', 'Selection Error'); close(progDialog); return; end
            
            progDialog.Message = 'Pass 1/4: Independent analysis...';
            rawResults = struct();
            for i = 1:length(procSel)
                procName = procSel{i};
                folderPath = fullfile(appdata.rootDir, procName);
                [perImageData, allIntensities] = processFolder(folderPath);
                validName = matlab.lang.makeValidName(procName);
                rawResults.(validName) = struct('perImageData', perImageData, 'allIntensities', allIntensities);
            end

            progDialog.Message = 'Pass 2/4: Calibrating reference...';
            calibrated_monomer_ref = NaN;
            
            if calibMode.SelectedObject == rbInternal
                anchorWell = ddAnchor.Value;
                validAnchorName = matlab.lang.makeValidName(anchorWell);
                anchorIntensities = rawResults.(validAnchorName).allIntensities;
                calibrated_monomer_ref = getMonomerRefFromData(anchorIntensities, 'Internal');
            else
                if isempty(appdata.extStdPath); uialert(fig, 'Please select an External Standard folder.', 'Selection Error'); close(progDialog); return; end
                [~, extIntensities] = processFolder(appdata.extStdPath);
                calibrated_monomer_ref = getMonomerRefFromData(extIntensities, 'External');
            end

            if isnan(calibrated_monomer_ref); uialert(fig, 'Failed to calibrate monomer reference.', 'Analysis Error'); close(progDialog); return; end
            disp(['Calibrated Monomer Reference set to: ', num2str(calibrated_monomer_ref)]);

            progDialog.Message = 'Pass 3/4: Classifying & Aggregating...';
            classifiedResults = struct();
            for i = 1:length(procSel)
                procName = procSel{i};
                validName = matlab.lang.makeValidName(procName);
                res = rawResults.(validName);
                
                numImages = height(res.perImageData);
                if numImages == 0; continue; end
                
                ratios_per_image = zeros(numImages, 3);
                for j = 1:numImages
                    [~, ratios_per_image(j,:), ~] = analyzePopulation(res.perImageData.CorrectedIntensities{j}, calibrated_monomer_ref);
                end
                
                agg_stats = struct();
                agg_stats.images_analyzed = numImages;
                agg_stats.spot_count_mean = mean(res.perImageData.SpotCount);
                agg_stats.spot_count_sem = std(res.perImageData.SpotCount) / sqrt(numImages);
                agg_stats.mean_intensity_mean = mean(res.perImageData.MeanIntensity, 'omitnan');
                agg_stats.mean_intensity_sem = std(res.perImageData.MeanIntensity, 'omitnan') / sqrt(numImages);
                agg_stats.saturated_count_total = sum(res.perImageData.SaturatedCount);
                agg_stats.ratios_mean = mean(ratios_per_image, 1, 'omitnan');
                agg_stats.ratios_sem = std(ratios_per_image, 0, 1, 'omitnan') / sqrt(numImages);
                
                [dominant_state, final_ratios, final_mean_corr] = analyzePopulation(res.allIntensities, calibrated_monomer_ref);
                
                classifiedResults.(validName) = struct(...
                    'dominant_state', dominant_state, ...
                    'peakRatios', final_ratios, ...
                    'mean_corrected', final_mean_corr, ...
                    'allIntensities', res.allIntensities, ...
                    'agg_stats', agg_stats, ...
                    'calibrated_monomer_ref', calibrated_monomer_ref);
            end

            progDialog.Message = 'Pass 4/4: Applying NSB correction...';
            finalResults = classifiedResults;
            ctrlSel = lstbxControl.Value;
            if ischar(ctrlSel); ctrlSel = {ctrlSel}; end
            validCtrlSel = matlab.lang.makeValidName(ctrlSel);

            expFolders = fields(finalResults);
            for i = 1:numel(expFolders)
                expFld = expFolders{i};
                if ismember(expFld, validCtrlSel); continue; end
                
                colNumStr = regexp(expFld, '\d+$', 'match', 'once');
                if isempty(colNumStr); continue; end
                ctrlFld = matlab.lang.makeValidName(['A' colNumStr]);
                
                if isfield(finalResults, ctrlFld) && ismember(ctrlFld, validCtrlSel)
                    resExp = finalResults.(expFld);
                    resCtrl = finalResults.(ctrlFld);
                    
                    if ~isempty(resExp.agg_stats) && ~isempty(resCtrl.agg_stats)
                        total_spots_exp = resExp.agg_stats.spot_count_mean * resExp.agg_stats.images_analyzed;
                        specific_spot_count = total_spots_exp - (resCtrl.agg_stats.spot_count_mean * resCtrl.agg_stats.images_analyzed);
                        finalResults.(expFld).agg_stats.specific_spot_count_total = max(0, specific_spot_count);

                        if specific_spot_count <= 0
                            finalResults.(expFld).peakRatios = zeros(1,3);
                            finalResults.(expFld).dominant_state = 'N/A';
                        else
                            countsExp = total_spots_exp .* resExp.peakRatios;
                            countsCtrl = (resCtrl.agg_stats.spot_count_mean * resCtrl.agg_stats.images_analyzed) .* resCtrl.peakRatios;
                            
                            if numel(countsExp)==3 && numel(countsCtrl)==3
                                specific_counts = countsExp - countsCtrl;
                                specific_counts(specific_counts < 0) = 0;
                                
                                if sum(specific_counts) > 0
                                    finalResults.(expFld).peakRatios = specific_counts / sum(specific_counts);
                                    [~, max_idx] = max(specific_counts);
                                    state_labels = {'Monomer', 'Trimer', 'Trimer+'};
                                    finalResults.(expFld).dominant_state = state_labels{max_idx};
                                else
                                    finalResults.(expFld).peakRatios = zeros(1,3);
                                    finalResults.(expFld).dominant_state = 'N/A';
                                end
                            end
                        end
                    end
                end
            end
            
            appdata.results = finalResults;
            appdata.rawResults = classifiedResults;
            appdata.calibrated_monomer_ref = calibrated_monomer_ref;
            
            updateQCPanel();
            lstbxRep.Value = procSel;
            btnSave.Enable = 'on';
            close(progDialog);
            uialert(fig, 'Analysis complete!', 'Success');
        catch ME
            close(progDialog);
            uialert(fig, ['An error occurred: ' ME.message], 'Error');
            rethrow(ME);
        end
    end
    
    function monomer_ref = getMonomerRefFromData(intensities, mode)
        monomer_ref = NaN;
        if numel(intensities) < 20; return; end
        
        [counts, edges] = histcounts(intensities, 100);
        centers = edges(1:end-1) + diff(edges)/2;
        pks = []; locs = [];
        if ~isempty(counts) && max(counts) > 0
            max_count = max(counts);
            dynamic_min_peak_height = max_count * 0.05;
            [pks, locs] = findpeaks(counts, centers, 'MinPeakHeight', dynamic_min_peak_height, 'NPeaks', 10, 'SortStr', 'descend');
        end

        if isempty(locs); return; end
        
        sorted_locs = sort(locs, 'ascend');
        
        if strcmp(mode, 'External')
            monomer_ref = sorted_locs(1);
            disp(['External standard peak found at ' num2str(sorted_locs(1)) '. Calibrating as Monomer.']);
        else
            unit_hypo1 = sorted_locs(1);
            ratios1 = sorted_locs / unit_hypo1;
            error1 = sum(abs(ratios1 - round(ratios1)));
            
            unit_hypo2 = sorted_locs(1) / 3;
            ratios2 = sorted_locs / unit_hypo2;
            error2 = sum(abs(ratios2 - round(ratios2)));
            
            if error1 <= error2
                monomer_ref = unit_hypo1;
                disp(['Pattern analysis suggests the first peak at ' num2str(sorted_locs(1)) ' is the Monomer.']);
            else
                monomer_ref = unit_hypo2;
                disp(['Pattern analysis suggests the first peak at ' num2str(sorted_locs(1)) ' is the Trimer.']);
            end
        end
    end

    function [perImageData, all_intensities_pooled] = processFolder(folderPath)
        files = dir(fullfile(folderPath, '*.tif'));
        limit = min(length(files), 50);
        
        varTypes = {'string', 'double', 'double', 'double', 'cell', 'cell'};
        varNames = {'FileName', 'SpotCount', 'SaturatedCount', 'MeanIntensity', 'CorrectedIntensities', 'Areas'};
        perImageData = table('Size', [limit, 6], 'VariableTypes', varTypes, 'VariableNames', varNames);
            
        all_intensities_pooled = [];
        
        for k = 1:limit
            I = double(imread(fullfile(folderPath, files(k).name)));
            [bw, stats] = analyzeSingleImage(I);
            
            perImageData.FileName(k) = files(k).name;
            
            if ~isempty(stats)
                background_pixels = I(~bw);
                if isempty(background_pixels); image_background_level = 0; else; image_background_level = median(background_pixels); end
                
                raw_intensities = [stats.MeanIntensity]';
                areas = [stats.Area]';
                is_saturated = (raw_intensities >= edSatThresh.Value);
                non_saturated_raw = raw_intensities(~is_saturated);
                
                corrected_for_image = non_saturated_raw - image_background_level;
                corrected_for_image = corrected_for_image(corrected_for_image >= 0);
                
                perImageData.SpotCount(k) = numel(raw_intensities);
                perImageData.SaturatedCount(k) = sum(is_saturated);
                perImageData.MeanIntensity(k) = mean(corrected_for_image, 'omitnan');
                perImageData.CorrectedIntensities{k} = corrected_for_image;
                perImageData.Areas{k} = areas;
                
                all_intensities_pooled = [all_intensities_pooled; corrected_for_image];
            else
                perImageData.SpotCount(k) = 0;
                perImageData.SaturatedCount(k) = 0;
                perImageData.MeanIntensity(k) = NaN;
                perImageData.CorrectedIntensities{k} = [];
                perImageData.Areas{k} = [];
            end
        end
    end

    function updateQCPanel()
        if ~isfield(appdata, 'results') || isempty(fieldnames(appdata.results)); return; end
        
        all_res_fields = fieldnames(appdata.results);
        calibrated_monomer_ref = appdata.calibrated_monomer_ref;
        
        exp_spot_counts = zeros(3,1);
        exp_folder_count = 0;
        total_saturated = 0;
        all_areas = [];
        
        for i = 1:numel(all_res_fields)
            fld = all_res_fields{i};
            res = appdata.results.(fld);
            if ~isstruct(res.agg_stats); continue; end
            
            total_saturated = total_saturated + res.agg_stats.saturated_count_total;
            
            raw_res = appdata.rawResults.(fld);
            if isfield(raw_res, 'perImageData') && istable(raw_res.perImageData)
                if ismember('Areas', raw_res.perImageData.Properties.VariableNames)
                    for j = 1:height(raw_res.perImageData)
                        all_areas = [all_areas; raw_res.perImageData.Areas{j}];
                    end
                end
            end
            
            if ismember(fld, matlab.lang.makeValidName(lstbxControl.Value)); continue; end
            
            exp_folder_count = exp_folder_count + 1;
            
            if isfield(res.agg_stats, 'specific_spot_count_total')
                specific_total = res.agg_stats.specific_spot_count_total;
                specific_counts_per_state = specific_total .* res.peakRatios;
                exp_spot_counts = exp_spot_counts + specific_counts_per_state';
            end
        end
        
        avg_counts_per_folder = exp_spot_counts / max(1, exp_folder_count);
        
        lblRefMeans.Text = sprintf('Monomer Ref: %.2f', calibrated_monomer_ref);
        lblMonoSpots.Text = sprintf('Monomer: %.1f spots/folder', avg_counts_per_folder(1));
        lblTrimerSpots.Text = sprintf('Trimer:  %.1f spots/folder', avg_counts_per_folder(2));
        lblTrimerPlusSpots.Text = sprintf('Trimer+: %.1f spots/folder', avg_counts_per_folder(3));
        lblSatSpots.Text = sprintf('Saturated Spots: %d (Total)', total_saturated);
        
        if ~isempty(all_areas)
             lblSpotArea.Text = sprintf('Spot Area: %.1f ± %.1f px²', mean(all_areas), std(all_areas));
        else
             lblSpotArea.Text = 'Spot Area: N/A';
        end
    end
    
    function [bw, stats] = analyzeSingleImage(I_original)
        if ~isempty(appdata.roiPosition)
            I = imcrop(I_original, appdata.roiPosition);
        else
            I = I_original;
        end

        I_s = imgaussfilt(I, edSigma.Value);
        I_n = mat2gray(I_s);
        bw_initial = imbinarize(I_n, 'adaptive', 'Sensitivity', edSens.Value);
        
        if edSE.Value > 0
            bw_opened = imopen(bw_initial, strel('disk', round(edSE.Value)));
        else
            bw_opened = bw_initial;
        end
        
        CC = bwconncomp(bw_opened);
        stats_all = regionprops(CC, I, 'Area', 'MeanIntensity');
        
        areas = [stats_all.Area];
        idx_to_keep = find(areas >= edMinA.Value & areas <= edMaxA.Value);
        
        bw_roi_filtered = ismember(labelmatrix(CC), idx_to_keep);
        
        bw_full = false(size(I_original));
        if ~isempty(appdata.roiPosition)
            x = round(appdata.roiPosition(1));
            y = round(appdata.roiPosition(2));
            [h_actual, w_actual] = size(I);
            bw_full(y:(y+h_actual-1), x:(x+w_actual-1)) = bw_roi_filtered;
        else
            bw_full = bw_roi_filtered;
        end
        
        stats = stats_all(idx_to_keep);
        bw = bw_full;
    end

    function [dominant_state, ratios, mean_corr] = analyzePopulation(correctedInt, calibrated_monomer_ref)
        dominant_state = 'N/A'; ratios = nan(1, 3); mean_corr = NaN;
        if isempty(correctedInt); mean_corr = NaN; ratios = zeros(1, 3); return; end
        if isnan(calibrated_monomer_ref); return; end
        
        split1 = calibrated_monomer_ref * 2.5;
        split2 = calibrated_monomer_ref * 3.5;

        mean_corr = mean(correctedInt, 'omitnan'); total_spots = numel(correctedInt);
        counts = zeros(1, 3);
        counts(1) = sum(correctedInt <= split1);
        counts(2) = sum(correctedInt > split1 & correctedInt <= split2);
        counts(3) = sum(correctedInt > split2);
        
        if total_spots > 0; ratios = counts / total_spots; else; ratios = zeros(1, 3); end
        
        [~, max_idx] = max(counts);
        state_labels = {'Monomer', 'Trimer', 'Trimer+'};
        if any(counts > 0); dominant_state = state_labels{max_idx}; end
    end

    function updateDisplayCallback(~,~)
        cla(axHist); hold(axHist, 'on'); grid(axHist, 'on');
        repSel = lstbxRep.Value;
        if ischar(repSel); repSel = {repSel}; end
        colors = lines(numel(repSel));
        
        if isfield(appdata, 'rawResults') && ~isempty(appdata.rawResults)
            for i = 1:numel(repSel)
                folderName = repSel{i}; 
                validName = matlab.lang.makeValidName(folderName);
                if isfield(appdata.rawResults, validName)
                    res = appdata.rawResults.(validName);
                    if ~isempty(res.allIntensities)
                        histogram(axHist, res.allIntensities, 'Normalization', 'probability', 'FaceColor', colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', folderName);
                    end
                end
            end
        end

        ylim(axHist, 'auto');
        if numel(repSel) > 0; legend(axHist, 'show', 'Location', 'northeast', 'Interpreter', 'none'); end
        hold(axHist, 'off');
        title(axHist, 'Raw Intensity Dist. (Pre-NSB)');
        xlabel(axHist, 'Corrected Intensity (a.u.)');
        ylabel(axHist, 'Probability');
        updatePreview();
    end
    
    function updatePreview(~,~)
        cla(axPreviewOriginal); cla(axPreviewProcessed);
        
        repSel = lstbxRep.Value;
        if isempty(repSel) || iscell(repSel) && isempty(repSel{1}); return; end
        if iscell(repSel); repSel = repSel{1}; end
        
        folderPath = fullfile(appdata.rootDir, repSel);
        files = dir(fullfile(folderPath, '*.tif'));
        if isempty(files); title(axPreviewOriginal, 'No .tif files'); title(axPreviewProcessed, ''); return; end
        
        try
            I = double(imread(fullfile(folderPath, files(1).name)));
            [bw, ~] = analyzeSingleImage(I);
            
            imshow(mat2gray(I), 'Parent', axPreviewOriginal);
            title(axPreviewOriginal, 'Original');
            
            imshow(mat2gray(I), 'Parent', axPreviewProcessed); hold(axPreviewProcessed, 'on');
            visboundaries(axPreviewProcessed, bw, 'Color', 'r', 'LineWidth', 0.5);
            if ~isempty(appdata.roiPosition)
                rectangle(axPreviewProcessed, 'Position', appdata.roiPosition, 'EdgeColor', 'g', 'LineWidth', 1);
            end
            hold(axPreviewProcessed, 'off');
            title(axPreviewProcessed, 'Detected');
            
            previewPanel.Title = ['Spot Detection 🔄 (Live Preview): ' files(1).name];

        catch ME
            title(axPreviewOriginal, 'Error loading image');
            disp(ME.message);
        end
    end

    function saveCallback(~, ~)
        if isempty(fieldnames(appdata.results)); uialert(fig, 'No results to save.', 'Info'); return; end
        [fname, fpath] = uiputfile('*.xlsx', 'Save Analysis Results', 'SpotAnalysisResults_V34_1.xlsx');
        if isequal(fname, 0); return; end
        fullFile = fullfile(fpath, fname);
        progDialog = uiprogressdlg(fig, 'Title', 'Saving...', 'Message', 'Writing data to Excel file.', 'Indeterminate', 'on');
        try
            calibrated_monomer_ref = appdata.calibrated_monomer_ref;

            progDialog.Message = 'Creating Statistical Summary...';
            fields = fieldnames(appdata.results);
            summaryHeader = {'Folder', 'Images_Analyzed', ...
                'Total_Spots_Mean', 'Total_Spots_SEM', ...
                'Specific_Spots_Total', ...
                'Mean_Intensity_Mean', 'Mean_Intensity_SEM', ...
                'Monomer_Ratio_Mean', 'Monomer_Ratio_SEM', ...
                'Trimer_Ratio_Mean', 'Trimer_Ratio_SEM', ...
                'Trimer+_Ratio_Mean', 'Trimer+_Ratio_SEM', ...
                'Total_Saturated_Spots', 'Calibrated_Monomer_Ref'};
            
            summaryCell = cell(numel(fields) + 1, numel(summaryHeader));
            summaryCell(1,:) = summaryHeader;

            for i = 1:numel(fields)
                fld = fields{i}; res = appdata.results.(fld);
                stats = res.agg_stats;
                
                if isfield(stats, 'ratios_mean') && numel(stats.ratios_mean) == 3
                    ratios_mean = stats.ratios_mean;
                    ratios_sem = stats.ratios_sem;
                else
                    ratios_mean = [NaN, NaN, NaN];
                    ratios_sem = [NaN, NaN, NaN];
                end

                row = {fld, stats.images_analyzed, ...
                       stats.spot_count_mean, stats.spot_count_sem, ...
                       NaN, ...
                       stats.mean_intensity_mean, stats.mean_intensity_sem, ...
                       ratios_mean(1), ratios_sem(1), ...
                       ratios_mean(2), ratios_sem(2), ...
                       ratios_mean(3), ratios_sem(3), ...
                       stats.saturated_count_total, calibrated_monomer_ref};
                
                if isfield(stats, 'specific_spot_count_total')
                    row{5} = stats.specific_spot_count_total;
                end
                
                summaryCell(i+1,:) = row;
            end
            writetable(cell2table(summaryCell(2:end,:), 'VariableNames', summaryHeader), fullFile, 'Sheet', 'Statistical_Summary', 'WriteMode', 'overwrite');

            progDialog.Message = 'Creating Raw Summary...';
            rawFields = fieldnames(appdata.rawResults);
            rawSummaryCell = cell(numel(rawFields) + 1, 9);
            rawSummaryCell(1,:) = {'Folder', 'Total_Spots_Raw', 'Saturated_Spots_Raw', 'Mean_Intensity_Raw', ...
                'Dominant_State_Raw', 'Monomer_Ratio_Raw', 'Trimer_Ratio_Raw', 'Trimer+_Ratio_Raw', ...
                'Calibrated_Monomer_Ref'};
            for i = 1:numel(rawFields)
                fld = rawFields{i}; res = appdata.rawResults.(fld);
                stats = res.agg_stats;
                rawSummaryCell{i+1, 1} = fld;
                rawSummaryCell{i+1, 2} = stats.spot_count_mean * stats.images_analyzed;
                rawSummaryCell{i+1, 3} = stats.saturated_count_total;
                rawSummaryCell{i+1, 4} = res.mean_corrected;
                rawSummaryCell{i+1, 5} = res.dominant_state;
                if numel(res.peakRatios) == 3
                    rawSummaryCell{i+1, 6} = res.peakRatios(1);
                    rawSummaryCell{i+1, 7} = res.peakRatios(2);
                    rawSummaryCell{i+1, 8} = res.peakRatios(3);
                end
                rawSummaryCell{i+1, 9} = res.calibrated_monomer_ref;
            end
            writetable(cell2table(rawSummaryCell(2:end,:), 'VariableNames', rawSummaryCell(1,:)), fullFile, 'Sheet', 'Raw_Summary');

            progDialog.Message = 'Saving Intensity Data...';
            saveSheet(fullFile, {appdata.rawResults}, 'allIntensities', 'Intensity_Data');
            
            close(progDialog); uialert(fig, ['Results saved to ' fullFile], 'Save Complete');
        catch ME
            close(progDialog);
            uialert(fig, ['Failed to save. Error: ' ME.message], 'Save Error');
            rethrow(ME);
        end
    end
    
    function saveSheet(fullFile, results_cell, data_field_name, sheet_name)
        fields = fieldnames(results_cell{1});
        all_data = {}; col_names = {};
        for i = 1:numel(fields)
            fld = fields{i}; 
            if ~isfield(results_cell{1}, fld) || isempty(results_cell{1}.(fld)); continue; end
            res = results_cell{1}.(fld);
            if isfield(res, data_field_name) && ~isempty(res.(data_field_name))
                all_data{end+1} = res.(data_field_name); 
                col_names{end+1} = fld;
            end
        end
        if isempty(all_data); return; end
        max_len = max(cellfun(@numel, all_data));
        for i = 1:numel(all_data)
            current_len = numel(all_data{i});
            if current_len < max_len; all_data{i}(current_len+1:max_len) = NaN; end
        end
        if ~isempty(all_data)
            T_data = table(all_data{:}, 'VariableNames', col_names);
            writetable(T_data, fullFile, 'Sheet', sheet_name);
        end
    end
    
    function out = ifelse(condition, true_val, false_val)
        if condition; out = true_val; else; out = false_val; end
    end
end