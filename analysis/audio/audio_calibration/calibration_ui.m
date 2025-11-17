function calibration_ui()
% =========================================================================
% FUNCTION: calibration_ui
%
% Description:
%   Interactive UI for manual scoring of tone calibration segments.
%   Allows user to click tone start/end points and associate them with
%   known dB SPL values (94, 104, 114). Computes RMS for each segment and
%   checks for intra-tone consistency. Saves scored intervals and data
%   for use in downstream calibration.
% =========================================================================

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'audio_calibration';
dataDir = fullfile(path_utils.dataset_root(), out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% --- USER SETTINGS ---
audioFiName = 'calibration_audio_20250417_134132.wav';
uiOutFiName    = 'calibration_ui_output.mat';
knownLevels   = [94, 104, 114];  % Expected calibration tone levels (dB SPL)

% --- LOAD AUDIO ---
[audioData, fs] = audioread(fullfile(dataDir, audioFiName));
t = (0:length(audioData)-1) / fs;

% --- INTERNAL STATE ---
clickBuffer = [];
storedWindows = [];
storedLevels  = [];
rawRMS        = [];
rawdB         = [];

% --- SETUP UI ---
fig = figure('Name', 'Calibration Tone Selector', ...
    'Position', [100, 100, 1400, 600], ...
    'NumberTitle', 'off', ...
    'WindowButtonDownFcn', @click_callback);

movegui(fig, 'center');

ax = axes('Parent', fig, ...
    'Units', 'normalized', ...
    'Position', [0.05, 0.25, 0.9, 0.7]);
plot(ax, t, audioData, 'k');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Amplitude');
title(ax, 'Click tone start/end: Green = Start, Red = End');
grid on;

% --- UI CONTROLS ---
uicontrol('Style', 'text', 'String', 'Select Tone Level', ...
    'Units', 'normalized', 'Position', [0.05, 0.15, 0.2, 0.05], ...
    'FontSize', 11, 'HorizontalAlignment', 'left');

btnGroup = uibuttongroup('Visible','on',...
    'Units','normalized',...
    'Position',[0.05, 0.05, 0.4, 0.1], ...
    'SelectionChangedFcn', @tone_select_callback);

for i = 1:length(knownLevels)
    uicontrol(btnGroup, 'Style','togglebutton',...
        'String', sprintf('%d dB', knownLevels(i)), ...
        'Units','normalized',...
        'Position', [0.01 + 0.33*(i-1), 0.05, 0.3, 0.9], ...
        'Tag', num2str(knownLevels(i)));
end

btnGroup.SelectedObject = [];

storeBtn = uicontrol('Style','pushbutton','String','Store Level',...
    'Units','normalized','Position',[0.7, 0.08, 0.12, 0.06],...
    'Callback', @store_level_callback, 'FontSize', 11);

doneBtn = uicontrol('Style','pushbutton','String','Done',...
    'Units','normalized','Position',[0.85, 0.08, 0.1, 0.06],...
    'Callback', @done_callback, 'FontSize', 11);

% --- TRACK ACTIVE STATE ---
currentLevel = [];

% ------------------------- CALLBACKS -------------------------------------

    function tone_select_callback(~, event)
        currentLevel = str2double(event.NewValue.Tag);
        clickBuffer = [];  % Reset buffer on level switch
        log_utils.printf('[INFO] Selected level: %d dB\n', currentLevel);
    end

    function click_callback(~, ~)
        if isempty(currentLevel)
            log_utils.printf('[WARN] Select a tone level first.\n');
            return;
        end

        cp = get(ax, 'CurrentPoint');
        t_click = cp(1,1);
        clickBuffer(end+1) = t_click;

        if mod(length(clickBuffer), 2) == 1
            xline(ax, t_click, 'g', 'LineWidth', 1.5);
        else
            xline(ax, t_click, 'r', 'LineWidth', 1.5);
        end
    end

    function store_level_callback(~, ~)
        if length(clickBuffer) ~= 8
            log_utils.printf('[WARN] You must click 4 start/end tone pairs (8 clicks total).\n');
            return;
        end
        ranges = reshape(clickBuffer, 2, [])';

        log_utils.printf('\n[INFO] Scoring tone: %d dB\n', currentLevel);
        segRMS = zeros(4, 1);
        for i = 1:4
            idx_start = max(1, round(ranges(i,1) * fs));
            idx_end   = min(length(audioData), round(ranges(i,2) * fs));
            seg = audioData(idx_start:idx_end);
            segRMS(i) = rms(seg);
            log_utils.printf('  Segment %d: RMS = %.6f\n', i, segRMS(i));
        end

        % Consistency check
        dB_vals = 20 * log10(segRMS);
        dB_range = max(dB_vals) - min(dB_vals);
        dB_std   = std(dB_vals);
        log_utils.printf('  dB Range: %.2f dB   |  Std Dev: %.2f dB\n', dB_range, dB_std);
        if dB_range > 1.5
            log_utils.printf('  [WARNING] High variability detected — check segment placement.\n');
        else
            log_utils.printf('  [OK] Segment RMS levels are consistent.\n');
        end

        % Store data
        storedWindows = [storedWindows; ranges];
        storedLevels  = [storedLevels; repmat(currentLevel, 4, 1)];
        rawRMS        = [rawRMS; segRMS];
        rawdB         = [rawdB; dB_vals];

        clickBuffer = [];

        % Delete all red and green xlines
        lines = findall(ax, 'Type', 'ConstantLine');
        for i = 1:numel(lines)
            if isequal(lines(i).Color, [1 0 0]) || isequal(lines(i).Color, [0 1 0])
                delete(lines(i));
            end
        end
    end

    function done_callback(~, ~)
        close(fig);
        if isempty(storedWindows)
            error('No tone windows were stored.');
        end

        % Save data to file
        save(fullfile(dataDir, uiOutFiName), ...
            'storedWindows', 'storedLevels', 'rawRMS', 'rawdB');

        log_utils.printf('[DONE] Manual scoring saved to: %s\n', fullfile(dataDir, uiOutFiName));
    end

log_utils.close()
end
