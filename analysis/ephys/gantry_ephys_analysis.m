function gantry_ephys_analysis()
% =========================================================================
% FUNCTION: gantry_ephys_analysis
%
% Description:
%   Analyze electrophysiological responses to gantry movement events.
%   Adapts gate_ephys_analysis to use the same event extraction and timing
%   offsets as gantry_audio_analysis.
%   Compares all movement events against baseline without splitting by direction.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'gantry_ephys_test';
figDir = path_utils.results_figures(out_dir_name);
tabDir = path_utils.results_tables(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% --- USER SETTINGS ---

% Path to MATLAB .mat file with CSC data and ROS .bag file
matPath = fullfile(path_utils.dataset_root(), 'NC40023\20250806_170156\processed\matlab_csc\20250806_170156_csc_export.mat');
bagPath = fullfile(path_utils.dataset_root(), 'NC40023\20250806_170156\raw\ROS\ros_session_20250806_170207.bag');

% Timing parameters (matching gantry_audio_analysis)
windowDuration      = 1.5;   % analysis window duration (s)
baselineStartOffset = 2.5;   % seconds after move for baseline window
plotTimeRange       = [3, 23];

plotCol = [0.5 0.5 0.5;   % Baseline
           1.0 0.3 0.3];  % Move (unified)

channelsExclude = [3, 4, 6, 28, 30, 31, 60, 62];

posF1 = [100, 500, 800, 300];
posF2 = [100, 100, 250, 300];
posF3 = [400, 100, 450, 300];

% --- LOAD DATA ---
mat = load(matPath);
cscData = mat.csc_data * mat.gain_to_uv;
cscRosTS = mat.ros_ts;
cscFs = mat.sampling_rate;

% --- EXCLUDE CHANNELS ---
if ~isempty(channelsExclude)
    cscData(:, channelsExclude) = [];  % remove columns for excluded channels
end
nChIncluded = size(cscData, 2);

% --- PRINT VARIABLE VALUES ---
log_utils.printf('gain_to_uv: %.6f\n', mat.gain_to_uv);
log_utils.printf('mean(mat.csc_data): %.4f (raw units)\n', mean(mean(mat.csc_data)));
log_utils.printf('mean(cscData): %.4f µV\n', mean(mean(cscData)));
log_utils.printf('Channels included: %d (excluded: %s)\n', nChIncluded, mat2str(channelsExclude));

% --- ANALYSIS ---
gantryEvents = load_ros_events(bagPath);
[~, moveTimes, timeZero] = extract_gantry_movement_events(gantryEvents);
dataWindow = analyze_windows(cscData, moveTimes, windowDuration, cscFs, baselineStartOffset, cscRosTS);

% --- PLOTS ---
plot_time_series_with_windows_flat(posF1, plotCol, cscData, cscRosTS, dataWindow, windowDuration, figDir, plotTimeRange, timeZero);
plot_amplitude_bar_strip_flat(posF2, plotCol, dataWindow, figDir);
plot_spectral_comparison_flat(posF3, plotCol, cscData, dataWindow, windowDuration, cscFs, figDir);

% --- STATS ---
compute_ttest(dataWindow);
print_stats_flat(dataWindow);

% --- SAVE ---
save(fullfile(tabDir, 'gantry_test_data_flat.mat'), 'gantryEvents', 'cscRosTS', 'cscFs', 'dataWindow');
log_utils.printf('[DONE] Results and plots saved to disk.\n');

log_utils.close()
end % gantry_ephys_analysis

% ------------------- HELPER FUNCTIONS -------------------

function gantryEvents = load_ros_events(bagPath)
bag = rosbag(bagPath);
gantryBag = select(bag, 'Topic', '/gantry_testing');
msgs = readMessages(gantryBag, 'DataFormat', 'struct');
gantryEvents(length(msgs)) = struct('Time', [], 'Event', '');
for i = 1:length(msgs)
    gantryEvents(i).Time  = gantryBag.MessageList.Time(i);
    gantryEvents(i).Event = msgs{i}.Data;
end
end

% -------------------------------------------------------------------------
function [moveEvents, moveTimes, timeZero] = extract_gantry_movement_events(gantryEvents)
% Collect all gantry_test_move_* events, anchor timeZero to the first one,
% then drop that first move from analysis and keep the next 20.

moveEvents = {};
moveTimes  = [];
timeZero   = [];

% Collect all move events, remember the first as timeZero
for i = 1:numel(gantryEvents)
    ev = gantryEvents(i).Event;
    if startsWith(ev, 'gantry_test_move_')
        if isempty(timeZero)
            timeZero = gantryEvents(i).Time;   % keep t=0 behavior
        end
        moveEvents{end+1} = 'move';            %#ok<AGROW>
        moveTimes(end+1)  = gantryEvents(i).Time; %#ok<AGROW>
    end
end

if isempty(timeZero)
    error('No gantry_test_move_* events found to anchor timeZero.');
end
if numel(moveTimes) < 2
    error('Fewer than 2 move events. Cannot drop the first.');
end

% Drop the first, keep up to 20 following (i.e., indices 2..21)
keep_last_idx = min(21, numel(moveTimes));
moveEvents = moveEvents(2:keep_last_idx);
moveTimes  = moveTimes(2:keep_last_idx);

log_utils.printf('[INFO] Using %d move events (dropped the first; timeZero anchored at first move).\n', numel(moveEvents));
end

% -------------------------------------------------------------------------
function dataWindow = analyze_windows(cscData, moveTimes, windowDuration, cscFs, baselineStartOffset, cscRosTS)
% =========================================================================
% FUNCTION: analyze_windows
%
% Description:
%   Extract and compute stats for windows anchored directly to ROS move
%   event times (no sessionStartOffset). For each move event:
%     - Move window starts at moveTimes(i) (event onset).
%     - Baseline window starts at moveTimes(i) + baselineStartOffset.
%
% RMS policy:
%   Per-window, per-channel de-mean (AC RMS), then aggregate across channels.
%
% Outputs per window (one struct per window):
%   Type          : 'move' or 'baseline'
%   Time          : ROS timestamp at window start (seconds)
%   SampleIdx     : starting sample index into cscData for this window
%   MeanAmplitude : mean of per-channel RMS (µV)
%   MinAmplitude  : min across channels of per-channel RMS (µV)
%   MaxAmplitude  : max across channels of per-channel RMS (µV)
%   StdAmplitude  : std across channels of per-channel RMS (µV)
% =========================================================================

% Preliminaries
nSamples = size(cscData, 1);
winLen   = max(1, round(windowDuration * cscFs));

% Output container
dataWindow = struct('Type', {}, 'Time', {}, 'SampleIdx', {}, ...
    'MeanAmplitude', {}, 'MinAmplitude', {}, 'MaxAmplitude', {}, 'StdAmplitude', {});

% Iterate over move events
for i = 1:numel(moveTimes)
    % ---------------- Movement window ----------------
    move_start_time = moveTimes(i);  % anchor directly to event onset
    move_start_idx  = find(cscRosTS >= move_start_time, 1, 'first');

    if isempty(move_start_idx) || (move_start_idx + winLen - 1 > nSamples)
        log_utils.printf('[WARN] Skipping MOVE %d: start sample missing or window exceeds data bounds.\n', i);
    else
        seg = cscData(move_start_idx : move_start_idx + winLen - 1, :);

        % AC RMS per channel
        seg0    = seg - mean(seg, 1, 'omitnan');      % remove per-channel DC
        rms_ch  = sqrt(mean(seg0.^2, 1, 'omitnan'));  % 1 x nCh
        rms_avg = mean(rms_ch, 'omitnan');

        dataWindow(end+1) = struct( ... %#ok<AGROW>
            'Type',          'move', ...
            'Time',          move_start_time, ...
            'SampleIdx',     move_start_idx, ...
            'MeanAmplitude', rms_avg, ...
            'MinAmplitude',  min(rms_ch), ...
            'MaxAmplitude',  max(rms_ch), ...
            'StdAmplitude',  std(rms_ch, 0, 'omitnan'));
    end

    % ---------------- Baseline window ----------------
    base_start_time = moveTimes(i) + baselineStartOffset;  % anchor directly to offset after event
    base_start_idx  = find(cscRosTS >= base_start_time, 1, 'first');

    if isempty(base_start_idx) || (base_start_idx + winLen - 1 > nSamples)
        log_utils.printf('[WARN] Skipping BASELINE %d: start sample missing or window exceeds data bounds.\n', i);
    else
        seg = cscData(base_start_idx : base_start_idx + winLen - 1, :);

        % AC RMS per channel
        seg0    = seg - mean(seg, 1, 'omitnan');
        rms_ch  = sqrt(mean(seg0.^2, 1, 'omitnan'));
        rms_avg = mean(rms_ch, 'omitnan');

        dataWindow(end+1) = struct( ... %#ok<AGROW>
            'Type',          'baseline', ...
            'Time',          base_start_time, ...
            'SampleIdx',     base_start_idx, ...
            'MeanAmplitude', rms_avg, ...
            'MinAmplitude',  min(rms_ch), ...
            'MaxAmplitude',  max(rms_ch), ...
            'StdAmplitude',  std(rms_ch, 0, 'omitnan'));
    end
end
end

% -------------------------------------------------------------------------
function dataGroup = compute_group_stats_flat(dataWindow)
types = {'move'};
fields = {'MeanAmplitude','MinAmplitude','MaxAmplitude','StdAmplitude'};
for i = 1:numel(types)
    group_data = dataWindow(strcmp({dataWindow.Type}, types{i}));
    for j = 1:numel(fields)
        vals = [group_data.(fields{j})];
        dataGroup.(types{i}).(fields{j}).mean = mean(vals);
        dataGroup.(types{i}).(fields{j}).std  = std(vals);
        dataGroup.(types{i}).(fields{j}).min  = min(vals);
        dataGroup.(types{i}).(fields{j}).max  = max(vals);
    end
end
end

% -------------------------------------------------------------------------
function plot_time_series_with_windows_flat(posF1, plotCol, cscData, cscRosTS, dataWindow, windowDuration, outputDir, plotTimeRange, timeZero)
patchAlpha = 0.2;
ylimVals = [-1000 1000];
t = cscRosTS - timeZero;
avg_csc = mean(cscData, 2);
figure;
set(gcf, 'Position', posF1);
plot(t, avg_csc, 'k'); hold on;

for i = 1:length(dataWindow)
    t_start = cscRosTS(dataWindow(i).SampleIdx) - timeZero;
    t_end   = t_start + windowDuration;
    if strcmpi(dataWindow(i).Type, 'baseline')
        col = plotCol(1,:);
    else
        col = plotCol(2,:);
    end
    patch([t_start t_end t_end t_start], ...
        [ylimVals(1) ylimVals(1) ylimVals(2) ylimVals(2)], ...
        col, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
end

set(gca, 'FontName', 'Arial', 'FontSize', 15);
xlabel('Time (s)'); ylabel('Amplitude (µV)');
ylim(ylimVals); box off
if ~isempty(plotTimeRange)
    xlim(plotTimeRange);
    ax = gca;
    ticks  = plotTimeRange(1):5:plotTimeRange(2);
    labels = ticks - plotTimeRange(1);  % e.g., 0:5:20 for [5 25]
    set(ax, 'XTick', ticks, 'XTickLabel', string(labels));
end
print(fullfile(outputDir, 'gantry_ephys_time_series.svg'), '-dsvg', '-r300');
end

% -------------------------------------------------------------------------
function plot_amplitude_bar_strip_flat(posF2, plotCol, dataWindow, outputDir)
types = {'baseline', 'move'};
amps = cellfun(@(t) [dataWindow(strcmp({dataWindow.Type}, t)).MeanAmplitude], ...
    types, 'UniformOutput', false);
means = cellfun(@mean, amps);
stds  = cellfun(@std, amps);
x_pos = 1:numel(types);

figure; hold on;
set(gcf, 'Position', posF2);
for i = 1:numel(types)
    bar(x_pos(i), means(i), 0.6, 'FaceColor', plotCol(i,:));
end
errorbar(x_pos, means, stds, 'k.', 'CapSize', 6, 'LineWidth', 1.5);
set(gca, 'XTick', x_pos, 'XTickLabel', {'Baseline', 'Move'}, 'FontName', 'Arial', 'FontSize', 15);
ylabel('RMS amplitude (µV)');
xlabel('Condition');
xlim([0.5, 2.5]); 
ylim([0, 300]); 
box off
print(fullfile(outputDir, 'gantry_ephys_bar.svg'), '-dsvg', '-r300');
end

% -------------------------------------------------------------------------
function plot_spectral_comparison_flat(posF3, plotCol, cscData, dataWindow, windowDuration, cscFs, outputDir)
types = {'baseline','move'};
patchAlpha = 0.1;
ylimVals = [-10, 50];
freqRange = [1, 500];
fftLen = 8192;
win = hann(4096);
olap = numel(win)/2;

figure; hold on;
set(gcf, 'Position', posF3);
for k = 1:2
    idx = strcmp({dataWindow.Type}, types{k});
    allPSD = [];

    for j = find(idx)
        s0 = dataWindow(j).SampleIdx;
        s1 = s0 + round(windowDuration * cscFs) - 1;
        if s1 > size(cscData, 1)
            continue;
        end
        seg = cscData(s0:s1, :);
        for ch = 1:size(seg,2)
            [pxx, f] = pwelch(seg(:,ch), win, olap, fftLen, cscFs);
            allPSD = [allPSD; pxx'];  
        end
    end

    if isempty(allPSD)
        warning('No PSD data found for type %s', types{k});
        continue;
    end

    f_vec = f;
    freqMask = (f_vec >= freqRange(1)) & (f_vec <= freqRange(2));
    f_plot = f_vec(freqMask);

    mean_lin = mean(allPSD(:, freqMask), 1);
    std_lin  = std(allPSD(:, freqMask), [], 1);

    meanPSD = 10 * log10(mean_lin);
    semPSD  = 10 * log10(mean_lin + std_lin) - meanPSD;

    f_plot = f_plot(:);
    meanPSD = meanPSD(:);
    semPSD = semPSD(:);

    fill([f_plot; flipud(f_plot)], ...
        [meanPSD + semPSD; flipud(meanPSD - semPSD)], ...
        plotCol(k,:), 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');

    plot(f_plot', meanPSD', 'Color', plotCol(k,:), 'LineWidth', 1.6);
end

xlim(freqRange); ylim(ylimVals);
xlabel('Frequency (Hz)'); ylabel('PSD (dB re 1 µV^2/Hz)');
set(gca, 'FontName', 'Arial', 'FontSize', 15);
box off;
print(fullfile(outputDir, 'gantry_ephys_spectrum.svg'), '-dsvg', '-r300');
end

% -------------------------------------------------------------------------
function compute_ttest(dataWindow)
% COMPUTE_TTEST
% Paired two-tailed t test comparing baseline vs move mean amplitudes.
% Reports: t(df), p value, Cohen's d, 95% CI of the mean difference.

baseline = [dataWindow(strcmp({dataWindow.Type}, 'baseline')).MeanAmplitude];
move     = [dataWindow(strcmp({dataWindow.Type}, 'move')).MeanAmplitude];

% Ensure paired lengths
n = min(length(baseline), length(move));
baseline_paired = baseline(1:n);
move_paired     = move(1:n);

% Paired t test
[~, p, ci, stats] = ttest(baseline_paired, move_paired, 'Tail', 'both');

% Effect size (Cohen’s d for paired samples)
diffs = move_paired - baseline_paired;
d = mean(diffs) / std(diffs);

% Summary stats (means ± SD)
m_base = mean(baseline_paired);
sd_base = std(baseline_paired);
m_move = mean(move_paired);
sd_move = std(move_paired);

% Print results
log_utils.printf('\n=========== T-Test Summary ===========\n');
log_utils.printf('Two-tailed paired t test comparing Baseline vs Move amplitudes:\n');
log_utils.printf('Baseline = %.2f ± %.2f µV, Move = %.2f ± %.2f µV (mean ± SD, n = %d)\n', ...
    m_base, sd_base, m_move, sd_move, n);
log_utils.printf('t(%d) = %.2f, p = %.3f, Cohen''s d = %.2f, 95%% CI [%.2f, %.2f]\n', ...
    stats.df, stats.tstat, p, d, ci(1), ci(2));
end


% -------------------------------------------------------------------------
function print_stats_flat(dataWindow)
types = {'baseline','move'};
fields = {'MeanAmplitude','MinAmplitude','MaxAmplitude'};
log_utils.printf('\n=========== Summary Stats ===========\n\n');
for f = 1:numel(fields)
    log_utils.printf('%s:\n', fields{f});
    for t = 1:numel(types)
        vals = [dataWindow(strcmp({dataWindow.Type}, types{t})).(fields{f})];
        log_utils.printf('  %s : %0.2f ± %0.2f µV\n', upper(types{t}), mean(vals), std(vals));
    end
    log_utils.printf('\n');
end
end
