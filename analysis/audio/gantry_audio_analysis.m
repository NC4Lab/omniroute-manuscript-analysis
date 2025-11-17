function gantry_audio_analysis()
% =========================================================================
% FUNCTION: gantry_audio_analysis
%
% Description:
%   Modular script to analyze gantry movement audio + ROS bag data.
%   Computes per-move audio stats, treating all movements identically.
%   Baseline and movement windows are plotted and analyzed uniformly.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'gantry_audio_test';
dataDir = fullfile(path_utils.dataset_root(), out_dir_name);
figDir = path_utils.results_figures(out_dir_name);
calibPath = fullfile(path_utils.dataset_root(), 'audio_calibration', 'calibration_params.mat');

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% --- USER SETTINGS ---
bagFiName   = 'ros_session_20250729_161226.bag';
audioFiName = 'audio_recording_20250729_161226.wav';

audioBandFilter      = [200, 90000];  % Hz
sessionStartOffset   = 2.0;           % Audio started this many seconds before first move
windowDuration       = 1.5;           % Analysis window duration (s)
baselineStartOffset  = 2.5;           % Baseline window starts this many seconds after move
plotTimeRange        = [5, 25];      % Limit on plot x-axis (set [] for full trace)

plotCol = [0.5 0.5 0.5;    % Baseline
    1.0 0.3 0.3];   % Move

posF1 = [100, 500, 800, 300];
posF2 = [100, 100, 250, 300];
posF3 = [400, 100, 450, 300];

% ---------------------- PIPELINE EXECUTION -------------------------------
[gantryEvents] = load_ros_events(dataDir, bagFiName);
[moveEvents, moveTimes, timeZero] = extract_gantry_movement_events(gantryEvents);
[audioData, audioFs] = load_audio(dataDir, audioFiName);
audioData = bandpass_filter_audio(audioData, audioBandFilter, audioFs);

calibrationGain_G = load_calibration(calibPath);
audioPa = apply_audio_calibration(audioData, calibrationGain_G);
audioDB = compute_audio_dB_envelope(audioPa, audioFs);

dataWindow = analyze_windows(audioDB, audioFs, moveEvents, moveTimes, timeZero, sessionStartOffset, windowDuration, baselineStartOffset);
[dataWindowCorrected, audioDBCorrected, baselineMean, splBaselineFloor] = compute_baseline_correction(dataWindow, audioDB, windowDuration, audioFs);

plot_time_series_with_windows(posF1, plotCol, audioDBCorrected, audioFs, dataWindowCorrected, windowDuration, figDir, plotTimeRange);
plot_amplitude_bar_strip(posF2, plotCol, dataWindow, figDir);
plot_spectral_comparison(posF3, plotCol, dataWindow, audioPa, audioFs, figDir, audioBandFilter);

print_stats(dataWindow, dataWindowCorrected, audioFs, baselineMean, splBaselineFloor, calibrationGain_G, audioBandFilter);
log_utils.printf('[DONE] Gantry audio analysis complete. Results saved.\n');

log_utils.close()
end % gantry_audio_analysis

% -------------------------------------------------------------------------
function [gantryEvents] = load_ros_events(dataDir, bagFiName)
bag_path = fullfile(dataDir, bagFiName);
log_utils.printf('[INFO] Loading ROS bag file: %s\n', bagFiName);
bag = rosbag(bag_path);
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

% 1) Collect all move events and mark the very first as timeZero
for i = 1:length(gantryEvents)
    ev = gantryEvents(i).Event;
    if startsWith(ev, 'gantry_test_move_')
        if isempty(timeZero)
            timeZero = gantryEvents(i).Time;  % keep t=0 behavior for plots
        end
        moveEvents{end+1} = 'move';           %#ok<AGROW>
        moveTimes(end+1)  = gantryEvents(i).Time; %#ok<AGROW>
    end
end

if isempty(timeZero)
    error('No gantry_test_move_* events found to anchor timeZero.');
end
if numel(moveTimes) < 2
    error('Fewer than 2 move events found. Cannot drop the first.');
end

% 2) Drop the first move entirely; keep the next 20 (2..21) if available
keep_last_idx = min(21, numel(moveTimes));  % index 21 => 20 kept after dropping #1
moveEvents = moveEvents(2:keep_last_idx);
moveTimes  = moveTimes(2:keep_last_idx);

log_utils.printf('[INFO] Extracted %d gantry movement events.\n', numel(moveEvents));
end

% -------------------------------------------------------------------------
function [audioData, audioFs] = load_audio(dataDir, audioFiName)
audio_path = fullfile(dataDir, audioFiName);
log_utils.printf('[INFO] Loading audio file: %s\n', audioFiName);
[audioData, audioFs] = audioread(audio_path);
log_utils.printf('[INFO] Loaded %d samples @ %.1f Hz\n', length(audioData), audioFs);
end

% -------------------------------------------------------------------------
function audioFiltered = bandpass_filter_audio(audioData, band, Fs)
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', band(1), ...
    'HalfPowerFrequency2', band(2), ...
    'SampleRate', Fs);
audioFiltered = filtfilt(bpFilt, audioData);
log_utils.printf('[INFO] Applied bandpass filter: %.0f–%.0f Hz\n', band(1), band(2));
end

% -------------------------------------------------------------------------
function calibrationGain_G = load_calibration(calibrationFilePath)
if ~isfile(calibrationFilePath)
    error('Calibration file not found: %s', calibrationFilePath);
end
data = load(calibrationFilePath);
if ~isfield(data, 'calibrationGain_G')
    error('Missing field "calibrationGain_G" in calibration file.');
end
calibrationGain_G = data.calibrationGain_G;
log_utils.printf('[INFO] Loaded calibration gain: G = %.6f\n', calibrationGain_G);
end

% -------------------------------------------------------------------------
function audioPa = apply_audio_calibration(audioData, G)
ref_Pa = 20e-6;
audioPa = audioData .* G .* ref_Pa;
end

% -------------------------------------------------------------------------
function audioDB = compute_audio_dB_envelope(audioPa, Fs)
rmsWindow_ms = 20;
ref_Pa = 20e-6;
winSamples = round(Fs * rmsWindow_ms / 1000);
rmsEnvelope = sqrt(movmean(audioPa.^2, winSamples));
audioDB = 20 * log10(rmsEnvelope / ref_Pa);
end

% -------------------------------------------------------------------------
function dataWindow = analyze_windows(audioData, Fs, moveEvents, moveTimes, timeZero, sessionStartOffset, winDur, baseOffset)
dataWindow = struct('Type', {}, 'Time', {}, 'SampleIdx', {}, ...
    'MeanAmplitude', {}, 'MinAmplitude', {}, 'MaxAmplitude', {}, 'StdAmplitude', {});

for i = 1:length(moveEvents)
    rel_time = (moveTimes(i) - timeZero) + sessionStartOffset;
    start_sample = round(rel_time * Fs) + 1;
    end_sample   = start_sample + round(winDur * Fs) - 1;

    if start_sample < 1 || end_sample > length(audioData)
        log_utils.printf('[WARN] Skipping event %d: out-of-bounds\n', i);
        continue;
    end

    window = audioData(start_sample:end_sample);

    dataWindow(end+1).Type          = 'move';
    dataWindow(end).Time            = moveTimes(i);
    dataWindow(end).SampleIdx       = start_sample;
    dataWindow(end).MeanAmplitude   = mean(window);
    dataWindow(end).MinAmplitude    = min(window);
    dataWindow(end).MaxAmplitude    = max(window);
    dataWindow(end).StdAmplitude    = std(window);

    baseline_start_sample = start_sample + round(baseOffset * Fs);
    baseline_end_sample   = baseline_start_sample + round(winDur * Fs) - 1;

    if baseline_end_sample > length(audioData)
        log_utils.printf('[WARN] Skipping baseline for event %d: out-of-bounds\n', i);
        continue;
    end

    baseline_window = audioData(baseline_start_sample:baseline_end_sample);

    dataWindow(end+1).Type          = 'baseline';
    dataWindow(end).Time            = moveTimes(i) + baseOffset;
    dataWindow(end).SampleIdx       = baseline_start_sample;
    dataWindow(end).MeanAmplitude   = mean(baseline_window);
    dataWindow(end).MinAmplitude    = min(baseline_window);
    dataWindow(end).MaxAmplitude    = max(baseline_window);
end
end

% -------------------------------------------------------------------------
function [dataWindowCorrected, audioDBCorrected, baselineMean, splBaselineFloor] = compute_baseline_correction(dataWindow, audioDB, windowDuration, audioFs)
baselineAmps = [dataWindow(strcmp({dataWindow.Type}, 'baseline')).MeanAmplitude];
baselineMean = mean(baselineAmps);
dataWindowCorrected = dataWindow;

for i = 1:length(dataWindow)
    dataWindowCorrected(i).MeanAmplitude = dataWindow(i).MeanAmplitude - baselineMean;
    dataWindowCorrected(i).MinAmplitude  = dataWindow(i).MinAmplitude  - baselineMean;
    dataWindowCorrected(i).MaxAmplitude  = dataWindow(i).MaxAmplitude  - baselineMean;
end

audioDBCorrected = audioDB - baselineMean;

baselineSamples = [];
for i = 1:length(dataWindow)
    if strcmp(dataWindow(i).Type, 'baseline')
        idx = dataWindow(i).SampleIdx;
        duration = round(windowDuration * audioFs);  % 1s window
        segment = audioDB(idx:(idx + duration - 1));
        baselineSamples = [baselineSamples; segment];
    end
end

splBaselineFloor = min(baselineSamples);
end

% -------------------------------------------------------------------------
function plot_time_series_with_windows(posF1, plotCol, audioDB, audioFs, dataWindow, windowDuration, figDir, plotTimeRange)
patchAlpha = 0.2;
ylimVals   = [-10, 50];

t = (0:length(audioDB)-1) / audioFs;

figure;
plot(t, audioDB, 'k'); hold on;
set(gcf, 'Position', posF1);

for i = 1:length(dataWindow)
    t_start = dataWindow(i).SampleIdx / audioFs;
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
xlabel('Time (s)');
ylabel('Δ SPL (dB rel. baseline)');
ylim(ylimVals);
box off

if ~isempty(plotTimeRange)
    xlim(plotTimeRange);
    ticks  = plotTimeRange(1):5:plotTimeRange(2);   % 5:5:25
    labels = ticks - plotTimeRange(1);              % 0:5:20
    set(gca, 'XTick', ticks, 'XTickLabel', string(labels));
end

FiName = fullfile(figDir, 'gantry_audio_time_series.svg');
print(FiName, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved time-series waveform plot with event windows: %s\n', FiName);
end

% -------------------------------------------------------------------------
function plot_amplitude_bar_strip(posF2, plotCol, dataWindow, figDir)
types = {'baseline', 'move'};
x_pos = 1:numel(types);

amps = cellfun(@(t) [dataWindow(strcmp({dataWindow.Type}, t)).MeanAmplitude], ...
    types, 'UniformOutput', false);
means = cellfun(@mean, amps);
stds  = cellfun(@std, amps);

figure; hold on;
set(gcf, 'Position', posF2);

barWidth = 0.6;

for i = 1:numel(types)
    bar(x_pos(i), means(i), barWidth, 'FaceColor', plotCol(i,:));
end

errorbar(x_pos, means, stds, 'k.', 'CapSize', 6, 'LineWidth', 1.5);

set(gca, 'XTick', x_pos, 'XTickLabel', {'Baseline', 'Move'});
set(gca, 'FontName', 'Arial', 'FontSize', 15);
ylabel('Mean SPL (dB SPL)');
xlabel('Condition');
ylim([0, 125]);
xlim([0.5, numel(types) + 0.5]);

FiName = fullfile(figDir, 'gantry_audio_bar.svg');
print(FiName, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved bar+strip plot: %s\n', FiName);
end

% -------------------------------------------------------------------------
function plot_spectral_comparison(posF3, plotCol, winMeta, audioPa, Fs, figDir, band)
ref_Pa = 20e-6;
fftLen = 8192;
win    = hann(4096);
olap   = numel(win)/2;

types = {'baseline', 'move'};
leg   = {'Baseline', 'Move'};

figure; hold on;
set(gcf,'Position', posF3);

for k = 1:numel(types)
    idx = strcmp({winMeta.Type}, types{k});
    allPSD = [];
    for j = find(idx)
        s0 = winMeta(j).SampleIdx;
        s1 = s0 + round(Fs) - 1;
        if s1 > length(audioPa), continue; end
        seg = audioPa(s0:s1);
        [pxx, f] = pwelch(seg, win, olap, fftLen, Fs);
        allPSD = [allPSD; pxx'];
    end
    if isempty(allPSD), continue; end
    pxx_dB = 10*log10(mean(allPSD,1) / ref_Pa^2);
    semilogx(f, smoothdata(pxx_dB,'movmean',5), ...
        'Color', plotCol(k,:), 'LineWidth', 1.6);
end

xlim(band);
ylim([0, 100]);
xlabel('Frequency (Hz, log scale)');
ylabel('PSD (dB re (20 µPa)^2/Hz)');

xticks = [20000, 40000, 80000];
xticklabels = arrayfun(@(x) sprintf('%dk', x/1000), xticks, 'UniformOutput', false);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

set(gca, 'FontName', 'Arial', 'FontSize', 15);
box off;

FiName = fullfile(figDir, 'gantry_audio_spectrum.svg');
print(FiName, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved spectral comparison plot: %s\n', FiName);
end

% -------------------------------------------------------------------------
function print_stats(dataWindow, dataWindowCorrected, Fs, baselineMean, splBaselineFloor, G, band)
types = {'baseline', 'move'};
fields = {'MeanAmplitude', 'MinAmplitude', 'MaxAmplitude'};
field_labels = {'Mean Amplitude', 'Min Amplitude', 'Max Amplitude'};

log_utils.printf('\n=========== Audio Analysis Parameters ===========\n');
log_utils.printf('Sampling Rate (Fs):           %d Hz\n', Fs);
log_utils.printf('Baseline Mean SPL:            %.2f dB SPL\n', baselineMean);
log_utils.printf('SPL Floor (Baseline Periods): %.2f dB SPL\n', splBaselineFloor);
log_utils.printf('Calibration Gain (G):         %.6f\n', G);
log_utils.printf('Audio Band Filter:            [%d, %d] Hz\n', band(1), band(2));

log_utils.printf('\n=========== Absolute Audio Summary ===========\n\n');
for i = 1:numel(types)
    count = sum(strcmp({dataWindow.Type}, types{i}));
    log_utils.printf('  %s   : %d events\n', upper(types{i}), count);
end
log_utils.printf('\n');
for f = 1:numel(fields)
    log_utils.printf('%s:\n', field_labels{f});
    for i = 1:numel(types)
        vals = [dataWindow(strcmp({dataWindow.Type}, types{i})).(fields{f})];
        log_utils.printf('  %s   : %0.2f ± %0.2f dB SPL\n', ...
            upper(types{i}), mean(vals), std(vals));
    end
    log_utils.printf('\n');
end

log_utils.printf('=========== Baseline-Corrected ΔL Summary ===========\n\n');
for f = 1:numel(fields)
    log_utils.printf('%s:\n', field_labels{f});
    for i = 1:numel(types)
        vals = [dataWindowCorrected(strcmp({dataWindowCorrected.Type}, types{i})).(fields{f})];
        log_utils.printf('  %s   : %+0.2f ± %0.2f dB\n', ...
            upper(types{i}), mean(vals), std(vals));
    end
    log_utils.printf('\n');
end
end
