function gate_audio_analysis()
% =========================================================================
% FUNCTION: gate_audio_analysis
%
% Description:
%   Analyzes gate movement audio + ROS bag data.
%   All movements are treated as a single class ('move') regardless of direction.
%   Computes per-event stats, generates plots, and prints summary.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'gate_audio_test';
dataDir = fullfile(path_utils.dataset_root(), out_dir_name);
figDir = path_utils.results_figures(out_dir_name);
calibPath = fullfile(path_utils.dataset_root(), 'audio_calibration', 'calibration_params.mat');

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% --- USER SETTINGS ---
bagFiName   = 'ros_session_20250424_190556.bag';
audioFiName = 'audio_recording_20250424_190556.wav';

audioBandFilter      = [200, 90000];   % Hz
sessionStartOffset   = 2.0;            % Seconds before first movement event
windowDuration       = 1.0;            % Event window duration (s)
baselineStartOffset  = 3.0;            % Baseline starts this many seconds after movement
plotTimeRange        = [0, 25];        % X-axis range for time series plot

plotCol = [0.5 0.5 0.5;   % Baseline
    0.2 0.6 1.0];  % Move (unified)

posF1 = [100, 500, 800, 300];
posF2 = [100, 100, 250, 300];
posF3 = [400, 100, 450, 300];

% --- PIPELINE EXECUTION ---
[gateEvents] = load_ros_events(dataDir, bagFiName);
[moveEvents, moveTimes, timeZero] = extract_gate_movement_events(gateEvents);
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
log_utils.printf('[DONE] Results and plots saved to disk.\n');

log_utils.close()
end % gate_audio_analysis

% -------------------------------------------------------------------------
function [gateEvents] = load_ros_events(dataDir, bagFilename)
% Load ROS bag and extract /gate_testing messages
bag_path = fullfile(dataDir, bagFilename);
log_utils.printf('[INFO] Loading ROS bag file: %s\n', bagFilename);
bag = rosbag(bag_path);
gateBag = select(bag, 'Topic', '/gate_testing');
msgs = readMessages(gateBag, 'DataFormat', 'struct');

gateEvents(length(msgs)) = struct('Time', [], 'Event', '');
for i = 1:length(msgs)
    gateEvents(i).Time  = gateBag.MessageList.Time(i);
    gateEvents(i).Event = msgs{i}.Data;
end
end

% -------------------------------------------------------------------------
function [moveEvents, moveTimes, timeZero] = extract_gate_movement_events(gateEvents)
% Extract all movement events (up/down), treat them as 'move'

moveEvents = {};
moveTimes  = [];
timeZero   = [];

for i = 1:length(gateEvents)
    ev = gateEvents(i).Event;
    if contains(ev, 'gate_test_wall_move_up') || contains(ev, 'gate_test_wall_move_down')
        moveEvents{end+1} = 'move';
        moveTimes(end+1)  = gateEvents(i).Time;

        if isempty(timeZero) && contains(ev, 'gate_test_wall_move_up')
            timeZero = gateEvents(i).Time;
        end
    end
end

if isempty(timeZero)
    error('No "gate_test_wall_move_up" event found to anchor timeZero.');
end
end

% -------------------------------------------------------------------------
function [audioData, audioFs] = load_audio(dataDir, audioFilename)
% Load .wav audio file
audio_path = fullfile(dataDir, audioFilename);
log_utils.printf('[INFO] Loading audio file: %s\n', audioFilename);
[audioData, audioFs] = audioread(audio_path);
log_utils.printf('[INFO] Loaded %d samples @ %.1f Hz\n', length(audioData), audioFs);
end

function audioFiltered = bandpass_filter_audio(audioData, audioBandFilter, audioFs)
% Applies a bandpass filter to the raw audio signal
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', audioBandFilter(1), ...
    'HalfPowerFrequency2', audioBandFilter(2), ...
    'SampleRate', audioFs);

audioFiltered = filtfilt(bpFilt, audioData);
log_utils.printf('[INFO] Applied bandpass filter: %.0f–%.0f Hz\n', ...
    audioBandFilter(1), audioBandFilter(2));
end

% -------------------------------------------------------------------------
function calibrationGain_G = load_calibration(calibrationFilePath)
% Loads audio calibration parameters from the specified .mat file
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
function audioPa = apply_audio_calibration(audioData, calibrationGain_G)
% Converts raw audio samples to Pascals
ref_Pa   = 20e-6;
audioPa  = audioData .* calibrationGain_G .* ref_Pa;
end

% -------------------------------------------------------------------------
function audioDB = compute_audio_dB_envelope(audioPa, audioFs)
% Computes dB SPL envelope from calibrated audio in Pascals
rmsWindow_ms = 20;
ref_Pa = 20e-6;
winSamples = round(audioFs * rmsWindow_ms / 1000);
rmsEnvelope = sqrt(movmean(audioPa.^2, winSamples));
audioDB = 20 * log10(rmsEnvelope / ref_Pa);
end

% -------------------------------------------------------------------------
function [dataWindowCorrected, audioDBCorrected, baselineMean, splBaselineFloor] = compute_baseline_correction(dataWindow, audioDB, windowDuration, audioFs)
% Applies baseline correction to amplitude stats and audio trace

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
function dataWindow = analyze_windows(audioData, audioFs, moveEvents, moveTimes, timeZero, sessionStartOffset, windowDuration, baselineStartOffset)
% Compute amplitude stats for each window

dataWindow = struct('Type', {}, 'Time', {}, 'SampleIdx', {}, ...
    'MeanAmplitude', {}, 'MinAmplitude', {}, 'MaxAmplitude', {}, 'StdAmplitude', {});

for i = 1:length(moveEvents)
    rel_time = (moveTimes(i) - timeZero) + sessionStartOffset;
    start_sample = round(rel_time * audioFs) + 1;
    end_sample   = start_sample + round(windowDuration * audioFs) - 1;

    if start_sample < 1 || end_sample > length(audioData)
        log_utils.printf('[WARN] Skipping event %d: out-of-bounds\n', i);
        continue;
    end

    window = audioData(start_sample:end_sample);
    dataWindow(end+1).Type        = 'move';
    dataWindow(end).Time          = moveTimes(i);
    dataWindow(end).SampleIdx     = start_sample;
    dataWindow(end).MeanAmplitude = mean(window);
    dataWindow(end).MinAmplitude  = min(window);
    dataWindow(end).MaxAmplitude  = max(window);
    dataWindow(end).StdAmplitude  = std(window);

    baseline_start_sample = start_sample + round(baselineStartOffset * audioFs);
    baseline_end_sample   = baseline_start_sample + round(windowDuration * audioFs) - 1;

    if baseline_end_sample > length(audioData)
        log_utils.printf('[WARN] Skipping baseline for event %d: out-of-bounds\n', i);
        continue;
    end

    baseline_window = audioData(baseline_start_sample:baseline_end_sample);
    dataWindow(end+1).Type        = 'baseline';
    dataWindow(end).Time          = moveTimes(i) + baselineStartOffset;
    dataWindow(end).SampleIdx     = baseline_start_sample;
    dataWindow(end).MeanAmplitude = mean(baseline_window);
    dataWindow(end).MinAmplitude  = min(baseline_window);
    dataWindow(end).MaxAmplitude  = max(baseline_window);
    dataWindow(end).StdAmplitude  = std(baseline_window);
end
end

% -------------------------------------------------------------------------
function plot_time_series_with_windows(posF1, plotCol, audioDB, audioFs, dataWindow, windowDuration, outputDir, plotTimeRange)
% Plot full audio time series with transparent patches for event windows

% --- Settings ---
patchAlpha = 0.2;             
ylimVals   = [-10, 50];       

% Time axis
t = (0:length(audioDB)-1) / audioFs;

% Create figure
figure;
plot(t, audioDB, 'k'); hold on;
set(gcf, 'Position', posF1);

% Draw transparent patches
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

% Axis formatting
set(gca, 'FontName', 'Arial', 'FontSize', 15);
xlabel('Time (s)');
ylabel('Δ SPL (dB rel. baseline)');
ylim(ylimVals);
box off

% X limits
if ~isempty(plotTimeRange)
    xlim(plotTimeRange);
end

% Save plot
filename = fullfile(outputDir, 'gate_audio_time_series.svg');
print(filename, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved time-series waveform plot with event windows: %s\n', filename);
end

% -------------------------------------------------------------------------
function plot_amplitude_bar_strip(posF2, plotCol, dataWindow, outputDir)
types = {'baseline', 'move'};
amps = cellfun(@(t) [dataWindow(strcmp({dataWindow.Type}, t)).MeanAmplitude], ...
    types, 'UniformOutput', false);
means = cellfun(@mean, amps);
stds  = cellfun(@std, amps);
x_pos = 1:numel(types);

figure; hold on;
set(gcf, 'Position', posF2);

barWidth = 0.6;
for i = 1:numel(types)
    bar(x_pos(i), means(i), barWidth, 'FaceColor', plotCol(i,:));
end
errorbar(x_pos, means, stds, 'k.', 'CapSize', 6, 'LineWidth', 1.5);

set(gca, 'XTick', x_pos, 'XTickLabel', {'Baseline', 'Move'});
set(gca, 'FontName', 'Arial', 'FontSize', 15)
ylabel('Mean SPL (dB SPL)');
xlabel('Condition');
ylim([0, 125]);
xlim([0.5, numel(types) + 0.5]);

filename = fullfile(outputDir, 'gate_audio_bar.svg');
print(filename, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved bar+strip plot: %s\n', filename);
end

% -------------------------------------------------------------------------
function plot_spectral_comparison(posF3, pltCol, winMeta, audioPa, Fs, outputDir, band)
ref_Pa = 20e-6;
fftLen = 8192;
win    = hann(4096);
olap   = numel(win)/2;

types = {'baseline','move'};

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
        'Color', pltCol(k,:), 'LineWidth', 1.6);
end

xlim(band);
ylim([0, 100]);
xlabel('Frequency (Hz, log scale)');
ylabel('PSD (dB re (20 µPa)^2/Hz)');
set(gca, 'XTick', [20000, 40000, 80000], ...
    'XTickLabel', {'20k','40k','80k'}, ...
    'FontName', 'Arial', 'FontSize', 15);
box off;

filename = fullfile(outputDir, 'gate_audio_spectrum.svg');
print(filename, '-dsvg', '-r300');
end

% -------------------------------------------------------------------------
function print_stats(dataWindow, dataWindowCorrected, audioFs, baselineMean, splBaselineFloor, calibrationGain_G, audioBandFilter)
types = {'baseline', 'move'};
fields = {'MeanAmplitude', 'MinAmplitude', 'MaxAmplitude'};
field_labels = {'Mean Amplitude', 'Min Amplitude', 'Max Amplitude'};

log_utils.printf('\n=========== Audio Analysis Parameters ===========\n');
log_utils.printf('Sampling Rate (Fs):           %d Hz\n', audioFs);
log_utils.printf('Baseline Mean SPL:            %.2f dB SPL\n', baselineMean);
log_utils.printf('SPL Floor (Baseline Periods): %.2f dB SPL\n', splBaselineFloor);
log_utils.printf('Calibration Gain (G):         %.6f\n', calibrationGain_G);
log_utils.printf('Audio Band Filter:            [%d, %d] Hz\n', audioBandFilter(1), audioBandFilter(2));

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