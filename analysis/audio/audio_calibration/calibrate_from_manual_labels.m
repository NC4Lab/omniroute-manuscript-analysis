function calibrate_from_manual_labels()
% =========================================================================
% FUNCTION: calibrate_from_manual_labels
%
% Description:
%   Loads manually labeled tone segments and known SPL levels.
%   Computes calibration parameters to convert raw RMS audio values
%   into calibrated dB SPL using:
%       dB_SPL = 20*log10(RMS * G)  ← Linear gain version
%       or
%       dB_SPL = 20*log10(RMS) + C  ← Additive offset version
%
%   Saves both C and G to the calibration file.
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
uiInFiName = 'calibration_ui_output.mat';
calOutFiName = 'calibration_params.mat';

% --- LOAD LABELED DATA ---
in = load(fullfile(dataDir, uiInFiName));
storedLevels = in.storedLevels;   % Nx1 known SPL values (e.g., [94; 94; ...])
rawdB        = in.rawdB;          % Nx1 raw dB = 20*log10(RMS)
rawRMS       = in.rawRMS;

if length(storedLevels) ~= length(rawdB)
    error('Mismatch between number of tone levels and dB values.');
end

% --- COMPUTE OFFSET AND GAIN ---
calibrationConstant_C = mean(storedLevels - rawdB);
calibrationGain_G = 10^(calibrationConstant_C / 20);  % New: linear multiplier

% --- PRINT RESULTS ---
log_utils.printf('\n[RESULTS] Calibration Complete\n');
log_utils.printf('  Calibration Offset C = %+0.4f dB\n', calibrationConstant_C);
log_utils.printf('  Linear Gain G        = %.6f\n', calibrationGain_G);
log_utils.printf('  Applied as: dB_SPL = 20*log10(RMS * G)\n');
log_utils.printf('  Raw dB Range: %0.2f to %0.2f\n', min(rawdB), max(rawdB));
log_utils.printf('  SPL Target Range: %0.2f to %0.2f\n', min(storedLevels), max(storedLevels));

% --- SAVE BOTH VERSIONS ---
calibrationDate  = datestr(now);
calibratedDevice = 'Pettersson M500-384';
calibrationTone  = 'Sinusoidal (94/104/114 dB)';

save(fullfile(dataDir, calOutFiName), ...
    'calibrationConstant_C', ...
    'calibrationGain_G', ...
    'storedLevels', 'rawRMS', 'rawdB', ...
    'calibrationDate', 'calibratedDevice', 'calibrationTone');

log_utils.printf('[DONE] Saved calibration parameters to: %s\n', ...
    fullfile(dataDir, calOutFiName));

log_utils.close()
end
