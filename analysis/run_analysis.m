function run_analysis()
% RUN_ANALYSIS
% Master script to run all analysis pipelines for audio, behavior, ephys,
% and system latency/movement.

close all;

%% SETUP ANALYSIS PATHS

% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'analysis', 'utils');
addpath(utils_dir);
root_dir = path_utils.repo_root();

% Add script file paths
addpath(fullfile(root_dir, 'analysis', 'audio'));
addpath(fullfile(root_dir, 'analysis', 'behavior'));
addpath(fullfile(root_dir, 'analysis', 'ephys'));
addpath(fullfile(root_dir, 'analysis', 'system'));

% Clean out existing data in the results
path_utils.clean_results_dirs();

%% AUDIO ANALYSIS
gantry_audio_analysis;
gate_audio_analysis;

%% BEHAVIOR ANALYSIS
behavior_all_analysis;
behavior_phase3_analysis;
behavior_trajectory_analysis;

%% EPHYS ANALYSIS
gantry_ephys_analysis;
gate_ephys_analysis;

%% SYSTEM-LEVEL ANALYSIS
gantry_movement_analysis;
system_latency_analysis;

clc; close all;
end
