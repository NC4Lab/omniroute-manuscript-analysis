function behavior_phase3_analysis()
% =========================================================================
% FUNCTION: behavior_phase3_analysis
%
% Description:
%   Stand-alone analysis of Phase 3 trials for two rats (23, 24). For each
%   rat, include only Phase 3 sessions that contain trials from all four
%   start chambers (ChamberNumber ∈ {1,3,5,7}). Compute per-session success
%   rates by start, then aggregate across sessions to produce a bar plot
%   (mean ± SD across sessions) styled to match success_rate_bar_plot.
%
%   Uses only trial_times.csv in each session Data_Dir. No fallbacks.
%
% Outputs:
%   - One SVG bar plot per rat: success rate by start (N, E, S, W), mean ± SD.
%   - Console logs with basic stats similar in style to other scripts.
%
% Conventions:
%   ChamberNumber → Start index → Label → Color
%      1 → 1 → N (Top)    → orange  [1.00 0.65 0.00]
%      3 → 2 → E (Right)  → teal    [0.00 0.65 0.65]
%      5 → 3 → S (Bottom) → green   [0.35 0.70 0.30]
%      7 → 4 → W (Left)   → purple  [0.60 0.40 0.80]
%
% Save dirs:
%   Uses the same directories as behavior_performance_analysis.m for
%   session summary input and figure outputs.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'behavior_phase3_test';
figDir = path_utils.results_figures(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% ------------------------- USER SETTINGS ---------------------------------
rats = [23, 24]; % rat numbers

% Session summary data
sesTabInPath = fullfile(path_utils.dataset_root(), 'processed_behavior', 'session_summary.csv');

% Bar plot formatting (carry over where applicable)
posF = [650, 500, 225, 300];    % compact bar figure, similar scale to prior plots
tickLabels = {'N','E','S','W'};
yLim = [0 100];
fontName = 'Arial';
fontSize = 15;

% Colors (N, E, S, W) in the established order
plotCol = [ ...
    1.00 0.65 0.00;  % N (Top)    orange
    0.00 0.65 0.65;  % E (Right)  teal
    0.35 0.70 0.30;  % S (Bottom) green
    0.60 0.40 0.80]; % W (Left)   purple

% Chamber mapping and order
validChambers = [1 3 5 7];   % must all be present to include a session
% Map ChamberNumber -> start index (1..4)
%   1->1 (N), 3->2 (E), 5->3 (S), 7->4 (W)
chamberToStartIdx = containers.Map({1,3,5,7}, {1,2,3,4});

% --------------------------- LOAD SUMMARY --------------------------------
if ~isfile(sesTabInPath)
    error('Missing session_summary.csv at: %s', sesTabInPath);
end
S = readtable(sesTabInPath, 'VariableNamingRule','preserve');

reqCols = {'Rat','Phase','Session_Number','Performance_All','Data_Dir'};
missing = setdiff(reqCols, S.Properties.VariableNames);
if ~isempty(missing)
    error('session_summary.csv missing required column(s): %s', strjoin(missing, ', '));
end

% ---------------------------- PER RAT ------------------------------------
for ri = 1:numel(rats)
    ratID = rats(ri);

    % Filter to Phase 3 rows for this rat
    isRat   = toNum(S.Rat) == ratID;
    isP3    = toNum(S.Phase) == 3;
    keepIdx = isRat & isP3;
    if ~any(keepIdx)
        log_utils.printf('[WARN] No Phase 3 entries in session_summary for Rat %d\n', ratID);
        continue;
    end

    dirsP3 = S.Data_Dir(keepIdx);
    sessNumsP3 = toNum(S.Session_Number(keepIdx));

    % Sort by session number for reproducible order
    [sessNumsP3, ord] = sort(sessNumsP3(:));
    dirsP3 = dirsP3(ord);

    % Containers for included sessions (must have all 4 start chambers)
    ratesPerSession = [];   % rows = sessions, cols = 4 starts (N,E,S,W), values in %
    trialsPerSession = [];  % same shape, trial counts per start
    includedSessNums = [];
    includedSessDirs = {};

    % Also track total trials and successes per start for reporting
    totalTrialsByStart  = zeros(1,4);
    totalSuccessByStart = zeros(1,4);

    % Loop sessions
    for s = 1:numel(dirsP3)
        sessPath = strtrim(asChar(dirsP3{s}));
        trialFile = fullfile(sessPath, 'trial_times.csv');

        if ~isfile(trialFile)
            log_utils.printf('[WARN][RAT %02d] Missing trial_times.csv in %s (skipping session)\n', ratID, sessPath);
            continue;
        end

        T = readtable(trialFile, 'VariableNamingRule','preserve');

        % Validate required columns
        needCols = {'Result','ChamberNumber'};
        if ~all(ismember(needCols, T.Properties.VariableNames))
            log_utils.printf('[WARN][RAT %02d] %s missing required columns (%s). Skipping.\n', ...
                ratID, sessPath, strjoin(setdiff(needCols, T.Properties.VariableNames), ', '));
            continue;
        end

        % Filter valid trials: Result in {SUCCESS, ERROR} and ChamberNumber ∈ {1,3,5,7}
        resultStr = upper(string(T.Result));
        chamberVal = toNum(T.ChamberNumber);
        isValidRes = resultStr == "SUCCESS" | resultStr == "ERROR";
        isValidCh  = ismember(chamberVal, validChambers);
        keep = isValidRes & isValidCh;

        if ~any(keep)
            log_utils.printf('[WARN][RAT %02d] No valid trials after filtering in %s\n', ratID, sessPath);
            continue;
        end

        resultStr = resultStr(keep);
        chamberVal = chamberVal(keep);

        % Count per start for this session
        sessCounts  = zeros(1,4);
        sessCorrect = zeros(1,4);

        for k = 1:numel(chamberVal)
            ch = chamberVal(k);
            st = chamberToStartIdx(ch); % 1..4
            sessCounts(st) = sessCounts(st) + 1;
            if resultStr(k) == "SUCCESS"
                sessCorrect(st) = sessCorrect(st) + 1;
            end
        end

        % Require all 4 starts present in this session
        if any(sessCounts == 0)
            continue;
        end

        % Session-level success (%) per start
        sessRates = 100 * (sessCorrect ./ max(1, sessCounts));

        % Accumulate included session
        ratesPerSession = [ratesPerSession; sessRates];      %#ok<AGROW>
        trialsPerSession = [trialsPerSession; sessCounts];   %#ok<AGROW>
        includedSessNums(end+1) = sessNumsP3(s);             %#ok<AGROW>
        includedSessDirs{end+1} = sessPath;                  %#ok<AGROW>

        % Totals for reporting
        totalTrialsByStart  = totalTrialsByStart  + sessCounts;
        totalSuccessByStart = totalSuccessByStart + sessCorrect;
    end

    % Aggregate across included sessions
    nIncluded = size(ratesPerSession, 1);
    nTotalP3  = numel(dirsP3);

    if nIncluded == 0
        log_utils.printf('[WARN] Rat %02d: No Phase 3 sessions with all four starts. No plot produced.\n', ratID);
        continue;
    end

    mu = mean(ratesPerSession, 1, 'omitnan');     % mean across sessions, %
    sd = std(ratesPerSession, 0, 1, 'omitnan');   % SD across sessions
    nSessPerStart = sum(~isnan(ratesPerSession), 1); % should all equal nIncluded

    % ---------------------------- PLOT -----------------------------------
    F = figure('Color','w', 'Position', posF); ax = axes('Parent',F); hold(ax,'on'); box(ax,'off');
    set(ax, 'FontName', fontName, 'FontSize', fontSize);

    x_pos = 1:4; barWidth = 0.6;

    for j = 1:4
        if ~isnan(mu(j))
            bar(ax, x_pos(j), mu(j), barWidth, 'FaceColor', plotCol(j,:));
        else
            bar(ax, x_pos(j), 0, barWidth, 'FaceColor', plotCol(j,:), ...
                'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
        end
    end

    % SD error bars
    errorbar(ax, x_pos, mu, sd, 'k.', 'CapSize', 6, 'LineWidth', 1.5);

    % Axes formatting
    set(ax, 'XTick', x_pos, 'XTickLabel', tickLabels);
    ylim(ax, yLim);
    xlim(ax, [0.5, 4.5]);
    yline(ax, 50, 'k--', 'LineWidth', 1.2);
    ylabel(ax, 'Success rate (%)');
    xlabel(ax, 'Start condition');

%     % Annotate n = sessions above bars
%     for j = 1:4
%         if ~isnan(mu(j))
%             text(ax, x_pos(j), mu(j) + 3, sprintf('n=%d', nSessPerStart(j)), ...
%                 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11);
%         end
%     end

    % Save
    outName = fullfile(figDir, sprintf('phase3_success_by_start_avg_RAT%02d.svg', ratID));
    print(F, outName, '-dsvg', '-r300');
    log_utils.printf('[INFO] Saved Phase 3 success-by-start bar (mean ± SD) for Rat %02d: %s\n', ratID, outName);

    % -------------------------- LOG STATS --------------------------------
    stats = struct();
    stats.ratID = ratID;
    stats.nPhase3Total = nTotalP3;
    stats.nIncludedAllFour = nIncluded;
    stats.includedSessionNumbers = includedSessNums(:)';
    stats.includedSessionDirs = includedSessDirs;
    stats.mu_by_start = mu;
    stats.sd_by_start = sd;
    stats.nSess_by_start = nSessPerStart;
    stats.totalTrials_by_start = totalTrialsByStart;
    stats.totalSuccess_by_start = totalSuccessByStart;
    stats.labels = tickLabels;

    print_phase3_start_stats(stats);
end

log_utils.printf('[DONE] Phase 3 success-by-start (per rat) complete.\n');

log_utils.close()
end % behavior_phase3_analysis

% ========================================================================
% Helpers
% ========================================================================
function num = toNum(x)
% Robustly convert table column to numeric
if isnumeric(x)
    num = x;
elseif iscell(x)
    num = cellfun(@(v) str2double(asChar(v)), x);
elseif isstring(x)
    num = str2double(x);
elseif iscategorical(x)
    num = double(x);
else
    num = str2double(string(x));
end
end

function s = asChar(v)
% Convert value to char row vector safely
if ischar(v)
    s = v;
elseif isstring(v)
    s = char(v);
else
    s = char(string(v));
end
end

function print_phase3_start_stats(stats)
% Console summary similar in spirit to print_behavior_stats

labels = stats.labels;
mu  = stats.mu_by_start;
sd  = stats.sd_by_start;
ns  = stats.nSess_by_start;
Ntr = stats.totalTrials_by_start;
Ktr = stats.totalSuccess_by_start;

log_utils.printf('\n================== PHASE 3 START-SIDE SUMMARY (RAT %02d) ==================\n', stats.ratID);
log_utils.printf('Phase 3 sessions in summary: %d\n', stats.nPhase3Total);
log_utils.printf('Included sessions with all 4 starts: %d\n', stats.nIncludedAllFour);

if ~isempty(stats.includedSessionNumbers)
    log_utils.printf('Included Session Numbers: %s\n', sprintf('%d ', stats.includedSessionNumbers));
else
    log_utils.printf('Included Session Numbers: none\n');
end

log_utils.printf('\nPer-start aggregates (mean ± SD across sessions; n = sessions contributing that start)\n');
for i = 1:numel(labels)
    lbl = labels{i};
    log_utils.printf('  %s: %6.2f%% ± %5.2f  (n = %d sessions)  | total trials = %5d, total correct = %5d\n', ...
        lbl, mu(i), sd(i), ns(i), Ntr(i), Ktr(i));
end
log_utils.printf('Reference: 50%% line drawn in figure.\n');
log_utils.printf('==========================================================================\n\n');
end
