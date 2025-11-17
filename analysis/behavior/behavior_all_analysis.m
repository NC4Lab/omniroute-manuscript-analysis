function behavior_all_analysis()
% =========================================================================
% FUNCTION: behavior_all_analysis
%
% Description:
%   Analyze behavioral performance across ALL training phases (1–3) for
%   each rat. Uses session_summary.csv to determine session→phase mapping.
%
%   Produces, for EACH RAT:
%     (1) Session-level success plot:
%         - % correct per session plotted in black
%         - phase membership indicated by full-height colored patches
%     (2) Trial-level success plots:
%         - Phase 1 and Phase 3 figures using boxcharts, where each x=trial
%           index shows the distribution of per-session success rates (%)
%           computed across that session’s blocks
%         - Dashed mean line is the across-session mean (%) at each trial
%           index (session-centric, not block-centric)
%     (3) Phase 1 vs Phase 3 bar plot of session-level % correct
%
%   Notes:
%     - Session-level performance is taken from Performance_All in the
%       session_summary. If values are in 0–1, they are converted to %;
%       if already 0–100, they are left unchanged.
%     - Trial-level performance curves are recomputed from trial_times.csv
%       Result values. For plotting we aggregate per-session trialwise
%       success rates across that session’s blocks (continuous 0–100%),
%       which are suitable for boxplots.
%     - Backward compatibility: block-level matrices (.M) and aggregates
%       (.mu/.sem/.n) are retained, but plotting now consumes the new
%       per-session fields (.M_sess, .mu_sess/.sem_sess/.n_sess, .X, .Y).
%     - Results are cached to avoid repeated parsing.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'behavior_all_test';
figDir = path_utils.results_figures(out_dir_name);
tabDir = path_utils.results_tables(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% --- USER SETTINGS --------------------------------------------------------
rats = [23, 24]; % rat numbers

% Session summary data
sesTabInPath = fullfile(path_utils.dataset_root(), 'processed_behavior', 'session_summary.csv');

% Phase colors (P1,P2,P3)
plotCol = [ ...
    0.20 0.60 1.00;  % Blue
    1.00 0.80 0.20;  % Yellow
    1.00 0.20 0.20]; % Red

% Figure positions and axis limits
posF1 = [100, 100, 650, 300];  % Session-level figure (per rat)
posF2 = [100, 500, 200, 300];  % Trial-level (P1) figure (per rat)
posF3 = [350, 500, 200, 300];  % Trial-level (P3) figure (per rat)
posF4 = [650, 500, 125, 300]; % Bar plot figure (per rat)
posF5 = [900, 500, 125, 300]; % Bar plot figure (per rat)
plotXLimF1 = [];               % e.g., [0, 40]; empty => auto
plotXLimF2 = [1,4];            % e.g., [1, 60]; empty => auto
plotXLimF3 = [1,4];            % e.g., [1, 60]; empty => auto

% Cache behavior
doRecompute = true;

% ---------------------- PIPELINE EXECUTION -------------------------------
% Load session summary (phase membership and per-session performance)
summaryTbl = readtable(sesTabInPath, 'VariableNamingRule','preserve');

[sessionPerf, trialPerf] = load_or_process_data( ...
    summaryTbl, rats, tabDir, doRecompute);

% Plots: per-rat figures
if ~isfolder(figDir), mkdir(figDir); end
for i = 1:numel(rats)
    ratID = rats(i);

    % Session-level plots
    plot_success_by_session(sessionPerf(i), ratID, plotCol, posF1, plotXLimF1, figDir);

    % Phase 1: compute per-session trialwise success (%) and plot
    p1 = analyze_phase1_block_aligned(sessionPerf(i));
    if ~isempty(p1.M) || ~isempty(p1.M_sess)
        plot_success_by_trial(p1, ratID, 1, plotCol, posF2, plotXLimF2, figDir);
    end

    % Phase 3: compute per-session trialwise success (%) and plot
    p3 = analyze_phase3_block_aligned(sessionPerf(i));
    if ~isempty(p3.M) || ~isempty(p3.M_sess)
        plot_success_by_trial(p3, ratID, 3, plotCol, posF3, plotXLimF3, figDir);
    end

    % Phase 1: trial 1 vs later bar plot
    if ~isempty(p1.M)
        plot_trial1_vs_later_bar(p1, ratID, 1, posF4, plotCol, figDir);
    end

    % Phase 3: trial 1 vs later bar plot
    if ~isempty(p3.M)
        plot_trial1_vs_later_bar(p3, ratID, 3, posF5, plotCol, figDir);
    end
end

% Print summary
print_behavior_stats(sessionPerf, trialPerf, rats);

log_utils.printf('[DONE] Behavior performance (all phases) complete. Results saved.\n');

log_utils.close()
end % behavior_all_analysis

% -------------------------------------------------------------------------
function [sessionPerf, trialPerf] = load_or_process_data(summaryTbl, rats, ioDir, doRecompute)
cacheFile = fullfile(ioDir, 'behavior_performance_cache_allphases.mat');
if isfile(cacheFile) && ~doRecompute
    log_utils.printf('[INFO] Loading cached performance data: %s\n', cacheFile);
    S = load(cacheFile, 'sessionPerf', 'trialPerf');
    sessionPerf = S.sessionPerf;
    trialPerf   = S.trialPerf;
    return
end

log_utils.printf('[INFO] Processing session_summary and trial_times (no cache or recompute=true)...\n');

% Ensure output dir
if ~isfolder(ioDir), mkdir(ioDir); end

% Expected columns used (names must match exactly):
% 'Rat','Phase','Session_Number','Performance_All','Data_Dir'
reqCols = {'Rat','Phase','Session_Number','Performance_All','Data_Dir'};
missing = setdiff(reqCols, summaryTbl.Properties.VariableNames);
if ~isempty(missing)
    error('session_summary.csv missing required column(s): %s', strjoin(missing, ', '));
end

% Normalize/robustly extract columns
ratsCol      = toNum(summaryTbl.Rat);
phaseCol     = toNum(summaryTbl.Phase);
sessNumCol   = toNum(summaryTbl.Session_Number);
perfAllCol   = toNum(summaryTbl.Performance_All);
dataDirCol   = summaryTbl.Data_Dir;

% If Performance_All likely in 0–1, convert to 0–100
if all(perfAllCol(~isnan(perfAllCol)) <= 1)
    perfAllCol = perfAllCol * 100;
end

% Initialize outputs per rat and per phase
sessionPerf = struct([]); trialPerf = struct([]);
for r = 1:numel(rats)
    sessionPerf(r).ratID = rats(r); %#ok<AGROW>
    trialPerf(r).ratID   = rats(r); %#ok<AGROW>
    for p = 1:3
        sessionPerf(r).phase(p).phaseID         = p;
        sessionPerf(r).phase(p).sessionNumbers   = [];
        sessionPerf(r).phase(p).successPercents  = [];
        sessionPerf(r).phase(p).dataDirs         = {};
        trialPerf(r).phase(p).phaseID    = p;
        trialPerf(r).phase(p).allTrials  = []; % rows = sessions, cols = trial index (1..maxLen)
    end
end

% Aggregate session-level and collect per-session trial correctness
for r = 1:numel(rats)
    ratID = rats(r);

    % rows for this rat
    idxRat = ratsCol == ratID & ~isnan(phaseCol) & ~isnan(sessNumCol);

    for p = 1:3
        idxRP = idxRat & (phaseCol == p);

        sessNums = sessNumCol(idxRP);
        perfAll  = perfAllCol(idxRP);
        dirs     = dataDirCol(idxRP);

        % Sort by Session_Number for consistent plotting
        [sessNums, sortIdx] = sort(sessNums(:));
        perfAll = perfAll(sortIdx);
        if iscell(dirs)
            dirs = dirs(sortIdx);
        else
            dirs = dirs(sortIdx,:); %#ok<NASGU>
        end

        sessionPerf(r).phase(p).sessionNumbers  = sessNums(:)'; % row
        sessionPerf(r).phase(p).successPercents = perfAll(:)';  % row
        sessionPerf(r).phase(p).dataDirs        = dirs;

        % --- Build trial-level matrices per rat×phase by reading each session's trial_times.csv
        trialCell = {};
        for s = 1:numel(dirs)
            sessPath = strtrim(asChar(dirs{s}));
            T = load_trial_table(sessPath);
            if isempty(T), continue; end
            correct = strcmpi(T.Result, 'SUCCESS');
            trialCell{end+1} = correct(:)'; %#ok<AGROW>
        end

        if isempty(trialCell)
            trialPerf(r).phase(p).allTrials = [];
        else
            maxLen = max(cellfun(@length, trialCell));
            M = NaN(numel(trialCell), maxLen);
            for k = 1:numel(trialCell)
                lenk = length(trialCell{k});
                M(k,1:lenk) = trialCell{k};
            end
            trialPerf(r).phase(p).allTrials = M;
        end
    end
end

save(cacheFile, 'sessionPerf', 'trialPerf');
log_utils.printf('[INFO] Saved processed data to cache: %s\n', cacheFile);
end

% -------------------------------------------------------------------------
function T = load_trial_table(sessPath)
% LOAD_TRIAL_TABLE
%   Loads trial_times.csv if present, otherwise extracted_biases.csv.
%   Both file types must include a "Result" column ("SUCCESS"/"ERROR").
%   Returns [] if no valid file found.
trialFile = fullfile(sessPath, 'trial_times.csv');
biasFile  = fullfile(sessPath, 'extracted_biases.csv');

if isfile(trialFile)
    T = readtable(trialFile, 'VariableNamingRule','preserve');
elseif isfile(biasFile)
    T = readtable(biasFile, 'VariableNamingRule','preserve');
else
    log_utils.printf('[WARN] No trial_times.csv or extracted_biases.csv in %s\n', sessPath);
    T = [];
    return
end

% Validate presence of Result column
if ~ismember('Result', T.Properties.VariableNames)
    error('File in %s does not contain a "Result" column.', sessPath);
end
end

% -------------------------------------------------------------------------
function plot_success_by_session(sessionPerfRat, ratID, plotCol, posF1, plotXLimF1, figDir)
% One figure per rat.
figure('Color','w','Position', posF1); hold on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
box off;

% Collect all sessions for this rat across phases
allX = []; allY = [];
for p = 1:3
    allX = [allX, sessionPerfRat.phase(p).sessionNumbers]; %#ok<AGROW>
    allY = [allY, sessionPerfRat.phase(p).successPercents]; %#ok<AGROW>
end

% Determine x-limits
if ~isempty(plotXLimF1)
    xlim(plotXLimF1);
else
    if isempty(allX)
        xlim([0 1]);
    else
        xlim([min(allX) max(allX)]);
    end
end
ylim([0 100]);

% Draw phase patches first (so they are behind the data)
yl = ylim;
for p = 1:3
    xvec = sessionPerfRat.phase(p).sessionNumbers;
    if isempty(xvec), continue; end
    xStart = min(xvec);
    xEnd   = max(xvec);
    % expand slightly to encompass markers centered on integer session numbers
    x0 = xStart - 0.5; x1 = xEnd + 0.5;
    patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], plotCol(p,:), ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
end

% Re-plot axes limits in case patches changed auto-scaling
if ~isempty(plotXLimF1)
    xlim(plotXLimF1);
else
    if isempty(allX)
        xlim([0 1]);
    else
        xlim([min(allX) max(allX)]);
    end
end
ylim([0 100]);

% Plot the session data in black, sorted by session number
[allXsorted, ord] = sort(allX);
allYsorted = allY(ord);

plot(allXsorted, allYsorted, '-o', ...
    'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'LineWidth', 1.4, 'MarkerSize', 5);

yline(50, 'k--', 'LineWidth', 1.2);

xlabel('Session #');
ylabel('Success rate (%)');

% Save
outName = fullfile(figDir, sprintf('performance_all_sssion_RAT%02d.svg', ratID));
print(outName, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved session-level plot for Rat %d: %s\n', ratID, outName);
end

% -------------------------------------------------------------------------
function out = analyze_phase1_block_aligned(sessionPerfRat)
% ANALYZE_PHASE1_BLOCK_ALIGNED
% Blocks are defined by FloorCue switches (Green <-> not-Green).
% Trial index resets at each block.
%
% Outputs (new session-centric fields for boxplots):
%   .M_sess     : sessions × trials matrix of per-session trialwise success (%)
%   .mu_sess    : mean across sessions at each trial index (%)
%   .sem_sess   : SEM across sessions at each trial index (%)
%   .n_sess     : number of sessions contributing at each trial index
%   .X, .Y      : long vectors for boxchart (trial index, per-session %)
%
% Backward compatibility (block-centric):
%   .M          : blocks × trials matrix of 0/1 correctness
%   .mu/.sem/.n : aggregates across blocks (0–1 scale)
%   .X_block/.Y_block : long-form from .M (0/1), not used for plotting now

out = init_out_struct(1);

dirs  = get_or_empty(sessionPerfRat, {'phase',1,'dataDirs'});
sessN = get_or_empty(sessionPerfRat, {'phase',1,'sessionNumbers'});
if isempty(dirs), return; end

block_list = {};     % all blocks across sessions
blocks     = [];     % block metadata (flat list)
sess_block_list = cell(numel(dirs),1); % per-session cell array of block vectors

for s = 1:numel(dirs)
    sessPath = strtrim(asCharLocal(dirs{s}));
    T = load_trial_table_block(sessPath);
    if isempty(T), continue; end

    if ~ismember('Result', T.Properties.VariableNames)
        warning('[P1] Missing Result in %s, skipping session.', sessPath);
        continue;
    end

    % Floor cue ON if explicitly "Green" (others treated as OFF)
    if ismember('FloorCue', T.Properties.VariableNames)
        floor_on = strcmpi(string(T.FloorCue), 'Green');
        floor_on(ismissing(T.FloorCue)) = false;
    else
        warning('[P1] Missing FloorCue in %s, treating entire session as single block.', sessPath);
        floor_on = false(height(T),1);
    end

    correct = strcmpi(string(T.Result), 'SUCCESS');

    % Block starts: first trial + any switch in floor_on
    idx_start = [1; find(floor_on(2:end) ~= floor_on(1:end-1)) + 1];
    idx_end   = [idx_start(2:end)-1; height(T)];

    for b = 1:numel(idx_start)
        ii = idx_start(b):idx_end(b);
        v  = double(correct(ii)); % 1/0 correctness
        block_list{end+1} = v; %#ok<AGROW>
        % record under session
        if isempty(sess_block_list{s}), sess_block_list{s} = {v};
        else, sess_block_list{s}{end+1} = v; %#ok<AGROW>
        end

        blk.sessionIdx        = s;
        blk.sessionNumber     = if_nonempty(sessN, s);
        blk.dir               = sessPath;
        blk.blockIdxInSession = b;
        blk.floorOnAtStart    = floor_on(idx_start(b));
        blk.len               = numel(ii);
        blk.trialIdxRange     = [idx_start(b), idx_end(b)];
        blocks = [blocks; blk]; %#ok<AGROW>
    end
end

% -------- Block-centric outputs (retained)
if ~isempty(block_list)
    maxLen = max(cellfun(@numel, block_list));
    M = NaN(numel(block_list), maxLen);
    for k = 1:numel(block_list)
        v = block_list{k};
        M(k,1:numel(v)) = v;
    end
    out.M = M;
    [out.mu, out.sem, out.n] = agg_curve(M);        % 0..1
    [Xb, Yb] = matrix_to_long_for_box(M);           % Xb in trial index, Yb in {0,1}
    out.X_block = Xb;
    out.Y_block = Yb;
end

out.blocks = blocks;
out.block_lengths = arrayfun(@(z) z.len, blocks);
out.meta.phase = 1;
if isfield(sessionPerfRat,'ratID'), out.meta.ratID = sessionPerfRat.ratID; end
out.meta.nBlocks = numel(block_list);

% -------- Session-centric aggregation for boxplots
[out.M_sess, out.mu_sess, out.sem_sess, out.n_sess] = per_session_trial_success(sess_block_list, true); % percent
[Xs, Ys] = session_matrix_to_long_percent(out.M_sess);
out.X = Xs;               % trial index per observation
out.Y = Ys;               % per-session % success
end

% -------------------------------------------------------------------------
function out = analyze_phase3_block_aligned(sessionPerfRat)
% ANALYZE_PHASE3_BLOCK_ALIGNED
% Blocks defined by orientation changes anchored to MazeRotation events.
% Trial index resets at each block.
%
% Outputs (new session-centric fields for boxplots):
%   .M_sess     : sessions × trials matrix of per-session trialwise success (%)
%   .mu_sess    : mean across sessions at each trial index (%)
%   .sem_sess   : SEM across sessions at each trial index (%)
%   .n_sess     : number of sessions contributing at each trial index
%   .X, .Y      : long vectors for boxchart (trial index, per-session %)
%
% Backward compatibility (block-centric):
%   .M          : blocks × trials matrix of 0/1 correctness
%   .mu/.sem/.n : aggregates across blocks (0–1 scale)
%   .X_block/.Y_block : long-form from .M (0/1), not used for plotting now

out = init_out_struct(3);

dirs  = get_or_empty(sessionPerfRat, {'phase',3,'dataDirs'});
sessN = get_or_empty(sessionPerfRat, {'phase',3,'sessionNumbers'});
if isempty(dirs), return; end

block_list = {};
blocks     = [];
sess_block_list = cell(numel(dirs),1);

for s = 1:numel(dirs)
    sessPath = strtrim(asCharLocal(dirs{s}));
    T = load_trial_table_block(sessPath);
    if isempty(T), continue; end

    needCols = {'Result','StartTime'};
    if ~all(ismember(needCols, T.Properties.VariableNames))
        warning('[P3] %s missing required columns (%s). Skipping session.', ...
            sessPath, strjoin(setdiff(needCols, T.Properties.VariableNames), ', '));
        continue;
    end

    correct = strcmpi(string(T.Result), 'SUCCESS');
    ntr = height(T);

    % Orientation label per trial via MazeRotationTime <= StartTime
    haveMR  = ismember('MazeRotation',     T.Properties.VariableNames);
    haveMRT = ismember('MazeRotationTime', T.Properties.VariableNames);
    orient  = repmat(string("<unknown>"), ntr, 1);

    if haveMR && haveMRT
        rawLbl = string(T.MazeRotation);
        rawTim = double(T.MazeRotationTime);
        isEvt  = ~ismissing(rawLbl) & (strlength(strtrim(rawLbl))>0) & ~isnan(rawTim);
        evtLbl = rawLbl(isEvt);
        evtTim = rawTim(isEvt);
        [evtTim, ord] = sort(evtTim(:));
        evtLbl = evtLbl(ord);

        startT = double(T.StartTime);
        j = 1; lastLbl = string("<unknown>");
        for t = 1:ntr
            st = startT(t);
            while j <= numel(evtTim) && evtTim(j) <= st
                lastLbl = strtrim(evtLbl(j));
                j = j + 1;
            end
            orient(t) = lastLbl;
        end
    end

    % Segment blocks at orientation changes
    idx_start = 1 + [0; find(orient(2:end) ~= orient(1:end-1))];
    idx_end   = [idx_start(2:end)-1; ntr];
    start_orient_lbls = orient(idx_start);

    for b = 1:numel(idx_start)
        ii = idx_start(b):idx_end(b);
        v  = double(correct(ii));  % 1/0 correctness
        block_list{end+1} = v; %#ok<AGROW>
        % record under session
        if isempty(sess_block_list{s}), sess_block_list{s} = {v};
        else, sess_block_list{s}{end+1} = v; %#ok<AGROW>
        end

        blk.sessionIdx        = s;
        blk.sessionNumber     = if_nonempty(sessN, s);
        blk.dir               = sessPath;
        blk.blockIdxInSession = b;
        blk.startOrientation  = start_orient_lbls(b);
        blk.startChamberNum   = parse_chamber_num(start_orient_lbls(b));
        blk.len               = numel(ii);
        blk.trialIdxRange     = [idx_start(b), idx_end(b)];
        blocks = [blocks; blk]; %#ok<AGROW>
    end
end

% -------- Block-centric outputs (retained)
if ~isempty(block_list)
    maxLen = max(cellfun(@numel, block_list));
    M = NaN(numel(block_list), maxLen);
    for k = 1:numel(block_list)
        v = block_list{k};
        M(k,1:numel(v)) = v;
    end
    out.M = M;
    [out.mu, out.sem, out.n] = agg_curve(M);        % 0..1
    [Xb, Yb] = matrix_to_long_for_box(M);           % 0/1 (not plotted)
    out.X_block = Xb;
    out.Y_block = Yb;
end

out.blocks = blocks;
out.block_lengths = arrayfun(@(z) z.len, blocks);
out.meta.phase = 3;
if isfield(sessionPerfRat,'ratID'), out.meta.ratID = sessionPerfRat.ratID; end
out.meta.nBlocks = numel(block_list);

% -------- Session-centric aggregation for boxplots
[out.M_sess, out.mu_sess, out.sem_sess, out.n_sess] = per_session_trial_success(sess_block_list, true); % percent
[Xs, Ys] = session_matrix_to_long_percent(out.M_sess);
out.X = Xs;               % trial index per observation
out.Y = Ys;               % per-session % success
end

% -------------------------------------------------------------------------
function [Msess, mu_sess, sem_sess, n_sess] = per_session_trial_success(sess_block_list, asPercent)
% Compute per-session trialwise success by aggregating across that session’s blocks.
% INPUT:
%   sess_block_list : {nSess x 1} with each cell = {b1, b2, ...}, blocks are 0/1 vectors
%   asPercent       : if true, scale outputs to 0..100; else 0..1
% OUTPUT:
%   Msess  : nSessUsed × maxLen matrix, entries are per-session success per trial index
%   mu_sess/sem_sess/n_sess : across-session aggregates at each trial index

if nargin < 2, asPercent = true; end

% Filter to sessions that have at least one block
nonEmptyIdx = find(~cellfun(@isempty, sess_block_list));
nSessUsed   = numel(nonEmptyIdx);
if nSessUsed == 0
    Msess   = [];
    mu_sess = [];
    sem_sess= [];
    n_sess  = [];
    return
end

% Determine max length across all blocks in each session, then across sessions
maxLenPerSess = zeros(nSessUsed,1);
for i = 1:nSessUsed
    B = sess_block_list{nonEmptyIdx(i)};
    maxLenPerSess(i) = max(cellfun(@numel, B));
end
maxLenAll = max(maxLenPerSess);

% Build Msess
Msess = NaN(nSessUsed, maxLenAll);
for i = 1:nSessUsed
    B = sess_block_list{nonEmptyIdx(i)};
    L = cellfun(@numel, B);
    mx = max(L);
    row = NaN(1, mx);
    for t = 1:mx
        vals_t = [];
        for j = 1:numel(B)
            if L(j) >= t
                vals_t(end+1) = B{j}(t); %#ok<AGROW>
            end
        end
        if ~isempty(vals_t)
            row(t) = mean(vals_t, 'omitnan'); % 0..1
        end
    end
    Msess(i,1:mx) = row;
end

% Scale to percent if requested
if asPercent
    Msess = Msess * 100;
end

% Across-session aggregates per trial index
mu_sess = nanmean(Msess, 1);
n_sess  = sum(~isnan(Msess), 1);
sd_sess = nanstd(Msess, [], 1);
denom   = sqrt(n_sess);
denom(denom==0) = 1;
sem_sess = sd_sess ./ denom;
end

% -------------------------------------------------------------------------
function [X, Y] = session_matrix_to_long_percent(Msess)
% Convert sessions×trials matrix of per-session success (%) into long vectors
% X = trial index (1..maxLen) per observation; Y = % success (0..100)
if isempty(Msess)
    X = []; Y = [];
    return
end
[ns, nt] = size(Msess);
X = []; Y = [];
for t = 1:nt
    col = Msess(:,t);
    valid = ~isnan(col);
    if any(valid)
        X = [X; repmat(t, sum(valid), 1)]; %#ok<AGROW>
        Y = [Y; col(valid)];               %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function [X, Y] = matrix_to_long_for_box(M)
% Convert a blocks-by-trials matrix of 0/1/NaN into long vectors:
% X = trial index (1..maxLen) per observation, Y = 0 or 1
[~, nt] = size(M);
X = []; Y = [];
for j = 1:nt
    col = M(:,j);
    col = col(~isnan(col));
    if ~isempty(col)
        X = [X; repmat(j, numel(col), 1)]; %#ok<AGROW>
        Y = [Y; col(:)];                   %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function plot_success_by_trial(outStruct, ratID, phaseID, plotCol, posF2, plotXLimF2, figDir)
% PLOT_SUCCESS_BY_TRIAL (session-level mean ± SEM with scatter)
% For each trial index, plot:
%   - Individual session-level success rates (%)
%   - Mean (square marker) with ± SEM error bars
%   - Line connecting the means

F = figure('Color','w','Position', posF2);
ax = gca; hold(ax,'on'); box(ax,'off'); set(ax,'Color','w');
set(ax, 'FontName','Arial', 'FontSize', 15);

% Ensure we have session-level data
hasSess = isfield(outStruct,'M_sess') && ~isempty(outStruct.M_sess);
if ~hasSess
    text(0.5, 0.5, sprintf('No session-level data (Phase %d)', phaseID), 'Units','normalized', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Color',[0.2 0.2 0.2],'FontSize',12);
    ylim(ax,[0 100]);
    if ~isempty(plotXLimF2), xlim(ax,plotXLimF2); else, xlim(ax,[1 1]); end
    ylabel(ax,'Success rate (%)'); xlabel(ax,'Trial #');
    set(ax,'XTick',1:plotXLimF2(2));
    save_and_log(figDir, sprintf('performance_by_trial_plot_RAT%02d_P%d.svg',ratID,phaseID),F,ratID,phaseID);
    return
end

Msess = outStruct.M_sess; % sessions × trials (%)
[ns, nt] = size(Msess);

% Mean ± SEM
mu = nanmean(Msess,1);
sd = nanstd(Msess,[],1);
n  = sum(~isnan(Msess),1);
sem = sd ./ sqrt(n);

% Restrict to trial range
trialIdx = plotXLimF2(1):plotXLimF2(2);
mu = mu(trialIdx);
sem = sem(trialIdx);
valid = ~isnan(mu);
xvals = trialIdx(valid);
yvals = mu(valid);
semvals = sem(valid);

% Plot mean line
plot(ax, xvals, yvals, '-', 'Color', plotCol(phaseID,:), 'LineWidth', 1.5);

% Error bars (SEM)
errorbar(ax, xvals, yvals, semvals, 'LineStyle','none', ...
    'Color', [0 0 0], 'LineWidth',1.2);

% Mean markers (squares)
plot(ax, xvals, yvals, 's', 'MarkerSize',6, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor', plotCol(phaseID,:));

% Reference line
yline(ax,50,'k--','LineWidth',1.0);

% Axes limits
ylim(ax,[0 100]);
xlim(ax,[plotXLimF2(1)-0.25 plotXLimF2(2)+0.25]);
set(ax,'XTick',1:plotXLimF2(2));

ylabel(ax,'Success rate (%)');
xlabel(ax,'Trial #');

% Save
save_and_log(figDir, sprintf('performance_by_trial_plot_RAT%02d_P%d.svg',ratID,phaseID),F,ratID,phaseID);
end

% -------------------------------------------------------------------------
function save_and_log(figDir, fname, F, ratID, phaseID)
if ~isfolder(figDir), mkdir(figDir); end
outName = fullfile(figDir, fname);
print(F, outName, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved trial-level boxchart for Rat %d Phase %d: %s\n', ratID, phaseID, outName);
end

% -------------------------------------------------------------------------
function s = init_out_struct(phaseID)
s = struct( ...
    'M', [], 'mu', [], 'sem', [], 'n', [], ...                  % block-centric
    'M_sess', [], 'mu_sess', [], 'sem_sess', [], 'n_sess', [], ... % session-centric (percent)
    'X', [], 'Y', [], ...                                        % session-centric long form (percent)
    'X_block', [], 'Y_block', [], ...                             % block-centric long form (0/1)
    'blocks', [], 'block_lengths', [], ...
    'meta', struct('phase', phaseID, 'nBlocks', 0, 'ratID', []) ...
    );
end

% -------------------------------------------------------------------------
function v = get_or_empty(S, pathCells)
v = [];
try
    v = S.(pathCells{1})(pathCells{2}).(pathCells{3});
catch
    try
        v = S.(pathCells{1})(pathCells{2});
    catch
        v = [];
    end
end
end

% -------------------------------------------------------------------------
function y = if_nonempty(x, idx)
if isempty(x), y = []; else, y = x(idx); end
end

% -------------------------------------------------------------------------
function [mu, sem, n] = agg_curve(M)
mu = nanmean(M,1);
n  = sum(~isnan(M),1);
sd = nanstd(M,[],1);
denom = sqrt(n);
denom(denom==0) = 1;
sem = sd ./ denom;
end

% -------------------------------------------------------------------------
function num = parse_chamber_num(lbl)
% Extract numeric chamber ID from a label like "Chamber 5 selected"
% Returns NaN if not found.
lbl = string(lbl);
tokens = regexp(lbl, 'Chamber\s+(\d+)', 'tokens', 'once');
if isempty(tokens)
    num = NaN;
else
    num = str2double(tokens{1});
end
end

% -------------------------------------------------------------------------
function T = load_trial_table_block(sessPath)
trialFile = fullfile(sessPath, 'trial_times.csv');
biasFile  = fullfile(sessPath, 'extracted_biases.csv');

if isfile(trialFile)
    T = readtable(trialFile, 'VariableNamingRule','preserve');
elseif isfile(biasFile)
    T = readtable(biasFile, 'VariableNamingRule','preserve');
else
    log_utils.printf('[WARN] No trial_times.csv or extracted_biases.csv in %s\n', sessPath);
    T = [];
    return
end

% Standardize to string columns where needed
if ismember('Result', T.Properties.VariableNames)
    T.Result = string(T.Result);
end
if ismember('FloorCue', T.Properties.VariableNames)
    T.FloorCue = string(T.FloorCue);
end
if ismember('MazeRotation', T.Properties.VariableNames)
    T.MazeRotation = string(T.MazeRotation);
end
end

% -------------------------------------------------------------------------
function s = asCharLocal(v)
if ischar(v), s = v; elseif isstring(v), s = char(v); else, s = char(string(v)); end
end

% -------------------------------------------------------------------------
function F = plot_trial1_vs_later_bar(outStruct, ratID, phaseID, posF, plotCol, figDir)
% PLOT_TRIAL1_VS_LATER_BAR
% Bar plot of session-level success rate (%) for Trial 1 vs later trials.
% Uses the same per-session block-based aggregation as the stats function.
%
% INPUTS:
%   outStruct   - output of analyze_phase*_block_aligned (must have .blocks and .M)
%   ratID       - numeric rat identifier
%   phaseID     - phase number (1 or 3 supported)
%   posF        - figure position [x y w h]
%   plotCol     - Nx3 color matrix (RGB) for phases
%   figDir  - directory to save figure
%
% OUTPUT:
%   F           - figure handle

if isempty(outStruct) || ~isfield(outStruct,'blocks') || isempty(outStruct.blocks)
    warning('No data for Rat %d Phase %d in plot_trial1_vs_later_bar', ratID, phaseID);
    F = [];
    return;
end

sessNums = unique([outStruct.blocks.sessionIdx]);
t1_vals = nan(numel(sessNums),1);
later_vals = nan(numel(sessNums),1);

% --- Compute per-session means
for si = 1:numel(sessNums)
    idx = find([outStruct.blocks.sessionIdx] == sessNums(si));
    trial1_all = [];
    later_all  = [];
    for k = idx
        blk = outStruct.blocks(k);
        v   = outStruct.M(k, 1:blk.len); % 0/1 vector for this block
        if isempty(v), continue; end
        trial1_all(end+1) = v(1); %#ok<AGROW>
        if blk.len > 1
            later_all = [later_all, v(2:end)]; %#ok<AGROW>
        end
    end
    if ~isempty(trial1_all)
        t1_vals(si) = mean(trial1_all)*100;
    end
    if ~isempty(later_all)
        later_vals(si) = mean(later_all)*100;
    end
end

% --- Gather stats
means = [mean(t1_vals,'omitnan'), mean(later_vals,'omitnan')];
sems  = [std(t1_vals,'omitnan')/sqrt(sum(~isnan(t1_vals))), ...
         std(later_vals,'omitnan')/sqrt(sum(~isnan(later_vals)))];

% --- Plot
F = figure('Color','w','Position',posF); hold on;

x_pos = 1:2;
barWidth = 0.6;
colors = [plotCol(phaseID,:); plotCol(phaseID,:)]; 

for j = 1:2
    if isfinite(means(j))
        bar(x_pos(j), means(j), barWidth, 'FaceColor', colors(j,:));
    else
        bar(x_pos(j), 0, barWidth, 'FaceColor', colors(j,:), ...
            'FaceAlpha',0.15,'EdgeAlpha',0.15);
    end
end

% Error bars (SEM)
errorbar(x_pos, means, sems, 'k.', 'CapSize',6, 'LineWidth',1.5);

% Labels & formatting
set(gca,'XTick',x_pos,'XTickLabel',{'1','2+'});
set(gca,'FontName','Arial','FontSize',15);
xlabel('Trial #');
ylim([0 100]);
yline(50,'k--','LineWidth',1.2);
xlim([0.5, 2.5]);

% Remove y-axis fully
set(gca,'YColor','none');   % remove y-axis line & ticks
ylabel('');                 % clear label

% Save
fname = fullfile(figDir, sprintf('performance_by_trial_bar_RAT%02d_P%d.svg',ratID,phaseID));
print(F,fname,'-dsvg','-r300');
log_utils.printf('[INFO] Saved Trial1-vs-Later bar plot for Rat %d Phase %d: %s\n', ratID, phaseID, fname);
end

% -------------------------------------------------------------------------
function out = toNum(x)
% Robustly convert table column to numeric
if isnumeric(x)
    out = x;
elseif iscell(x)
    out = cellfun(@(v) str2double(asChar(v)), x);
elseif isstring(x)
    out = str2double(x);
elseif iscategorical(x)
    out = double(x);
else
    out = str2double(string(x));
end
end

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
function print_behavior_stats(sessionPerf, trialPerf, rats)
log_utils.printf('\n================== BEHAVIOR PERFORMANCE SUMMARY ==================\n');

for i = 1:numel(rats)
    ratID = sessionPerf(i).ratID;
    log_utils.printf('\n--------------------------- RAT %02d -----------------------------\n', ratID);

    phasePerfSessions = cell(1,3);

    for p = 1:3
        y    = sessionPerf(i).phase(p).successPercents;     % session-level % (0..100)
        Mraw = trialPerf(i).phase(p).allTrials;             % trials matrix (0/1/NaN)

        log_utils.printf('Phase %d\n', p);

        % ---------------- Session-level vs 50% ----------------
        if ~isempty(y)
            mu  = mean(y,'omitnan');
            sd  = std(y,'omitnan');
            n_s = sum(~isnan(y));
            log_utils.printf('  Session-level: mean = %.2f%% ± %.2f (n = %d)\n', mu, sd, n_s);

            [~, pSess, ci, statsSess] = ttest(y, 50, 'Tail','both');
            diffs = y - 50;
            d = mean(diffs) / std(diffs); % Cohen’s d
            log_utils.printf('  One-sample two-tailed t test vs 50%%: Δmean = %.2f pp, t(%d) = %.2f, p = %.3f, d = %.2f, 95%% CI [%.2f, %.2f]\n', ...
                mean(diffs), statsSess.df, statsSess.tstat, pSess, d, ci(1), ci(2));
            phasePerfSessions{p} = y;
        else
            log_utils.printf('  Session-level: n = 0\n');
            phasePerfSessions{p} = [];
        end

        % ---------------- Trial-level vs 50% ----------------
        if ~isempty(Mraw)
            valid = ~isnan(Mraw);
            k = sum(Mraw(valid));          
            N = sum(valid(:));             
            prop = 100 * (k / max(N,1));
            log_utils.printf('  Trial-level: overall = %.2f%% (%d/%d trials)\n', prop, k, N);

            pTrial = safe_binotest(k, N, 0.5);
            % For binomial, also approximate effect size with Cohen's h
            phat = k/N; h = 2*asin(sqrt(phat)) - 2*asin(sqrt(0.5));
            log_utils.printf('  Exact binomial test vs 50%%: p = %.3f, Cohen''s h = %.2f\n', pTrial, h);
        else
            log_utils.printf('  Trial-level: N = 0\n');
        end

        % ---------------- Within-phase improvement (first half vs last half of sessions)
        if numel(y) >= 4
            fh = y(1:floor(numel(y)/2));
            lh = y(floor(numel(y)/2)+1:end);

            [~, pImp, ci, statsImp] = ttest2(fh, lh, 'Vartype','unequal'); % two-sample, Welch’s by default

            % Effect size: Cohen’s d for two independent groups
            m1 = mean(fh); m2 = mean(lh);
            s1 = std(fh); s2 = std(lh);
            n1 = numel(fh); n2 = numel(lh);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
            d = (m2 - m1) / sp;

            log_utils.printf('  Within-phase change (first vs last half): Δmean = %.2f pp, t(%d) = %.2f, p = %.3f, d = %.2f, 95%% CI [%.2f, %.2f]\n', ...
                (m2 - m1), statsImp.df, statsImp.tstat, pImp, d, ci(1), ci(2));
        else
            log_utils.printf('  Within-phase change: insufficient data\n');
        end

        % ---------------- Trial 1 vs later ----------------
        if p==1
            outP = analyze_phase1_block_aligned(sessionPerf(i));
        elseif p==3
            outP = analyze_phase3_block_aligned(sessionPerf(i));
        else
            outP = [];
        end
        if ~isempty(outP) && isfield(outP,'blocks') && ~isempty(outP.blocks)
            sessIdxs = unique([outP.blocks.sessionIdx]);
            t1_vals = nan(numel(sessIdxs),1);
            later_vals = nan(numel(sessIdxs),1);

            for si = 1:numel(sessIdxs)
                idx = find([outP.blocks.sessionIdx] == sessIdxs(si));
                trial1_all = [];
                later_all  = [];
                for kblk = idx
                    blk = outP.blocks(kblk);
                    v   = outP.M(kblk, 1:blk.len); % 0/1
                    if isempty(v), continue; end
                    trial1_all(end+1) = v(1);
                    if blk.len > 1
                        later_all = [later_all, v(2:end)];
                    end
                end
                if ~isempty(trial1_all), t1_vals(si)    = mean(trial1_all)*100; end
                if ~isempty(later_all),  later_vals(si) = mean(later_all)*100;  end
            end

            valid = ~isnan(t1_vals) & ~isnan(later_vals);
            if nnz(valid) >= 2
                [~, pPair, ci, statsPair] = ttest(t1_vals(valid), later_vals(valid));
                diffs = later_vals(valid) - t1_vals(valid);
                d = mean(diffs) / std(diffs);
                log_utils.printf('  Trial 1 vs 2+: Δmean = %.2f pp, t(%d) = %.2f, p = %.3f, d = %.2f, 95%% CI [%.2f, %.2f]\n', ...
                    mean(diffs), statsPair.df, statsPair.tstat, pPair, d, ci(1), ci(2));
            else
                log_utils.printf('  Trial 1 vs 2+: insufficient data\n');
            end
        else
            log_utils.printf('  Trial 1 vs 2+: not computed\n');
        end

        log_utils.printf('\n');
    end

    % ---------------- Phase 1 vs Phase 3 comparison ----------------
    if ~isempty(phasePerfSessions{1}) && ~isempty(phasePerfSessions{3})
        y1 = phasePerfSessions{1};
        y3 = phasePerfSessions{3};

        [~, p13, ci, stats13] = ttest2(y1, y3, 'Vartype','unequal'); % Welch’s t test

        % Effect size: Cohen’s d (independent groups, pooled SD)
        m1 = mean(y1); m3 = mean(y3);
        s1 = std(y1);  s3 = std(y3);
        n1 = numel(y1); n3 = numel(y3);
        sp = sqrt(((n1-1)*s1^2 + (n3-1)*s3^2) / (n1+n3-2));
        d = (m3 - m1) / sp;

        log_utils.printf('Phase 1 vs Phase 3 (session %% correct): Δmean = %.2f pp, t(%d) = %.2f, p = %.3f, d = %.2f, 95%% CI [%.2f, %.2f]\n', ...
            (m3 - m1), stats13.df, stats13.tstat, p13, d, ci(1), ci(2));
    else
        log_utils.printf('Phase 1 vs Phase 3: insufficient data\n');
    end

end

log_utils.printf('\n==================================================================\n\n');
end

% ----------------------------- HELPERS -----------------------------------
function stars = p_to_stars(p)
if isnan(p)
    stars = '';
elseif p < 1e-3
    stars = '***';
elseif p < 1e-2
    stars = '**';
elseif p < 5e-2
    stars = '*';
else
    stars = 'ns';
end
end

function p = safe_binotest(k, n, p0)
% Two-sided exact binomial p-value vs p0. Uses MATLAB binotest when available,
% otherwise falls back to a conservative 2*min-tail method via binocdf.
if n <= 0
    p = NaN; return;
end
if exist('binotest','file') == 2
    % R2022a+: binotest(k,n,p0) returns a hypothesis-test object; get p-value.
    try
        ht = binotest(k, n, 'p', p0, 'Tail','two-sided');
        p  = ht.pValue;
        return;
    catch
        % fall through to manual
    end
end
% Fallback two-sided p-value
p_lower = binocdf(k, n, p0);
if k>0
    p_upper = 1 - binocdf(k-1, n, p0);
else
    p_upper = 1;
end
p = 2 * min(p_lower, p_upper);
p = min(max(p,0),1);
end

function s = phase_list_str(v)
if isempty(v)
    s = 'none';
else
    s = sprintf('P%s', strjoin(string(v), ', P'));
end
end



