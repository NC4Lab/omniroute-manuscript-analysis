function behavior_trajectory_analysis()
% =========================================================================
% FUNCTION: behavior_trajectory_analysis
%
% Description:
%   End-to-end script to:
%     (A) Plot ALL raw (un-normalized) trajectories overlaid in one figure
%         using per-start colors.
%     (B) Plot 4 separate figures with NORMALIZED trajectories, one per
%         start condition, with start mapped to LEFT and correct goal to RIGHT.
%     (C) Plot a BAR CHART of success rate by start condition (no error bars).
%
% Key refactor points:
%   - Normalization is modular: we augment 'trials' with Xn/Yn (do not
%     overwrite raw X/Y). Plotting chooses raw vs normalized via a flag.
%   - Plotting is a pure renderer: caller passes Normalize flag, StartIdx
%     filter, a single RGB base color, and a target axes handle.
%   - Saving is handled OUTSIDE the plotting function.
% =========================================================================

close all

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'behavior_trajectory_test';
figDir = path_utils.results_figures(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% ------------------------- USER SETTINGS ---------------------------------
% Pick one session
sessPath = fullfile(path_utils.dataset_root(), 'processed_behavior\NC40023\250628');
ttimFiPath = fullfile(sessPath, 'trial_times.csv');
posFiPath   = fullfile(sessPath, 'position_times.csv');

outRawName   = 'trajectories_raw_overlay.pdf';
outNormNames = { ...
    'trajectories_norm_start1_top.pdf', ...
    'trajectories_norm_start2_right.pdf', ...
    'trajectories_norm_start3_bottom.pdf', ...
    'trajectories_norm_start4_left.pdf' ...
};
outBarName = 'success_rate_by_start.svg';

% Colors (index: Top=1, Right=2, Bottom=3, Left=4)
plotCol = [ ...
    1.00 0.65 0.00;  % 1 - Top    (orange)
    0.00 0.65 0.65;  % 2 - Right  (teal)
    0.35 0.70 0.30;  % 3 - Bottom (green)
    0.60 0.40 0.80]; % 4 - Left   (purple)

% Visuals
posF1        = [100, 100, 279, 279];   % square figure
posF2        = [450, 100, 186, 279];   % square figure
posF3        = [700, 100, 279, 279];   % for bar chart (can reuse posF1)
lineWidth    = 0.5;
brightFactor = 1.25;   % HSV-V multiplier for correct (brighten)
dimFactor    = 0.55;   % HSV-V multiplier for error (dim)
alphaVal     = 1.0;    % line alpha (use <1 for translucency)

% Arena limits in meters (keep square and consistent with schematics)
plotXLimF1 = [0 0.9];
plotYLimF1 = [0 0.9];
plotXLimF2 = [0 0.6];
plotYLimF2 = [0 0.9];

% Start inference: use first N samples within the trial window
N_START_SAMPLES = 40;     % take earliest N samples to infer start
MIN_SAMPLES_TR  = 20;     % minimum samples required to plot a trial

% Centers of the four start chambers (Top, Right, Bottom, Left)
startCenters = [ ...
    0.45 0.75;   % Top
    0.75 0.45;   % Right
    0.45 0.15;   % Bottom
    0.15 0.45];  % Left

% --------------------------- PIPELINE ------------------------------------
T_trials = readtable(ttimFiPath, 'VariableNamingRule','preserve');
T_pos    = readtable(posFiPath,   'VariableNamingRule','preserve');

log_basic_info(T_trials, T_pos, sessPath);

% Robust time columns
tStartCol = pick_col(T_trials, {'StartTime','TrialStart','Start_TS','Start'});
tEndCol   = pick_col(T_trials, {'EndTrialTime','EndTime','TrialEnd','End'});
if isempty(tStartCol)
    error('No recognizable StartTime column found in trial_times.csv');
end
if isempty(tEndCol)
    log_utils.printf('[WARN] No EndTrialTime/EndTime found. Using last available sample per trial as end.\n');
end

% Extract raw trial segments (raw X/Y retained)
trials = extract_trial_segments(T_trials, T_pos, tStartCol, tEndCol, ...
    startCenters, N_START_SAMPLES, MIN_SAMPLES_TR);

% Augment with normalized coordinates (Xn/Yn) using trial metadata
trials = augment_trials_with_normalization(trials, T_trials, plotXLimF1, plotYLimF1);

% -------------------- (A) RAW OVERLAY: all starts ------------------------

figRaw = figure('Color','none'); set(figRaw, 'Position', posF1);
axRaw  = axes('Parent',figRaw); hold(axRaw, 'on');
for si = 1:4
    plot_trajs(trials, false, si, plotCol(si,:), axRaw, lineWidth, brightFactor, dimFactor, alphaVal);
end
style_axes(axRaw, plotXLimF1, plotYLimF1);
exportgraphics(axRaw, fullfile(figDir, outRawName), ...
    'BackgroundColor','none', 'ContentType','vector');
log_utils.printf('[INFO] Saved RAW overlay plot: %s\n', fullfile(figDir, outRawName));

% -------------------- (B) NORMALIZED: one fig per start ------------------
for si = 1:4
    figN = figure('Color','none'); set(figN, 'Position', posF2);
    axN  = axes('Parent',figN); hold(axN, 'on');
    plot_trajs(trials, true, si, plotCol(si,:), axN, lineWidth, brightFactor, dimFactor, alphaVal);
    style_axes(axN, plotXLimF2, plotYLimF2);
    exportgraphics(axN, fullfile(figDir, outNormNames{si}), ...
        'BackgroundColor','none', 'ContentType','vector');
    log_utils.printf('[INFO] Saved NORMALIZED plot (start=%s): %s\n', idx_to_name(si), fullfile(figDir, outNormNames{si}));
end

% -------------------- (C) BAR: success rate by start ---------------------
figBar = success_rate_bar_plot(trials, posF3, plotCol);
print(figBar, fullfile(figDir, outBarName), '-dsvg', '-r300');
log_utils.printf('[INFO] Saved success-rate bar plot: %s\n', fullfile(figDir, outBarName));

% -------------------- Summary log ----------------------------------------
summarize_trials(trials);
log_utils.printf('[DONE] All trajectory and bar plots complete.\n');

log_utils.close()
end % behavior_trajectory_analysis

% -------------------------------------------------------------------------
function colName = pick_col(T, candidates)
% Return the first matching column name from candidates (case-insensitive),
% or [] if none exist.
vnames = T.Properties.VariableNames;
lowerMap = containers.Map(lower(vnames), vnames);
colName = [];
for i = 1:numel(candidates)
    key = lower(candidates{i});
    if lowerMap.isKey(key)
        colName = lowerMap(key);
        return
    end
end
end

% -------------------------------------------------------------------------
function trials = extract_trial_segments(T_trials, T_pos, tStartCol, tEndCol, ...
    startCenters, N_START_SAMPLES, MIN_SAMPLES_TR)

nT = height(T_trials);
trials = repmat(struct('X',[],'Y',[],'Result',"UNKNOWN",'StartIdx',NaN,'StartName','', ...
                       't0',NaN,'t1',NaN), nT, 1);

% Pre-pull arrays for speed
tPos = T_pos.Time;
xPos = T_pos.X;
yPos = T_pos.Y;

for k = 1:nT
    % ---------------- timestamps ----------------
    t0 = T_trials.(tStartCol)(k);
    t1 = [];
    if ~isempty(tEndCol)
        t1 = T_trials.(tEndCol)(k);
    end
    if isnan(t0) || (~isempty(t1) && isnan(t1))
        continue
    end

    % If t1 missing, use last time before next trial start (or last pos time)
    if isempty(t1) || ~isfinite(t1)
        if k < nT && isfinite(T_trials.(tStartCol)(k+1))
            t1 = T_trials.(tStartCol)(k+1) - 1; % 1 tick before next start
        else
            t1 = max(tPos);
        end
    end
    if t1 <= t0
        % invalid window; skip
        continue
    end

    % ---------------- segment positions ----------------
    idx = tPos >= t0 & tPos <= t1;
    if ~any(idx), continue; end

    x = xPos(idx);
    y = yPos(idx);
    if numel(x) < MIN_SAMPLES_TR
        continue
    end

    % ---------------- result ----------------
    if ismember('Result', T_trials.Properties.VariableNames)
        trials(k).Result = string(T_trials.Result(k));
    end

    % ---------------- infer start side ----------------
    nTake = min(N_START_SAMPLES, numel(x));
    xs = x(1:nTake);
    ys = y(1:nTake);

    % nearest of the four start centers
    centers = startCenters;
    d2 = (xs(:) - centers(:,1)').^2 + (ys(:) - centers(:,2)').^2; % [nTake x 4]
    mean_d2 = mean(d2, 1);
    [~, startIdx] = min(mean_d2);

    trials(k).X = x;
    trials(k).Y = y;
    trials(k).StartIdx = startIdx;
    trials(k).StartName = idx_to_name(startIdx);
    trials(k).t0 = t0;
    trials(k).t1 = t1;
end

% prune empty entries
trials = trials(~cellfun(@isempty, {trials.X}));

end

% -------------------------------------------------------------------------
function name = idx_to_name(i)
switch i
    case 1, name = 'Top';
    case 2, name = 'Right';
    case 3, name = 'Bottom';
    case 4, name = 'Left';
    otherwise, name = 'Unknown';
end
end

% -------------------------------------------------------------------------
function plot_trajs(trials, Normalize, StartIdxFilter, baseColor, ax, lineWidth, brightFactor, dimFactor, alphaVal)
% PLOT_TRAJS
% Pure renderer. Caller controls:
%   - Normalize: false→use raw X/Y; true→use Xn/Yn
%   - StartIdxFilter: [] or 1..4 (plot only that start)
%   - baseColor: 1x3 RGB to be brightened/dimmed by Success/Error
%   - ax: target axes handle
%   - Saving handled by caller.

if nargin < 3 || isempty(StartIdxFilter)
    useMask = true(1, numel(trials));
else
    useMask = arrayfun(@(t) t.StartIdx==StartIdxFilter, trials);
end

if Normalize
    hasNorm = arrayfun(@(t) isfield(t,'Xn') && ~isempty(t.Xn), trials);
    useMask = useMask & hasNorm;
end

trials_to_plot = trials(useMask);
for k = 1:numel(trials_to_plot)
    tr = trials_to_plot(k);
    if Normalize
        x = tr.Xn; y = tr.Yn;
    else
        x = tr.X;  y = tr.Y;
    end

    % choose brightness by outcome
    if strcmpi(tr.Result, 'SUCCESS')
        c = brighten_rgb(baseColor, brightFactor);
    elseif strcmpi(tr.Result, 'ERROR')
        c = brighten_rgb(baseColor, dimFactor);
    else
        c = baseColor; % unknown
    end

    h = plot(ax, x, y, 'LineWidth', lineWidth);
    h.Color = [c, alphaVal];
end

end

% -------------------------------------------------------------------------
function c_out = brighten_rgb(c_in, factor)
% Adjust brightness via HSV V-channel scaling, with clamping.
hsv = rgb2hsv(c_in);
hsv(3) = min(1, max(0, hsv(3) * factor));
c_out  = hsv2rgb(hsv);
end

% -------------------------------------------------------------------------
function style_axes(ax, plotXLim, plotYLim)
% Format axes with transparent background, correct limits, and
% aspect ratio based on the ranges of x and y.
xlim(ax, plotXLim);
ylim(ax, plotYLim);
dx = diff(plotXLim);
dy = diff(plotYLim);
pbaspect(ax, [dx dy 1]);

% Add invisible bounding box to enforce export area
rectangle(ax, 'Position', [plotXLim(1), plotYLim(1), dx, dy], ...
    'EdgeColor', [0 0 0 0.01], 'FaceColor', 'none', 'HitTest', 'off');

axis(ax, 'off');
set(ax,'Color','none');        % transparent axes background
set(ax.Parent,'Color','none'); % transparent figure background
end

% -------------------------------------------------------------------------
function F = success_rate_bar_plot(trials, posF, plotCol)
% SUCCESS_RATE_BAR_PLOT
% Compute and plot success rate by StartIdx (1..4) as percentages.
%
% INPUTS
%   trials  - struct array with fields .Result ('SUCCESS'/'ERROR') and .StartIdx
%   posF    - figure position [x y w h]
%   plotCol - 4x3 matrix of RGB colors, one per start condition
%
% OUTPUT
%   F - figure handle

% --- Compute success rates ---
res    = upper(string({trials.Result}));
starts = [trials.StartIdx];

C = 4;                          % number of start conditions
succ_rate = nan(1,C);           % success rate per condition
n_per     = zeros(1,C);         % trial counts per condition

for c = 1:C
    mask = (starts == c);
    n = sum(mask);
    n_per(c) = n;
    if n > 0
        n_succ = sum(res(mask) == "SUCCESS");
        succ_rate(c) = (n_succ / n) * 100;  % convert to percentage
    end
end

% --- Plot ---
F = figure('Color','w'); 
set(F,'Position',posF); 
ax = axes('Parent',F); 
hold(ax,'on');

x_pos   = 1:C;
barWidth = 0.6;

for i = 1:C
    if ~isnan(succ_rate(i))
        bar(ax, x_pos(i), succ_rate(i), barWidth, 'FaceColor', plotCol(i,:));
    else
        % no data: draw a faint placeholder for spacing
        bar(ax, x_pos(i), 0, barWidth, 'FaceColor', plotCol(i,:), ...
            'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
    end
end

% --- Axes formatting ---
set(ax, 'XTick', x_pos, 'XTickLabel', {'N','E','S','W'});
set(ax, 'FontName', 'Arial', 'FontSize', 15);
ylabel(ax, 'Success rate (%)');
xlabel(ax, 'Start condition');
ylim(ax, [0, 100]); 
xlim(ax, [0.5, C + 0.5]);
box(ax,'off');
yline(50, 'k--', 'LineWidth', 1.2);

% --- Annotate counts above bars ---
for i = 1:C
    if isfinite(succ_rate(i))
        text(ax, x_pos(i), succ_rate(i) + 3, sprintf('n=%d', n_per(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11);
    end
end

end

% -------------------------------------------------------------------------
function trials_out = augment_trials_with_normalization(trials, T_trials, plotXLimF1, plotYLimF1)
% AUGMENT_TRIALS_WITH_NORMALIZATION (corrected for categorical cue labels)
% Compute normalized coordinates (Xn/Yn) per trial without overwriting X/Y.
% Normalization rules:
%   1) Rotate so StartIdx maps to LEFT (canonical start).
%   2) Ensure the CORRECT goal (from cues) is to the RIGHT; if cues imply
%      the correct goal is on the LEFT after rotation, apply a VERTICAL
%      mirror (about y = cy) to satisfy the convention.
%
% CSV format (as observed):
%   LeftCue, RightCue = 'Triangle' or 'No_Cue'
%   FloorCue          = 'Green' (floor ON) or 'Black' (floor OFF)

% center for rotation/mirroring
cx = mean(plotXLimF1); %#ok<NASGU> % kept for symmetry with helpers (not used in vertical flip)
cy = mean(plotYLimF1);

% column references
col.Start    = pick_col(T_trials, {'StartTime','TrialStart','Start_TS','Start'});
col.LeftCue  = pick_col(T_trials, {'LeftCue'});
col.RightCue = pick_col(T_trials, {'RightCue'});
col.Floor    = pick_col(T_trials, {'FloorCue'});

if isempty(col.Start)
    error('[augment_trials_with_normalization] Could not locate StartTime column in trial_times table.');
end
if isempty(col.LeftCue) || isempty(col.RightCue) || isempty(col.Floor)
    error('[augment_trials_with_normalization] Missing one or more required cue columns (LeftCue, RightCue, FloorCue).');
end

tstarts_tbl = T_trials.(col.Start);
trials_out  = trials;

for i = 1:numel(trials)
    % map to trial row (nearest StartTime within tolerance)
    r = find_row_by_start_time(trials(i).t0, tstarts_tbl);

    % Parse categorical cue strings directly (case-insensitive)
    r_left  = false;
    r_right = false;
    r_floor = false;
    if ~isnan(r)
        r_left  = strcmpi(string(T_trials.(col.LeftCue)(r)),  "Triangle");
        r_right = strcmpi(string(T_trials.(col.RightCue)(r)), "Triangle");
        r_floor = strcmpi(string(T_trials.(col.Floor)(r)),    "Green");
    end

    % rotation by start side to map start→LEFT
    theta = 0; % deg
    switch trials(i).StartIdx
        case 1 % Top    -> +90
            theta = +90;
        case 2 % Right  -> 180
            theta = 180;
        case 3 % Bottom -> -90
            theta = -90;
        case 4 % Left   -> 0
            theta = 0;
        otherwise
            theta = 0;
    end

    [xr, yr] = rot_about_center(trials(i).X, trials(i).Y, theta, mean(plotXLimF1), cy);

    % derive correct side (LEFT/RIGHT) from cues
    corrSide = "UNKNOWN";
    if ~isnan(r)
        if r_floor
            % floor ON -> go to triangles
            if r_right && ~r_left, corrSide = "RIGHT"; end
            if r_left  && ~r_right, corrSide = "LEFT";  end
        else
            % floor OFF -> go to dark (opposite of triangles)
            if r_right && ~r_left, corrSide = "LEFT";  end
            if r_left  && ~r_right, corrSide = "RIGHT"; end
        end
    end

    % VERTICAL flip if correct side is LEFT (so it becomes RIGHT in canonical frame)
    flipped = false;
    if corrSide == "LEFT"
        yr = (2*cy) - yr; % mirror across horizontal midline (vertical flip)
        flipped = true;
    end

    % persist normalized fields
    trials_out(i).Xn = xr;
    trials_out(i).Yn = yr;
    trials_out(i).NormThetaDeg = theta;
    trials_out(i).NormFlip = flipped;
    trials_out(i).CorrSide = corrSide;
end

end

% -------------------------------------------------------------------------
function r = find_row_by_start_time(t0, tstarts_tbl)
% Find nearest StartTime row to t0 within tolerance
if ~isfinite(t0)
    r = NaN; return
end
tstarts = double(tstarts_tbl);
[delta, idx] = min(abs(tstarts - double(t0)));
tol = 1.0; % seconds tolerance (epoch or seconds)
if isempty(idx) || ~isfinite(delta) || delta > tol
    r = NaN;
else
    r = idx;
end
end

% -------------------------------------------------------------------------
function [xo, yo] = rot_about_center(x, y, theta_deg, cx, cy)
% Rotate (x,y) by theta_deg about (cx,cy).
if theta_deg == 0
    xo = x; yo = y; return
end
th = deg2rad(theta_deg);
ct = cos(th); st = sin(th);
x0 = x - cx; y0 = y - cy;
xo =  x0*ct - y0*st + cx;
yo =  x0*st + y0*ct + cy;
end

% -------------------------------------------------------------------------
function log_basic_info(T_trials, T_pos, sessPath)
log_utils.printf('\n[INFO] Session: %s\n', sessPath);
log_utils.printf('[INFO] trial_times.csv columns:\n  ');
log_utils.printf('%s ', T_trials.Properties.VariableNames{:}); log_utils.printf('\n');
log_utils.printf('[INFO] position_times.csv columns:\n  ');
log_utils.printf('%s ', T_pos.Properties.VariableNames{:}); log_utils.printf('\n');

% Time ranges
if ismember('Time', T_pos.Properties.VariableNames)
    tmin = min(T_pos.Time); tmax = max(T_pos.Time);
    log_utils.printf('[INFO] Position time range: [%0.0f  ..  %0.0f]\n', tmin, tmax);
end

% Basic counts
nTrials = height(T_trials);
nPos    = height(T_pos);
log_utils.printf('[INFO] Rows: trials=%d, position=%d\n', nTrials, nPos);
end

% -------------------------------------------------------------------------
function summarize_trials(trials)
if isempty(trials)
    log_utils.printf('[WARN] No trials with position samples passed filters.\n');
    return
end

% Outcome counts
res = upper(string({trials.Result}));
nCorr = sum(res == "SUCCESS");
nErr  = sum(res == "ERROR");
nU    = sum(~(res == "SUCCESS" | res == "ERROR"));

% Start-side indices (1=Top/N, 2=Right/E, 3=Bottom/S, 4=Left/W)
idxs = [trials.StartIdx];
labels = {'N','E','S','W'};

log_utils.printf('\n=========== Trajectory Summary ===========\n');
log_utils.printf('Trials plotted: %d\n', numel(trials));
log_utils.printf('  Overall outcomes: Correct=%d | Error=%d | Unknown=%d | SuccessRate=%.2f\n', ...
    nCorr, nErr, nU, nCorr / max(1,(nCorr+nErr)));

% Per-start breakdown
for si = 1:4
    mask = (idxs == si);
    n_tot = sum(mask);
    if n_tot > 0
        res_si = res(mask);
        n_corr = sum(res_si == "SUCCESS");
        n_err  = sum(res_si == "ERROR");
        n_unk  = sum(~(res_si == "SUCCESS" | res_si == "ERROR"));
        rate   = n_corr / max(1,(n_corr+n_err));
        log_utils.printf('  Start %s: %3d trials  (Correct=%d | Error=%d | Unknown=%d | SuccessRate=%.2f)\n', ...
            labels{si}, n_tot, n_corr, n_err, n_unk, rate);
    else
        log_utils.printf('  Start %s:   0 trials\n', labels{si});
    end
end

% Trial time stats
t0s = [trials.t0]; t1s = [trials.t1];
log_utils.printf('  Trial time windows (ticks): median dur = %0.0f\n', median(t1s - t0s, 'omitnan'));

% Normalization info (optional)
if isfield(trials, 'CorrSide')
    cs = string({trials.CorrSide});
    log_utils.printf('  CorrSide counts  Right=%d | Left=%d | Unknown=%d\n', ...
        sum(cs=="RIGHT"), sum(cs=="LEFT"), sum(cs=="UNKNOWN"));
end
log_utils.printf('==========================================\n\n');
end

