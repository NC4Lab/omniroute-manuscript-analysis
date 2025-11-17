function gantry_movement_analysis()
% =========================================================================
% FUNCTION: gantry_movement_analysis
%
% Description:
%   Analyze gantry movement timing and tracking using ros_utils.
%   - Loads required topics from a ROS bag recorded by gantry_manuscript_testing.py
%   - Decodes commanded targets and aligns them to the pose frame
%   - Computes move latency and accuracy (via ros_utils.analyze_gantry_move)
%   - Generates a cleaned target-vs-pose time-series plot (X and Y)
%   - Saves all intermediate tables (.mat and .csv) to disk
% =========================================================================

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'gantry_movement_test';
bagPath = fullfile(path_utils.dataset_root(), out_dir_name, 'ros_session_20250815_211524.bag');
figDir = path_utils.results_figures(out_dir_name);
tabDir = path_utils.results_tables(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% ---------------------- USER SETTINGS (edit here) ------------------------

% Velocity plots
fs_resample = 200;           % Hz, uniform resampling rate for pose
smooth_win  = 5;             % samples, moving average on velocit

% X-axis plotting range (explicit, editable)
plotXLimF1 = [0, 40];
plotXLimF2 = [0.0, 2.0];
plotXLimF3 = [0, 20];

% Figure colors
plotColF1 = [...
    1.0 0.3 0.3;   % Pose
    0.5 0.5 0.5];  % Target
plotColF2 = [1.0 0.3 0.3];
plotColF3 = [ ...
    1.00 0.65 0.00;  % Move 1 - orange
    0.00 0.65 0.65;  % Move 2 - teal
    0.35 0.70 0.30;  % Move 3 - green
    0.60 0.40 0.80]; % Move 4 - purple
plotColF4 = [ ...
    1.00 0.65 0.00;  % Move 1 - orange
    0.00 0.65 0.65;  % Move 2 - teal
    0.35 0.70 0.30;  % Move 3 - green
    0.60 0.40 0.80]; % Move 4 - purple

% Figure layouts
posF1 = [100, 100, 300, 300];
posF2 = [400, 100, 300, 300];
posF3 = [100, 600, 450, 300];
posF4 = [600, 600, 300, 300];

% Left and right axis for scatter/boxchart
posAxLeft  = [0.06 0.1 0.575 0.9];
posAxRight = [0.575 0.1 0.45 0.9];

% -------------------------- LOAD ROS BAG ----------------------------------
log_utils.printf('[INFO] Loading ROS bag: %s\n', bagPath);
bag = rosbag(bagPath);

% --------------------------- PARSE TOPICS ---------------------------------
% Pose
gantry_pose_tbl = ros_utils.parse_pose_topic(bag, '/gantry_pose_in_maze');
if isempty(gantry_pose_tbl) || height(gantry_pose_tbl)==0
    error('No /gantry_pose_in_maze messages found.');
end
ros_utils.save_table(fullfile(tabDir, 'gantry_pose_tbl'), gantry_pose_tbl);

% EASE write/read for gantry
gantry_write_tbl = ros_utils.parse_ease_topic(bag, '/Esmacat_write_gantry_ease');
ros_utils.save_table(fullfile(tabDir, 'gantry_write_tbl'), gantry_write_tbl);

gantry_read_tbl  = ros_utils.parse_ease_topic(bag, '/Esmacat_read_gantry_ease');
ros_utils.save_table(fullfile(tabDir, 'gantry_read_tbl'), gantry_read_tbl);

% Testing strings (optional, used later for condition labeling if desired)
gantry_testing_tbl = ros_utils.parse_std_msgs_topic(bag, '/gantry_testing');
ros_utils.save_table(fullfile(tabDir, 'gantry_testing_tbl'), gantry_testing_tbl);

% ---------------------- DECODE COMMANDED TARGETS --------------------------
% Align decoded X/Y (m) to pose frame using offset from first pose row
gantry_move_tbl = ros_utils.decode_gantry_move_table(gantry_write_tbl, gantry_pose_tbl);
if isempty(gantry_move_tbl) || height(gantry_move_tbl)==0
    error('No decoded gantry move commands.');
end
ros_utils.save_table(fullfile(tabDir, 'gantry_move_tbl'), gantry_move_tbl);

% ---------------------- REMOVE FIRST MOVE EVENT ----------------------------
gantry_move_tbl(1,:)    = [];
gantry_testing_tbl(1,:) = [];

% -------------------- ANALYZE MOVE LATENCY/ACCURACY -----------------------
% NOTE: arrival tolerance is set inside ros_utils.analyze_gantry_move.
gantry_move_analysis_tbl = ros_utils.analyze_gantry_move(gantry_move_tbl, gantry_pose_tbl);
ros_utils.save_table(fullfile(tabDir, 'gantry_move_analysis_tbl'), gantry_move_analysis_tbl);

% Print brief summary
log_utils.printf('\n=========== Gantry Move Summary ===========\n');
nMoves = height(gantry_move_analysis_tbl);
nArrived = sum(~isnan(gantry_move_analysis_tbl.TimeArrived));
log_utils.printf('Moves: %d, Arrived within tol: %d\n', nMoves, nArrived);
if nArrived > 0
    lat = gantry_move_analysis_tbl.Latency(~isnan(gantry_move_analysis_tbl.Latency));
    acc = gantry_move_analysis_tbl.Accuracy(~isnan(gantry_move_analysis_tbl.Accuracy));
    log_utils.printf('Latency (s): mean=%.3f, sd=%.3f, min=%.3f, max=%.3f\n', mean(lat), std(lat), min(lat), max(lat));
    log_utils.printf('Accuracy (m): mean=%.3f, sd=%.3f, min=%.3f, max=%.3f\n', mean(acc), std(acc), min(acc), max(acc));
end

% --------------------------- FIG 1 AND 2 ------------------------------
% Cleaned overlay of Pose vs Target (X and Y) vs time
target_pose_time_series( ...
    posF1, plotColF1, ...
    gantry_pose_tbl, gantry_move_tbl, ...
    figDir, plotXLimF1);

% Plot mean±SD velocity profiles, return stats struct
vel_out = move_velocity_time_series(posF2, plotColF2, ...
    gantry_pose_tbl, gantry_move_analysis_tbl, ...
    figDir, plotXLimF2, fs_resample, smooth_win);

% --------------------------- FIG 3 ------------------------------
lat_all  = gantry_move_analysis_tbl.Latency(:);         % seconds
cond_all = local_assign_conditions(gantry_move_analysis_tbl, gantry_testing_tbl);

valid = isfinite(lat_all) & isfinite(cond_all);
lat   = lat_all(valid);
cond  = cond_all(valid);

% command indices on the x-axis (preserve original ordering)
idx_all = (1:height(gantry_move_analysis_tbl)).';
idx     = idx_all(valid);

% split into 4 groups (x then y)
x1 = idx(cond==1); y1 = lat(cond==1);
x2 = idx(cond==2); y2 = lat(cond==2);
x3 = idx(cond==3); y3 = lat(cond==3);
x4 = idx(cond==4); y4 = lat(cond==4);

F = figure('Color','w','Position',posF3);

% plot scatter and boxchart
ylim_both = [1.375, 1.53];
ax1 = gantry_latency_scatter_plot4(x1,y1,x2,y2,x3,y3,x4,y4, posAxLeft,  plotColF3, plotXLimF3, ylim_both);
ax2 = gantry_latency_boxchart_plot4(y1,   y2,   y3,   y4, posAxRight, plotColF3, [0.5 4.5],  ylim_both);

% keep identical Y across axes
ylim(ax1, ylim_both); ylim(ax2, ylim_both); ax2.YAxis.Visible = 'off';

print(F, fullfile(figDir, 'gantry_move_latency_scatter_boxchart.svg'), '-dsvg', '-r300');
log_utils.printf('[INFO] Saved figure: gantry_move_latency_scatter_boxchart\n');

% --------------------------- FIG 4: Accuracy by condition --------------------
gantry_accuracy_bar_plot( ...
    gantry_move_analysis_tbl, ...
    gantry_testing_tbl, ...
    posF4, plotColF4, figDir);

% -------------------- ANALYZE MOVE FEATURES / REPORT ----------------------
gantry_stats_out = gantry_compute_and_print_stats( ...
    gantry_pose_tbl, ...
    gantry_move_analysis_tbl, ...
    gantry_testing_tbl, ...
    plotXLimF2, fs_resample, smooth_win);

log_utils.close()
end % gantry_movement

%% =========================================================================
function target_pose_time_series(posF1, plotColF1, pose_tbl, move_tbl, figDir, plotXLimF1)
% TARGET_POSE_TIME_SERIES
%   Clean overlay of measured pose and commanded targets (X and Y) vs time.
%   Based directly on the debug example in test_omni_ros_utils, with:
%     * relative time axis (t0 = first command time)
%     * explicit axis limits (plotTimeRange)
%     * improved styling and labeling
%
% Inputs:
%   posF1         : figure position [L B W H]
%   plotColF1     : struct with fields poseX, targX, poseY, targY (RGB)
%   pose_tbl      : table from ros_utils.parse_pose_topic
%   move_tbl      : table from ros_utils.decode_gantry_move_table
%   figDir    : output directory for figures
%   plotTimeRange : [tmin tmax] in seconds (relative to first command)
%   figExt        : 'svg' or 'png'
%   figDPI        : integer DPI for raster export
%
% Output:
%   Saves figure as: figDir/gantry_target_pose_time_series.<ext>

% Relative time base
t0 = min(move_tbl.Time);

% Pose traces (continuous)
PoseTime = pose_tbl.Time - t0;
PoseX    = pose_tbl.X;
PoseY    = pose_tbl.Y;

% Target traces (stairs)
TargTime = move_tbl.Time - t0;
TargX    = move_tbl.X;
TargY    = move_tbl.Y;

% Window mask
in_win_pose = PoseTime >= plotXLimF1(1) & PoseTime <= plotXLimF1(2);

% Figure
figure('Color','w'); set(gcf, 'Position', posF1);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- X subplot ---
nexttile; hold on
plot(PoseTime(in_win_pose), PoseX(in_win_pose), '-', 'Color', plotColF1(1,:), 'LineWidth', 1.6);
stairs(TargTime, TargX, '--', 'Color', plotColF1(2,:), 'LineWidth', 1.6);
xlim(plotXLimF1);
ylim([0,1.1]);
ylabel('X (m)');
box off
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XTick', [], 'XTickLabel', []);

% --- Y subplot ---
nexttile; hold on
plot(PoseTime(in_win_pose), PoseY(in_win_pose), '-', 'Color', plotColF1(1,:), 'LineWidth', 1.5, 'DisplayName', 'Pose Y');
stairs(TargTime, TargY, '--', 'Color', plotColF1(2,:), 'LineWidth', 1.6, 'DisplayName', 'Target Y');
xlim(plotXLimF1);
ylim([0,1.1]);
xlabel('Time (s)'); ylabel('Y (m)');
box off
set(gca, 'FontName', 'Arial', 'FontSize', 15);

% Save
print(fullfile(figDir, 'gantry_target_pose_time_series.svg'), '-dsvg', '-r300');
log_utils.printf('[INFO] Saved figure: gantry_target_pose_time_series\n');
end

%% =========================================================================
function vel_out = move_velocity_time_series(posF2, plotColF2, ...
    pose_tbl, analysis_tbl, ...
    figDir, plotXLimF2, fs_resample, smooth_win)
% Averaged velocity profile across all moves, aligned to command time.
%   - Preprocess pose to a uniform grid
%   - Zero-phase low-pass on position before differentiation
%   - Centered derivative for velocity
%   - Extract equal-length windows around each TimeWrite
%   - Plot single mean (black) with SD patch (gray)
%
% Inputs:
%   posF2        : figure position [L B W H]
%   plotColF2    : 1x3 color matrix
%   pose_tbl     : table from ros_utils.parse_pose_topic
%   analysis_tbl : table from ros_utils.analyze_gantry_move (uses TimeWrite)
%   figDir   : output directory
%   plotXLimF2   : [tmin tmax] seconds relative to command time (e.g., [0 2])
%   fs_resample  : target Hz for uniform grid (will not upsample above native)
%   smooth_win   : movmean window on velocity (samples), 1 disables
%
% Output:
%   vel_out struct with fields:
%       fs_work, fc_hz, t_rel, mean_v, sd_v, count, segments

% ==================== preprocess pose (uniform grid) ======================
t_pose = pose_tbl.Time(:);
x = pose_tbl.X(:);
y = pose_tbl.Y(:);

% sort and collapse duplicate timestamps (average ties)
[tp_sorted, ord] = sort(t_pose, 'ascend');
x_sorted = x(ord);
y_sorted = y(ord);
[tpu, ~, ic] = unique(tp_sorted, 'stable');
xu = accumarray(ic, x_sorted, [], @mean);
yu = accumarray(ic, y_sorted, [], @mean);

if numel(tpu) < 3
    error('Not enough pose samples to compute velocity.');
end

% native sample rate estimate
dt_native = median(diff(tpu));
Fs_native = 1 / dt_native;

% choose working sample rate: do not upsample beyond native
Fs_work = min(fs_resample, Fs_native);
Fs_work = max(20, round(Fs_work));  % keep reasonable floor
dt = 1 / Fs_work;

% build uniform time grid spanning all pose samples
tq = (tpu(1):dt:tpu(end)).';

% linear interpolation to uniform grid (filtered next)
xq = interp1(tpu, xu, tq, 'linear', 'extrap');
yq = interp1(tpu, yu, tq, 'linear', 'extrap');

% ==================== low-pass position then derivative ===================
fc_hz = min(6, 0.45*Fs_work);     % ~6 Hz or 45% of Nyquist
[b,a] = butter(4, fc_hz/(Fs_work/2), 'low');
xf = filtfilt(b, a, xq);
yf = filtfilt(b, a, yq);

vx = gradient(xf, dt);
vy = gradient(yf, dt);
v  = hypot(vx, vy);

if exist('smooth_win','var') && isnumeric(smooth_win) && smooth_win > 1
    v = movmean(v, smooth_win);
end

% ==================== build aligned velocity windows ======================
t_rel = (plotXLimF2(1) : dt : plotXLimF2(2)).';
Nwin  = numel(t_rel);

nMoves = height(analysis_tbl);
segments = [];

for i = 1:nMoves
    tw = analysis_tbl.TimeWrite(i);

    t_start = tw + plotXLimF2(1);
    i0 = find(tq >= t_start, 1, 'first');
    if isempty(i0), continue; end
    i1 = i0 + Nwin - 1;
    if i1 > numel(tq), continue; end

    seg_v = v(i0:i1).';
    segments = [segments; seg_v]; %#ok<AGROW>
end

count = size(segments,1);
if count == 0
    warning('No complete velocity windows found within requested plotXLim.');
    mean_v = nan(1, Nwin);
    sd_v   = nan(1, Nwin);
else
    mean_v = mean(segments, 1, 'omitnan');
    sd_v   = std(segments,  0, 1, 'omitnan');
end

% ============================== plot ======================================
figure('Color','w'); set(gcf, 'Position', posF2); hold on
% gSD patch
if count > 0
    fill([t_rel; flipud(t_rel)], [mean_v(:)+sd_v(:); flipud(mean_v(:)-sd_v(:))], ...
        plotColF2, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
end
% Mean line
plot(t_rel, mean_v, '-', 'Color', plotColF2, 'LineWidth', 1.6);

xlim(plotXLimF2);
ylim([0, 0.5]);
xlabel('Time from command (s)');
ylabel('Velocity (m/s)');
set(gca, 'FontName','Arial', 'FontSize', 15); box off

% Save
print(fullfile(figDir, 'gantry_velocity_time_series.svg'), '-dsvg', '-r300');
log_utils.printf('[INFO] Saved figure: gantry_velocity_time_series\n');

% ============================ outputs =====================================
vel_out = struct();
vel_out.fs_work  = Fs_work;
vel_out.fc_hz    = fc_hz;
vel_out.t_rel    = t_rel;
vel_out.mean_v   = mean_v;
vel_out.sd_v     = sd_v;
vel_out.count    = count;
vel_out.segments = segments;
end

% =========================================================================
function ax = gantry_latency_scatter_plot4(x1,y1,x2,y2,x3,y3,x4,y4, pos, colors, xlimv, ylimv)
% Scatter plot of gantry latencies (4 conditions) vs command index
% colors = [c1; c2; c3; c4]  (rows 1..4)

ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');

% scatter points (preserve 1..4 ordering)
plot(ax, x1, y1, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:), 'LineWidth', 1);
plot(ax, x2, y2, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:), 'LineWidth', 1);
plot(ax, x3, y3, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', colors(3,:), 'LineWidth', 1);
plot(ax, x4, y4, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(4,:), 'MarkerEdgeColor', colors(4,:), 'LineWidth', 1);

% mean lines (per condition)
ym = mean(y1, 'omitnan'); plot(ax, xlimv, [ym ym], '--', 'Color', colors(1,:), 'LineWidth', 1);
ym = mean(y2, 'omitnan'); plot(ax, xlimv, [ym ym], '--', 'Color', colors(2,:), 'LineWidth', 1);
ym = mean(y3, 'omitnan'); plot(ax, xlimv, [ym ym], '--', 'Color', colors(3,:), 'LineWidth', 1);
ym = mean(y4, 'omitnan'); plot(ax, xlimv, [ym ym], '--', 'Color', colors(4,:), 'LineWidth', 1);

% axes
xlim(ax, xlimv);
ylim(ax, ylimv);

set(ax, 'FontName', 'Arial', 'FontSize', 15);

% margins (exactly as gate helper)
ti = get(ax,'TightInset');
li = max(ti, [0.06 0.06 0.04 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);

xlabel(ax, 'Move event');
ylabel(ax, 'Latency (s)');
end

% =========================================================================
function ax = gantry_latency_boxchart_plot4(y1,y2,y3,y4, pos, colors, xlimv, ylimv)
% Violin plot of gantry latency distribution for 4 conditions
% colors = [c1; c2; c3; c4]

ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');
set(ax, 'FontName', 'Arial', 'FontSize', 15);

axes(ax);

% Plot boxcar
if ~isempty(y1), boxchart(1*ones(size(y1)), y1, 'BoxFaceColor', colors(1,:), 'BoxEdgeColor', colors(1,:), 'MarkerColor', colors(1,:)); end
if ~isempty(y2), boxchart(2*ones(size(y2)), y2, 'BoxFaceColor', colors(2,:), 'BoxEdgeColor', colors(2,:), 'MarkerColor', colors(2,:)); end
if ~isempty(y3), boxchart(3*ones(size(y3)), y3, 'BoxFaceColor', colors(3,:), 'BoxEdgeColor', colors(3,:), 'MarkerColor', colors(3,:)); end
if ~isempty(y4), boxchart(4*ones(size(y4)), y4, 'BoxFaceColor', colors(4,:), 'BoxEdgeColor', colors(4,:), 'MarkerColor', colors(4,:)); end

xlim(ax, xlimv);
ylim(ax, ylimv);

% margins (unchanged from gate helper)
ti = get(ax,'TightInset');
li = max(ti, [0.03 0.06 0.06 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);

xlabel(ax, 'Move condition');
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'1','2','3','4'};
ax.YAxis.Visible = 'off';
end

% ==========================================================================
function F = gantry_accuracy_bar_plot(analysis_tbl, testing_tbl, posF4, plotColF4, figDir)
% Color-coded bar plot of accuracy (cm) by move condition with SD error bars.
% - Conditions from local_assign_conditions(analysis_tbl, testing_tbl)
% - Bars: mean accuracy per condition
% - Error bars: SD per condition
% - Colors: rows of plotColF4 (1..4)

% ---- Extract accuracy and condition labels ----
acc_all = analysis_tbl.Accuracy(:) * 100;      % convert m → cm
cond_idx = local_assign_conditions(analysis_tbl, testing_tbl);  % 1..4

valid = isfinite(acc_all) & isfinite(cond_idx);
acc = acc_all(valid);
cond = cond_idx(valid);

% ---- Aggregate by condition ----
C = 4;
means = nan(1,C);
stds  = nan(1,C);
for c = 1:C
    v = acc(cond == c);
    if ~isempty(v)
        means(c) = mean(v, 'omitnan');
        stds(c)  = std(v, 0, 'omitnan');
    end
end

% ---- Build plot ----
F = figure('Color','w'); hold on;
set(F, 'Position', posF4);

x_pos   = 1:C;
barWidth = 0.6;

for i = 1:C
    if ~isnan(means(i))
        bar(x_pos(i), means(i), barWidth, 'FaceColor', plotColF4(i,:));
    else
        % draw a transparent placeholder if no data (keeps spacing)
        bar(x_pos(i), 0, barWidth, 'FaceColor', plotColF4(i,:), 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
    end
end

% Error bars (SD)
eb_means = means; eb_stds = stds;
eb_means(~isfinite(eb_means)) = 0;
eb_stds(~isfinite(eb_stds))   = 0;
errorbar(x_pos, eb_means, eb_stds, 'k.', 'CapSize', 6, 'LineWidth', 1.5);

% Axes/labels
set(gca, 'XTick', x_pos, 'XTickLabel', {'1','2','3','4'});
set(gca, 'FontName', 'Arial', 'FontSize', 15);
ylabel('Accuracy (cm)');
xlabel('Move condition');

% Limits
ylim([0, 1.5]);
xlim([0.5, C + 0.5]);

% Save
filename = fullfile(figDir, 'gantry_accuracy_bar.svg');
print(F, filename, '-dsvg', '-r300');
log_utils.printf('[INFO] Saved accuracy bar plot: %s\n', filename);
end

% ==========================================================================
function out = gantry_compute_and_print_stats( ...
    pose_tbl, analysis_tbl, testing_tbl, ...
    plotXLimF2, fs_resample, smooth_win)
% Console reporting of gantry move metrics using the SAME preprocessing as your F2 plot.

% -------------------------------- counts base --------------------------------
N = height(analysis_tbl); % moves analyzed (first already removed upstream)
cond_idx = local_assign_conditions(analysis_tbl, testing_tbl); % 1..4

% ----------------------- distance at command time (m) ------------------------
[dist_vals, n_dist, n_dist_excl] = local_distance_at_command(analysis_tbl, pose_tbl);

% --------------------- velocity window segments (for features) ---------------
[t_rel, segments, incl_mask] = local_velocity_segments( ...
    pose_tbl, analysis_tbl, plotXLimF2, fs_resample, smooth_win);
K_vel = size(segments,1);               % number of included moves
vel_excluded = N - K_vel;

% Per-move velocity features for included moves
[peak_v, t_peak, mean_v] = local_velocity_features(t_rel, segments);  % m/s, s, m/s

% ----------------------------- arrival latency ------------------------------
lat_all_s = analysis_tbl.Latency(:);                % seconds
valid_arrival = isfinite(lat_all_s);
lat_vals_ms = 1000 * lat_all_s(valid_arrival);      % convert to milliseconds for reporting
K_arr = numel(lat_vals_ms);

% ------------------------------- accuracy (cm) ------------------------------
acc_all_m = analysis_tbl.Accuracy(:);               % meters (from analyze_gantry_move)
acc_all_cm = 100 .* acc_all_m;                      % centimeters
valid_acc = isfinite(acc_all_cm);
K_acc = sum(valid_acc);

% ---------------------- arrival velocity (instant, m/s) ----------------------
[tq, vq] = local_velocity_timeseries(pose_tbl, fs_resample, smooth_win);  % tq in s, vq in m/s
arr_v = nan(N,1);
if ~isempty(tq)
    t_arr = analysis_tbl.TimeArrived(:);
    ok = isfinite(t_arr);
    % nearest-bin index
    idx = discretize(t_arr(ok), [-inf; (tq(1:end-1)+tq(2:end))/2; inf]);
    idx = max(1, min(numel(vq), idx));
    arr_v(ok) = vq(idx);
end
arr_v = arr_v(valid_arrival); % align with latency count

% ================================ LOGGING ====================================
log_utils.printf('\n=========== Gantry Movement: Inclusion Counts ===========\n');
log_utils.printf('Moves analyzed (after drop-first): %d\n', N);
log_utils.printf('Velocity windows included: %d (excluded: %d)\n', K_vel, vel_excluded);
log_utils.printf('Arrival latency computed: %d (excluded: %d)\n', K_arr, N - K_arr);
log_utils.printf('Distance computed: %d (excluded: %d)\n', n_dist, n_dist_excl);
log_utils.printf('Accuracy computed: %d (excluded: %d)\n', K_acc, N - K_acc);

% ============================== AGGREGATION ==================================
% Build aligned vectors over N for condition slicing
vec_dist   = nan(N,1);  % m
vec_pkv    = nan(N,1);  % m/s
vec_tpk    = nan(N,1);  % s
vec_mnv    = nan(N,1);  % m/s
vec_lat    = nan(N,1);  % ms
vec_arrv   = nan(N,1);  % m/s
vec_acc_cm = nan(N,1);  % cm

% Distance at command time: all N (already NaN-marked where not computable)
vec_dist = dist_vals;

% Velocity-derived features: only where incl_mask
vec_pkv(incl_mask) = peak_v;
vec_tpk(incl_mask) = t_peak;
vec_mnv(incl_mask) = mean_v;

% Arrival-based: only where valid_arrival
vec_lat(valid_arrival)  = lat_vals_ms;
vec_arrv(valid_arrival) = arr_v;

% Accuracy (cm): all N
vec_acc_cm(valid_acc) = acc_all_cm(valid_acc);

% ------------------------- Metric print blocks -------------------------------
% Each local_print_block prints overall and per-condition stats and returns a summary struct.
log_utils.printf('\n================ Overall and By Condition ================\n');

summ_dist = local_print_block('Distance (m)',           vec_dist,   cond_idx, 'm',   [], []);
summ_pkv  = local_print_block('Peak v (m/s)',           vec_pkv,    cond_idx, 'm/s', [], []);
summ_tpk  = local_print_block('T_peak (s)',             vec_tpk,    cond_idx, 's',   [], []);
summ_mnv  = local_print_block('Mean v (m/s)',           vec_mnv,    cond_idx, 'm/s', [], []);
summ_lat  = local_print_block('Arrival latency (ms)',   vec_lat,    cond_idx, 'ms',  [25 50], vec_lat); % latency-only bands
summ_arrv = local_print_block('Arrival velocity (m/s)', vec_arrv,   cond_idx, 'm/s', [], []);
summ_acc  = local_print_block('Accuracy (cm)',          vec_acc_cm, cond_idx, 'cm',  [], []);           % no bands for accuracy

% -------------------------- Between-condition summary ------------------------
% For each metric: deltas vs overall, range of means and medians, and % range of mean
log_utils.printf('\n================== Between-condition Summary ==================\n');
local_print_between_summary('Distance (m)',           summ_dist);
local_print_between_summary('Peak v (m/s)',           summ_pkv);
local_print_between_summary('T_peak (s)',             summ_tpk);
local_print_between_summary('Mean v (m/s)',           summ_mnv);
local_print_between_summary('Arrival latency (ms)',   summ_lat);
local_print_between_summary('Arrival velocity (m/s)', summ_arrv);
local_print_between_summary('Accuracy (cm)',          summ_acc);   % accuracy between-condition

% --------------------------------- outputs -----------------------------------
out = struct();
out.counts = struct('moves_analyzed', N, ...
    'vel_included', K_vel, 'vel_excluded', vel_excluded, ...
    'arrival_computed', K_arr, 'arrival_excluded', N - K_arr, ...
    'distance_computed', n_dist, 'distance_excluded', n_dist_excl, ...
    'accuracy_computed', K_acc, 'accuracy_excluded', N - K_acc);

out.t_rel   = t_rel;
out.included_mask_velocity = incl_mask;

out.vectors = struct( ...
    'distance_m', vec_dist, ...
    'peak_v_ms',  vec_pkv, ...
    't_peak_s',   vec_tpk, ...
    'mean_v_ms',  vec_mnv, ...
    'latency_ms', vec_lat, ...
    'arrival_v_ms', vec_arrv, ...
    'accuracy_cm', vec_acc_cm, ...
    'cond_idx',   cond_idx);

out.summaries = struct( ...
    'distance', summ_dist, ...
    'peak_v',   summ_pkv, ...
    't_peak',   summ_tpk, ...
    'mean_v',   summ_mnv, ...
    'latency',  summ_lat, ...
    'arr_v',    summ_arrv, ...
    'accuracy', summ_acc);
end

% ==============================================================================
% PRINTING HELPERS
% ==============================================================================

function summ = local_print_block(label_str, vecN, cond_idx, unit_str, band_ms, lat_vec_ms)
% Prints overall and per-condition summaries for a single metric.
% Returns a struct used later for between-condition summary.

% Overall
overall = local_comp_stats(vecN);

% Per-condition
per_cond = repmat(struct('mean',NaN,'sd',NaN,'median',NaN,'iqr',NaN,'min',NaN,'max',NaN,'n',0), 1, 4);
means = nan(1,4); medians = nan(1,4);
for c = 1:4
    vc = vecN(cond_idx==c);
    sc = local_comp_stats(vc);
    per_cond(c) = sc;
    means(c)   = sc.mean;
    medians(c) = sc.median;
end

% Print
log_utils.printf('\n--- %s ---\n', label_str);
local_print_stats_line('Overall', overall, unit_str);
for c = 1:4
    local_print_stats_line(sprintf('Cond %d', c), per_cond(c), unit_str);
end

% Optional latency band percentages around overall median
bands = struct();
if ~isempty(band_ms) && ~isempty(lat_vec_ms) && overall.n > 0 && isfinite(overall.median)
    ref_med = overall.median;
    for k = 1:numel(band_ms)
        w = band_ms(k);
        key = sprintf('pm_%dms', w);
        bands.(key) = nan(1,4);
        for c = 1:4
            lv = lat_vec_ms(cond_idx==c);
            lv = lv(isfinite(lv));
            if isempty(lv)
                bands.(key)(c) = NaN;
            else
                pct = 100 * mean(abs(lv - ref_med) <= w);
                bands.(key)(c) = pct;
            end
        end
    end
    % Print band table
    keys = fieldnames(bands);
    for ki = 1:numel(keys)
        kname = keys{ki};
        vals = bands.(kname);
        log_utils.printf('%% within %s of overall median: C1=%.1f%%, C2=%.1f%%, C3=%.1f%%, C4=%.1f%%\n', ...
            strrep(kname,'pm_','±'), vals(1), vals(2), vals(3), vals(4));
    end
end

% Assemble summary
summ = struct();
summ.overall   = overall;
summ.per_cond  = per_cond;
summ.centers   = struct('means', means, 'medians', medians, ...
    'overall_mean', overall.mean, ...
    'overall_median', overall.median);
summ.bands     = bands;
end

% ==========================================================================
function st = local_comp_stats(v)
% Compute mean, sd, median, IQR, min, max, n on finite entries.
v = v(:); v = v(isfinite(v));
st = struct('mean',NaN,'sd',NaN,'median',NaN,'iqr',NaN,'min',NaN,'max',NaN,'n',0);
if isempty(v), return; end
st.mean   = mean(v);
st.sd     = std(v);
st.median = median(v);
q = quantile(v, [0.25 0.75]);
st.iqr    = q(2) - q(1);
st.min    = min(v);
st.max    = max(v);
st.n      = numel(v);
end

% ==========================================================================
function local_print_stats_line(prefix, st, unit_str)
% Print one compact line: prefix: mean ± SD [median; IQR] {min–max} (n=K)
if st.n < 1
    log_utils.printf('%-10s: n=0\n', prefix);
    return;
end

if st.n == 1
    % SD and IQR are not meaningful with n=1
    log_utils.printf('%-10s: mean=%.3f %s; median=%.3f %s; range=%.3f–%.3f %s (n=1)\n', ...
        prefix, st.mean, unit_str, st.median, unit_str, st.min, st.max, unit_str);
    return;
end

log_utils.printf('%-10s: %.3f ± %.3f %s; median=%.3f %s [IQR=%.3f %s]; min–max=%.3f–%.3f %s (n=%d)\n', ...
    prefix, st.mean, st.sd, unit_str, st.median, unit_str, st.iqr, unit_str, st.min, st.max, unit_str, st.n);
end

% ==========================================================================
function local_print_between_summary(label_str, summ)
% Print deltas vs overall and ranges across condition centers.
m_over = summ.centers.overall_mean;
med_over = summ.centers.overall_median;
means   = summ.centers.means;
medians = summ.centers.medians;

% Deltas vs overall
dmean = means   - m_over;
dmed  = medians - med_over;

% Ranges
rng_mean = max(means)   - min(means);
rng_median = max(medians) - min(medians);

% Percent range relative to overall mean (guard zero)
if isfinite(m_over) && m_over ~= 0
    rng_mean_pct = 100 * rng_mean / m_over;
else
    rng_mean_pct = NaN;
end

log_utils.printf('\n[%s] Between-condition deltas and ranges\n', label_str);
log_utils.printf('Δmean vs overall:   C1=%.3f, C2=%.3f, C3=%.3f, C4=%.3f\n', dmean(1), dmean(2), dmean(3), dmean(4));
log_utils.printf('Δmedian vs overall: C1=%.3f, C2=%.3f, C3=%.3f, C4=%.3f\n', dmed(1),  dmed(2),  dmed(3),  dmed(4));
log_utils.printf('Range of means: %.3f  (%.1f%% of overall mean)\n', rng_mean, rng_mean_pct);
log_utils.printf('Range of medians: %.3f\n', rng_median);
end


% ==============================================================================
% DATA PIPELINE HELPERS
% ==============================================================================

function [dist_vals, n_ok, n_excl] = local_distance_at_command(analysis_tbl, pose_tbl)
% Distance from pose at command time (TimeWrite) to commanded target (TargetX, TargetY).
tw = analysis_tbl.TimeWrite(:);
tx = analysis_tbl.TargetX(:);
ty = analysis_tbl.TargetY(:);

tp = pose_tbl.Time(:); px = pose_tbl.X(:); py = pose_tbl.Y(:);

[tp_u, iu] = unique(tp, 'stable');
px_u = px(iu); py_u = py(iu);

px_at = interp1(tp_u, px_u, tw, 'linear', 'extrap');
py_at = interp1(tp_u, py_u, tw, 'linear', 'extrap');

dist_vals = hypot(tx - px_at, ty - py_at);

finite_mask = isfinite(dist_vals);
n_ok   = sum(finite_mask);
n_excl = numel(dist_vals) - n_ok;

dist_vals(~finite_mask) = NaN;
end

% ==========================================================================
function [t_rel, segments, incl_mask] = local_velocity_segments( ...
    pose_tbl, analysis_tbl, plotXLimF2, fs_resample, smooth_win)
% Clone of the F2 pipeline (no plotting). Returns velocity windows around each command.
t_pose = pose_tbl.Time(:);
x = pose_tbl.X(:); y = pose_tbl.Y(:);

[tp_sorted, ord] = sort(t_pose, 'ascend');
x_sorted = x(ord); y_sorted = y(ord);
[tpu, ~, ic] = unique(tp_sorted, 'stable');
xu = accumarray(ic, x_sorted, [], @mean);
yu = accumarray(ic, y_sorted, [], @mean);

if numel(tpu) < 3
    segments = []; t_rel = []; incl_mask = false(height(analysis_tbl),1);
    return;
end

dt_native = median(diff(tpu));
Fs_native = 1 / dt_native;
Fs_work = min(fs_resample, Fs_native);
Fs_work = max(20, round(Fs_work));
dt = 1 / Fs_work;

tq = (tpu(1):dt:tpu(end)).';
xq = interp1(tpu, xu, tq, 'linear', 'extrap');
yq = interp1(tpu, yu, tq, 'linear', 'extrap');

fc_hz = min(6, 0.45*Fs_work);
[b,a] = butter(4, fc_hz/(Fs_work/2), 'low');
xf = filtfilt(b, a, xq);
yf = filtfilt(b, a, yq);
vx = gradient(xf, dt);
vy = gradient(yf, dt);
v  = hypot(vx, vy);
if exist('smooth_win','var') && isnumeric(smooth_win) && smooth_win > 1
    v = movmean(v, smooth_win);
end

t_rel = (plotXLimF2(1) : dt : plotXLimF2(2)).';
Nwin  = numel(t_rel);

N = height(analysis_tbl);
segments = zeros(0, Nwin);
incl_mask = false(N,1);

for i = 1:N
    tw = analysis_tbl.TimeWrite(i);
    t_start = tw + plotXLimF2(1);
    i0 = find(tq >= t_start, 1, 'first');
    if isempty(i0), continue; end
    i1 = i0 + Nwin - 1;
    if i1 > numel(tq), continue; end

    seg_v = v(i0:i1).';
    segments(end+1, :) = seg_v; %#ok<AGROW>
    incl_mask(i) = true;
end
end

% ==========================================================================
function [peak_v, t_peak, mean_v] = local_velocity_features(t_rel, segments)
% Peak speed, time to peak, and mean speed per included move.
if isempty(segments)
    peak_v = []; t_peak = []; mean_v = [];
    return;
end
[peak_v, idx] = max(segments, [], 2, 'omitnan');
t_peak = t_rel(idx);
mean_v = mean(segments, 2, 'omitnan');

peak_v = peak_v(:);
t_peak = t_peak(:);
mean_v = mean_v(:);
end

% ==========================================================================
function [tq, vq] = local_velocity_timeseries(pose_tbl, fs_resample, smooth_win)
% Uniform-grid speed timeline v(t) consistent with the F2 pipeline.
t_pose = pose_tbl.Time(:);
x = pose_tbl.X(:); y = pose_tbl.Y(:);

[tp_sorted, ord] = sort(t_pose, 'ascend');
x_sorted = x(ord); y_sorted = y(ord);
[tpu, ~, ic] = unique(tp_sorted, 'stable');
xu = accumarray(ic, x_sorted, [], @mean);
yu = accumarray(ic, y_sorted, [], @mean);

if numel(tpu) < 3
    tq = []; vq = [];
    return;
end

dt_native = median(diff(tpu));
Fs_native = 1 / dt_native;
Fs_work   = min(fs_resample, Fs_native);
Fs_work   = max(20, round(Fs_work));
dt        = 1 / Fs_work;

tq = (tpu(1):dt:tpu(end)).';
xq = interp1(tpu, xu, tq, 'linear', 'extrap');
yq = interp1(tpu, yu, tq, 'linear', 'extrap');

fc_hz = min(6, 0.45*Fs_work);
[b,a] = butter(4, fc_hz/(Fs_work/2), 'low');
xf = filtfilt(b,a,xq);
yf = filtfilt(b,a,yq);

vx = gradient(xf, dt);
vy = gradient(yf, dt);
vq = hypot(vx,vy);
if exist('smooth_win','var') && isnumeric(smooth_win) && smooth_win > 1
    vq = movmean(vq, smooth_win);
end
end

% ==========================================================================
function cond_idx = local_assign_conditions(analysis_tbl, testing_tbl)
% Prefer labels from /gantry_testing if present; fallback to 1..4 repeating.
N = height(analysis_tbl);
cond_idx = nan(N,1);
has_testing = ~isempty(testing_tbl) && any(contains(testing_tbl.Event, 'gantry_test_move_'));

if has_testing
    evt_mask  = contains(testing_tbl.Event, 'gantry_test_move_');
    evt_times = testing_tbl.Time(evt_mask);
    evt_labels= testing_tbl.Event(evt_mask);
    evt_cond  = nan(numel(evt_labels),1);
    for k = 1:numel(evt_labels)
        tok = regexp(evt_labels(k), 'gantry_test_move_(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            c = str2double(tok{1});
            if isfinite(c) && c>=1 && c<=4
                evt_cond(k) = c;
            end
        end
    end
    valid_evt = isfinite(evt_cond);
    evt_times = evt_times(valid_evt);
    evt_cond  = evt_cond(valid_evt);

    for i = 1:N
        tw = analysis_tbl.TimeWrite(i);
        dt = abs(evt_times - tw);
        in_win = (evt_times >= (tw - 2.0)) & (evt_times <= (tw + 0.5));
        if any(in_win)
            [~, jlocal] = min(dt + (~in_win)*1e6);
            cond_idx(i) = evt_cond(jlocal);
        end
    end
end

fallback = isnan(cond_idx);
cond_idx(fallback) = mod(find(fallback)-1,4) + 1;
end

