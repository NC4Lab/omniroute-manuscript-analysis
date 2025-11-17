function system_latency_analysis()
% =========================================================================
% FUNCTION: system_latency_analysis
%
% Description:
%   Compute and plot projection, audio, and gate movement latencies using
%   ros_utils class methods.
% =========================================================================

clc; close all;

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Paths
out_dir_name = 'system_latency_test';
bagPathProj = fullfile(path_utils.dataset_root(), out_dir_name, 'ros_session_20250815_214051.bag');
bagPathGate = fullfile(path_utils.dataset_root(), out_dir_name, 'ros_session_20250815_210959.bag');
figDir = path_utils.results_figures(out_dir_name);
tabDir = path_utils.results_tables(out_dir_name);

% Logging
logFile = path_utils.results_log(out_dir_name);
log_utils.start(logFile);

% ---------------------- USER SETTINGS (edit here) -------------------------

% Topics
topic_projection_image = '/projection_image';
topic_rosout           = '/rosout';

% X-axis plotting range
plotXLimF1 = [0 20];   % Projection
plotXLimF2 = [0 20];   % Audio
plotXLimF3 = [0 40];   % Gate

% Y-axis plotting range
plotYLimF1 = [320 400]; % Projection
plotYLimF2 = [700 750]; % Audio
plotYLimF3 = [850 1000];  % Gate (ms) - adjust as needed

% Figure colors
plotColF1 = [0.35 0.70 0.30]; % Projection
plotColF2 = [0.60 0.40 0.80]; % Audio
plotColF3 = [0.44 0.76 1.00;   % Up
    0.14 0.42 0.70];  % Down

% Figure layouts
posF1 = [100, 100, 425, 300];
posF2 = [600, 100, 425, 300];
posF3 = [1100, 100, 575, 300];

% Left and right axis shared for all plots
posAxLeftF1  = [0.06 0.1 0.75 0.9];
posAxRightF1 = [0.75 0.1 0.25 0.9];

posAxLeftF2  = [0.06 0.1 0.75 0.9];
posAxRightF2 = [0.75 0.1 0.25 0.9];

posAxLeftF3  = [0.06 0.1 0.7 0.9];
posAxRightF3 = [0.7 0.1 0.3 0.9];

%% -------------------------- PROJECTION -----------------------------------
log_utils.printf('[INFO] Loading projection ROS bag: %s\n', bagPathProj);
bag_proj = rosbag(bagPathProj);

image_cmd_tbl = ros_utils.parse_int32_multiarray(bag_proj, topic_projection_image);
ros_utils.save_table(fullfile(tabDir, 'image_cmd_tbl'), image_cmd_tbl);

sound_cmd_tbl = ros_utils.parse_std_msgs_topic(bag_proj, '/sound_cmd');
ros_utils.save_table(fullfile(tabDir, 'sound_cmd_tbl'), sound_cmd_tbl);

rosout_tbl = ros_utils.parse_rosout(bag_proj, topic_rosout);
ros_utils.save_table(fullfile(tabDir, 'rosout_tbl'), rosout_tbl);

proj_latency_tbl = ros_utils.compute_projection_last_render_latency(image_cmd_tbl, rosout_tbl);
ros_utils.save_table(fullfile(tabDir, 'proj_latency_tbl'), proj_latency_tbl);

audio_latency_tbl = ros_utils.compute_audio_start_latency(sound_cmd_tbl, rosout_tbl);
ros_utils.save_table(fullfile(tabDir, 'audio_latency_tbl'), audio_latency_tbl);

%% ---------------------------- GATE ---------------------------------------
log_utils.printf('[INFO] Loading gate ROS bag: %s\n', bagPathGate);
bag_gate = rosbag(bagPathGate);

maze_write_tbl   = ros_utils.parse_ease_topic(bag_gate, '/Esmacat_write_maze_ease');
maze_read_tbl    = ros_utils.parse_ease_topic(bag_gate, '/Esmacat_read_maze_ease');
gate_testing_tbl = ros_utils.parse_std_msgs_topic(bag_gate, '/gate_testing');

gate_latency_tbl = ros_utils.compute_gate_move_latency(maze_write_tbl, maze_read_tbl, gate_testing_tbl);
ros_utils.save_table(fullfile(tabDir, 'gate_latency_tbl'), gate_latency_tbl);

%% --------------------------- SUMMARIES -----------------------------------
system_compute_and_print_stats(proj_latency_tbl, audio_latency_tbl, gate_latency_tbl);

%% --------------------------- PLOTS ---------------------------------------
lat_img_s = proj_latency_tbl.LatencyLast_s(isfinite(proj_latency_tbl.LatencyLast_s));
lat_aud_s = audio_latency_tbl.LatencyStart_s(isfinite(audio_latency_tbl.LatencyStart_s));
lat_gate_s = gate_latency_tbl.Latency_s(isfinite(gate_latency_tbl.Latency_s));

% Projection plot
idx_img    = (1:numel(lat_img_s)).';
lat_img_ms = lat_img_s * 1000;
F = figure('Color','w','Position',posF1);
proj_latency_scatter_plot(idx_img, lat_img_ms, posAxLeftF1,  plotColF1, plotXLimF1, plotYLimF1, 'Projection event');
proj_latency_boxchart_plot( lat_img_ms, posAxRightF1, plotColF1, [0.5 1.5], plotYLimF1);
print(F, fullfile(figDir, 'proj_latency_scatter_boxchart.svg'), '-dsvg', '-r300');

% Audio plot
idx_aud    = (1:numel(lat_aud_s)).';
lat_aud_ms = lat_aud_s * 1000;
F2 = figure('Color','w','Position',posF2);
proj_latency_scatter_plot(idx_aud, lat_aud_ms, posAxLeftF2,  plotColF2, plotXLimF2, plotYLimF2, 'Sound event');
proj_latency_boxchart_plot( lat_aud_ms, posAxRightF2, plotColF2, [0.5 1.5], plotYLimF2);
print(F2, fullfile(figDir, 'audio_latency_scatter_boxchart.svg'), '-dsvg', '-r300');

% Gate plot (Up vs Down)
idx_up    = find(gate_latency_tbl.MoveDir == "up");
idx_down  = find(gate_latency_tbl.MoveDir == "down");

lat_up_ms   = gate_latency_tbl.Latency_s(idx_up)   * 1000;
lat_down_ms = gate_latency_tbl.Latency_s(idx_down) * 1000;

F3 = figure('Color','w','Position',posF3);
gate_latency_scatter_plot(idx_up, lat_up_ms, idx_down, lat_down_ms, ...
    posAxLeftF3, plotColF3, plotXLimF3, plotYLimF3);
gate_latency_boxchart_plot(lat_up_ms, lat_down_ms, ...
    posAxRightF3, plotColF3, [0.5 2.5], plotYLimF3);

print(F3, fullfile(figDir, 'gate_latency_scatter_boxchart.svg'), '-dsvg', '-r300');

log_utils.printf('\n[DONE] Tables saved to: %s\n', tabDir);

log_utils.close()
end % system_latency_analysis

% =========================================================================
function ax = proj_latency_scatter_plot(x, y, pos, color, xlimv, ylimv, xLabel)
ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');
plot(ax, x, y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'LineWidth', 1);
ymean = mean(y, 'omitnan');
plot(ax, xlimv, [ymean ymean], '--', 'Color', color, 'LineWidth', 1);
if isinf(xlimv(2)), xlimv(2) = max(x); end
xlim(ax, xlimv); 
ylim(ax, ylimv);
set(ax, 'FontName', 'Arial', 'FontSize', 15);
ti = get(ax,'TightInset');
li = max(ti, [0.06 0.06 0.04 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);
xlabel(ax, xLabel); ylabel(ax, 'Latency (ms)');
end

% =========================================================================
function ax = proj_latency_boxchart_plot(y, pos, color, xlimv, ylimv)
ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');
set(ax, 'FontName', 'Arial', 'FontSize', 15);
axes(ax);
boxchart(ones(size(y)), y, ...
    'BoxFaceColor', color, 'BoxEdgeColor', color, ...
    'WhiskerLineColor', [0 0 0], 'MarkerColor', color);
xlim(ax, xlimv); 
ylim(ax, ylimv);
ti = get(ax,'TightInset');
li = max(ti, [0.03 0.06 0.06 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);
ax.XTick = []; ax.YAxis.Visible = 'off';
end

% =========================================================================
function ax = gate_latency_scatter_plot(x_up, y_up, x_down, y_down, pos, colors, xlimv, ylimv)
% Helper: scatter plot of gate latencies (up vs down) vs index
% colors = [up_color; down_color]

ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');

% scatter points
plot(ax, x_up,   y_up,   'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:), 'LineWidth', 1);
plot(ax, x_down, y_down, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:), 'LineWidth', 1);

% mean lines
ymean_up = mean(y_up, 'omitnan');
plot(ax, xlimv, [ymean_up ymean_up], '--', 'Color', colors(1,:), 'LineWidth', 1);
ymean_down = mean(y_down, 'omitnan');
plot(ax, xlimv, [ymean_down ymean_down], '--', 'Color', colors(2,:), 'LineWidth', 1);

% axes
if isinf(xlimv(2)), xlimv(2) = max([x_up; x_down]); end
xlim(ax, xlimv);
ylim(ax, ylimv);

set(ax, 'FontName', 'Arial', 'FontSize', 15);

% margins (same style as latency_scatter_plot)
ti = get(ax,'TightInset');
li = max(ti, [0.06 0.06 0.04 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);

xlabel(ax, 'Move event');
ylabel(ax, 'Latency (ms)');
end

% =========================================================================
function ax = gate_latency_boxchart_plot(y_up, y_down, pos, colors, xlimv, ylimv)
% Helper: boxchart plot of gate latency distribution (up vs down)
% colors = [up_color; down_color]

ax = axes('Position', pos);
hold(ax, 'on'); box(ax, 'off'); set(ax, 'Color','w');
set(ax, 'FontName', 'Arial', 'FontSize', 15);

axes(ax);

% Plot boxcar
boxchart(ones(size(y_up)),   y_up,   'BoxFaceColor', colors(1,:), 'BoxEdgeColor', colors(1,:), 'MarkerColor', colors(1,:));
boxchart(2*ones(size(y_down)), y_down, 'BoxFaceColor', colors(2,:), 'BoxEdgeColor', colors(2,:), 'MarkerColor', colors(2,:));

xlim(ax, xlimv);
ylim(ax, ylimv);

% margins (similar to latency_boxchart_plot)
ti = get(ax,'TightInset');
li = max(ti, [0.03 0.06 0.06 0.08]);
ax.Position = [pos(1)+li(1) pos(2)+li(2) pos(3)-li(1)-li(3) pos(4)-li(2)-li(4)];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);

ax.XTick = [1 2];
ax.XTickLabel = {'Up', 'Down'};
ax.YAxis.Visible = 'off';
end

% =========================================================================
function out = system_compute_and_print_stats(proj_latency_tbl, audio_latency_tbl, gate_latency_tbl)
% SYSTEM_COMPUTE_AND_PRINT_STATS
% Prints summaries (mean ± SD, n) for projection, audio, gate latencies
% in milliseconds. Gate is reported overall and split by Up/Down.
%
% Expected columns:
%   proj_latency_tbl.LatencyLast_s
%   audio_latency_tbl.LatencyStart_s
%   gate_latency_tbl.Latency_s, gate_latency_tbl.MoveDir

% --------- Extract finite latencies (convert to ms) ---------
lat_proj  = 1000 * proj_latency_tbl.LatencyLast_s(isfinite(proj_latency_tbl.LatencyLast_s));
lat_audio = 1000 * audio_latency_tbl.LatencyStart_s(isfinite(audio_latency_tbl.LatencyStart_s));

lat_gate_all = 1000 * gate_latency_tbl.Latency_s(isfinite(gate_latency_tbl.Latency_s));

% robust parsing of MoveDir -> "up"/"down"
lat_gate_up = []; lat_gate_down = [];
if ismember('MoveDir', gate_latency_tbl.Properties.VariableNames)
    dir_raw = gate_latency_tbl.MoveDir;
    if iscategorical(dir_raw), dir_raw = string(dir_raw); end
    if iscellstr(dir_raw), dir_raw = string(dir_raw); end
    if isstring(dir_raw) || ischar(dir_raw)
        dir_str = lower(string(dir_raw));
        m_up    = isfinite(gate_latency_tbl.Latency_s) & (dir_str == "up");
        m_down  = isfinite(gate_latency_tbl.Latency_s) & (dir_str == "down");
        lat_gate_up   = 1000 * gate_latency_tbl.Latency_s(m_up);
        lat_gate_down = 1000 * gate_latency_tbl.Latency_s(m_down);
    end
end

% --------- Printed summaries (match other script style) ---------
log_utils.printf('\n=========== System Latency Summary ===========\n');
local_print_stat('Projection latency (ms)', lat_proj);
local_print_stat('Audio latency (ms)',      lat_audio);

% Gate overall + split
local_print_stat('Gate latency (ms)',        lat_gate_all);
local_print_stat('Gate latency — Up (ms)',   lat_gate_up);
local_print_stat('Gate latency — Down (ms)', lat_gate_down);

% --------- t-test for Up vs Down gate latency (paired) ---------
if ~isempty(lat_gate_up) && ~isempty(lat_gate_down)
    k = min(numel(lat_gate_up), numel(lat_gate_down));
    if k < 2
        log_utils.printf('\nUp vs Down gate latency paired t-test: n/a (need at least 2 paired samples)\n');
    else
        % Pair in recorded order (truncate to equal lengths if needed)
        x = lat_gate_up(1:k);
        y = lat_gate_down(1:k);

        % Paired, two-tailed t test
        [~, p_t, ci, stats] = ttest(x, y);

        % Cohen's d for paired samples (dz): mean(diff)/std(diff) == t/sqrt(n)
        d        = x - y;
        sd_diff  = std(d, 'omitnan');
        if sd_diff > 0 && isfinite(sd_diff)
            cohen_d = mean(d, 'omitnan') / sd_diff;
        else
            cohen_d = NaN;
        end

        log_utils.printf(['\nUp vs Down gate latency (paired two-tailed t test): ' ...
                 't(%d) = %.2f, p = %.4g, Cohen''s d = %.2f, 95%% CI [%.2f, %.2f] (n=%d pairs)\n'], ...
                 stats.df, stats.tstat, p_t, cohen_d, ci(1), ci(2), k);

        if numel(lat_gate_up) ~= numel(lat_gate_down)
            log_utils.printf('[NOTE] Paired test used the first %d up/down events to form pairs (truncated to equal lengths).\n', k);
        end
    end
else
    log_utils.printf('\nUp vs Down gate latency paired t-test: n/a (insufficient data)\n');
end

% --------- outputs (optional programmatic use) ---------------------------
out = struct();
out.projection = local_pack(lat_proj);
out.audio      = local_pack(lat_audio);
out.gate_all   = local_pack(lat_gate_all);
out.gate_up    = local_pack(lat_gate_up);
out.gate_down  = local_pack(lat_gate_down);
end

% ---------- local helpers ----------
% =========================================================================
function S = local_pack(v)
S = struct('n',numel(v), ...
    'mean',mean(v,'omitnan'), ...
    'sd',std(v,'omitnan'), ...
    'median',median(v,'omitnan'), ...
    'min',min(v,[],'omitnan'), ...
    'max',max(v,[],'omitnan'));
end

% =========================================================================
function local_print_stat(label_str, vec)
% Print "Label: mean ± SD (n=K)"; handles empty/NaN.
vec = vec(:); vec = vec(isfinite(vec));
if isempty(vec)
    log_utils.printf('%s: n=0\n', label_str);
else
    mu = mean(vec); sd = std(vec);
    log_utils.printf('%s: %.3f ± %.3f (n=%d)\n', label_str, mu, sd, numel(vec));
end
end


