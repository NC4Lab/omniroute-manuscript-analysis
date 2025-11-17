function inspect_ephys_channels()
% =========================================================================
% FUNCTION: inspect_ephys_channels
%
% Description:
%   Load CSC data and plot a specified time window for ALL channels using a
%   fixed y-axis range. Layout is 2 columns where:
%       Column 1 = channels 1..N/2 (top to bottom)
%       Column 2 = channels N/2+1..N (top to bottom)
%   For 64 channels: first column 1–32, second column 33–64.
% =========================================================================

% --- USER SETTINGS ---
matPath   = ''; % path to MAT file (e.g., '20250806_164955_csc_export.mat')

% Plot window relative to start of recording [t_start, t_end] in seconds
timeRange = [60, 120];                   % modify as needed
yLims     = [-1000, 1000];             % fixed ylim for ALL subplots (µV)

% --- LOAD DATA ---
mat     = load(matPath);
cscData = mat.csc_data * mat.gain_to_uv;   % convert to µV
cscTs   = mat.ros_ts;                      % seconds (absolute)

% --- SELECT TIME RANGE (relative to first sample) ---
t0_abs   = cscTs(1);
tStart   = t0_abs + timeRange(1);
tEnd     = t0_abs + timeRange(2);
idx      = (cscTs >= tStart) & (cscTs <= tEnd);

if ~any(idx)
    error('No samples found in requested timeRange [%.3f %.3f] s.', timeRange(1), timeRange(2));
end

t_rel       = cscTs(idx) - tStart;        % time relative to window start
x           = cscData(idx, :);
durationSec = diff(timeRange);

% --- PLOT CONFIG: 2 columns; first col = 1..nRow, second col = nRow+1..nCh ---
nCh  = size(x, 2);
nCol = 2;
nRow = ceil(nCh / nCol);                   % rows needed to stack first half

fig = figure('Position',[100 100 3000 1600]);
tl  = tiledlayout(nRow, nCol, 'TileSpacing','tight', 'Padding','compact');

for r = 1:nRow
    % Left column channel index
    chL = r;
    tileL = (r-1)*nCol + 1;               % tile index for left cell in row r
    if chL <= nCh
        nexttile(tl, tileL);
        plot(t_rel, x(:, chL), 'k', 'LineWidth', 0.5);
        ylim(yLims); xlim([0, durationSec]);
        ylabel(sprintf('Ch %d', chL));
        if r < nRow, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
    end

    % Right column channel index
    chR = r + nRow;
    tileR = (r-1)*nCol + 2;               % tile index for right cell in row r
    if chR <= nCh
        nexttile(tl, tileR);
        plot(t_rel, x(:, chR), 'k', 'LineWidth', 0.5);
        ylim(yLims); xlim([0, durationSec]);
        ylabel(sprintf('Ch %d', chR));
        if r < nRow, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
    end
end

title(tl, sprintf('CSC: t = [%g, %g] s; fixed ylim [%d, %d] µV', ...
    timeRange(1), timeRange(2), yLims(1), yLims(2)));

end
