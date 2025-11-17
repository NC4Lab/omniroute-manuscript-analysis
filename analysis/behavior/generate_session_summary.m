function generate_session_summary()
% GENERATE_SESSION_SUMMARY
% Build a meta-level session summary across all rats/sessions.
% Does not need to be run more than once.
%
% Output columns:
%   Rat, Phase, Session_Number, Date, Duration_Min, N_Trials,
%   N_Start_Cham_1, N_Start_Cham_3, N_Start_Cham_5, N_Start_Cham_7,
%   N_RightCue, N_LeftCue, N_FloorCue, N_BlackCue,
%   Performance_All, Performance_RightCue, Performance_LeftCue,
%   Performance_FloorCue, Performance_BlackCue, Data_Dir
%
% Dependencies/assumptions:
%   - Per-session CSV: trial_times.csv with fields:
%       StartTime, EndTrialTime, EndTime, Result, ChamberNumber,
%       LeftCue, RightCue, FloorCue (cue cols may be absent in some sessions)
%   - Phase cutoffs per rat in Excel at:
%       C:\Users\lester\OneDrive - UBC\lester2025omniroute\analysis\rat_behavior\io\<rat>_dates.xlsx
%     with headers (as in your file):
%       'P1 start date', 'P2 start date ', 'P3 start date ', 'Surgery',
%       'Training start after the surgery', 'First day of data Collection'
%

% ------------------- PATH & UTILS SET UP (edit here) ---------------------
% Configure path utils
utils_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
addpath(utils_dir);

% Path config:
procDir = fullfile(path_utils.dataset_root(), 'processed_behavior');

% ------------------------- USER SETTINGS ---------------------------------
rats    = [23 24 25]; % rat numbers

% ---- Load phase boundaries per rat (robust to header variants) ----
phaseInfo = containers.Map;
for r = rats
    xlsx = fullfile(procDir, sprintf('%d_dates.xlsx', r));
    if ~isfile(xlsx)
        error('Missing phase file: %s', xlsx);
    end
    T = readtable(xlsx, 'PreserveVariableNames', true); % keep headers as-is
    varMap = buildVarMap(T.Properties.VariableNames);   % normalized -> actual
    
    phases = struct();
    phases.P1 = getPhaseCutoff(T, varMap, 'p1startdate');
    phases.P2 = getPhaseCutoff(T, varMap, 'p2startdate');
    phases.P3 = getPhaseCutoff(T, varMap, 'p3startdate');
    % Surgery exists in your file, but we do NOT exclude here unless you later ask
    % phases.Surgery = getPhaseCutoff(T, varMap, 'surgery');
    
    phaseInfo(num2str(r)) = phases;
end

% ---- Walk sessions and build rows ----
rows = {};
for r = rats
    ratDir = fullfile(procDir, sprintf('NC4%04d', r));
    if ~isfolder(ratDir), warning('Missing rat dir: %s', ratDir); continue; end
    
    d = dir(ratDir);
    d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));
    if isempty(d), continue; end
    
    % sort sessions by folder name (assume YYMMDD)
    [~, idx] = sort({d.name});
    d = d(idx);
    
    phases = phaseInfo(num2str(r));
    
    for s = 1:numel(d)
        sessName = d(s).name;           % 'YYMMDD'
        sessDir  = fullfile(ratDir, sessName);
        trialCSV = fullfile(sessDir, 'trial_times.csv');
        if ~isfile(trialCSV), continue; end
        
        T = readtable(trialCSV, 'PreserveVariableNames', true);
        if isempty(T), continue; end
        
        % Keep only complete trials
        if any(strcmpi('Result', T.Properties.VariableNames))
            valid = ismember(string(T.Result), ["SUCCESS","ERROR"]);
            T = T(valid,:);
        else
            % No Result column -> skip session
            warning('No Result column in %s, skipping.', trialCSV);
            continue;
        end
        if isempty(T), continue; end
        
        % Metadata
        Rat            = r;
        Date           = sessName;      % keep YYMMDD string
        Session_Number = s;
        Data_Dir       = sessDir;
        
        % Phase (P1/P2/P3) from cutoffs
        Phase = assignPhaseFromCutoffs(sessName, phases);
        
        % Duration_Min: robust to EndTrialTime missing
        Duration_Min = computeDurationMinutes(T);
        
        % Counts (start chambers)
        N_Trials       = height(T);
        N_Start_Cham_1 = countEqual(T, 'ChamberNumber', 1);
        N_Start_Cham_3 = countEqual(T, 'ChamberNumber', 3);
        N_Start_Cham_5 = countEqual(T, 'ChamberNumber', 5);
        N_Start_Cham_7 = countEqual(T, 'ChamberNumber', 7);
        
        % ----- Cue counts & masks (match actual CSV values) -----
        n = height(T);

        % Right/Left goal cue: count when the triangle is on that side
        if any(strcmp('RightCue', T.Properties.VariableNames))
            valsR = string(T.RightCue);
            maskRight = strcmpi(valsR, 'Triangle');
        else
            maskRight = false(n,1);
        end

        if any(strcmp('LeftCue', T.Properties.VariableNames))
            valsL = string(T.LeftCue);
            maskLeft = strcmpi(valsL, 'Triangle');
        else
            maskLeft = false(n,1);
        end

        % Floor cue: 'Green' means ON; OFF is absence (NaN/empty)
        if any(strcmp('FloorCue', T.Properties.VariableNames))
            valsF = string(T.FloorCue);
            maskFloor = strcmpi(valsF, 'Green');
            maskBlack = ~maskFloor;                 % everything else = dark
            maskBlack(ismissing(valsF)) = true;     % explicitly include missing
        else
            maskFloor = false(n,1);
            maskBlack = true(n,1);                  % if no column, treat all as dark
        end

        % Counts
        N_RightCue = sum(maskRight);
        N_LeftCue  = sum(maskLeft);
        N_FloorCue = sum(maskFloor);
        N_BlackCue = sum(maskBlack);

        % Performance
        Performance_All        = mean(strcmp(T.Result, 'SUCCESS'));
        Performance_RightCue   = perfMasked(T, maskRight);
        Performance_LeftCue    = perfMasked(T, maskLeft);
        Performance_FloorCue   = perfMasked(T, maskFloor);
        Performance_BlackCue   = perfMasked(T, maskBlack);
        
        % Append row
        rows(end+1,:) = {Rat, Phase, Session_Number, Date, Duration_Min, ...
            N_Trials, N_Start_Cham_1, N_Start_Cham_3, N_Start_Cham_5, N_Start_Cham_7, ...
            N_RightCue, N_LeftCue, N_FloorCue, N_BlackCue, ...
            Performance_All, Performance_RightCue, Performance_LeftCue, ...
            Performance_FloorCue, Performance_BlackCue, Data_Dir}; %#ok<AGROW>
    end
end

% Assemble and save
colNames = {'Rat','Phase','Session_Number','Date','Duration_Min', ...
    'N_Trials','N_Start_Cham_1','N_Start_Cham_3','N_Start_Cham_5','N_Start_Cham_7', ...
    'N_RightCue','N_LeftCue','N_FloorCue','N_BlackCue', ...
    'Performance_All','Performance_RightCue','Performance_LeftCue', ...
    'Performance_FloorCue','Performance_BlackCue','Data_Dir'};

if isempty(rows)
    warning('No rows generated; check inputs and filters.');
    outTable = cell2table(cell(0, numel(colNames)), 'VariableNames', colNames);
else
    outTable = cell2table(rows, 'VariableNames', colNames);
end

outCSV = fullfile(procDir, 'session_summary.csv');
writetable(outTable, outCSV);
fprintf('Saved session summary to %s\n', outCSV);

end % generate_session_summary


% ======================= Helpers =======================

function varMap = buildVarMap(varNames)
% Map normalized names -> actual names (first occurrence wins).
% Normalization: lowercase, strip all non-alnum.
varMap = containers.Map;
for i = 1:numel(varNames)
    key = normalizeKey(varNames{i});
    if ~isKey(varMap, key)
        varMap(key) = varNames{i};
    end
end
end

function key = normalizeKey(s)
s = char(s);
s = lower(s);
key = regexprep(s, '[^a-z0-9]', '');
end

function cutoff = getPhaseCutoff(T, varMap, canonicalKey)
% canonicalKey: 'p1startdate' | 'p2startdate' | 'p3startdate' | 'surgery'
if ~isKey(varMap, canonicalKey)
    cutoff = NaN;
    return;
end
col = varMap(canonicalKey);
val = T.(col)(1);
cutoff = normalizeYYMMDD(val); % numeric YYMMDD (e.g., 250403), or NaN
end

function yymmdd = normalizeYYMMDD(val)
% Return numeric YYMMDD (e.g., 250403). NaN if cannot parse.

% datetime directly
if isdatetime(val)
    yymmdd = str2double(datestr(val, 'yymmdd'));
    return;
end

% string/char
if isstring(val) || ischar(val)
    s = regexprep(char(val), '\s+', '');
    % digits only?
    if all(isstrprop(s, 'digit'))
        if numel(s) == 6
            yymmdd = str2double(s);
            return;
        elseif numel(s) == 8 % yyyymmdd -> use last 6
            yymmdd = str2double(s(3:8));
            return;
        end
    end
    % try common date formats
    fmtList = {'yyyy-MM-dd','yyyy/MM/dd','MM/dd/yyyy','dd/MM/yyyy','yyMMdd','yyyyMMdd'};
    for k = 1:numel(fmtList)
        try
            dt = datetime(s, 'InputFormat', fmtList{k});
            yymmdd = str2double(datestr(dt, 'yymmdd'));
            return;
        catch
        end
    end
    % last resort: extract digits and interpret
    digitsOnly = regexprep(s, '[^\d]', '');
    if numel(digitsOnly) >= 6
        yymmdd = str2double(digitsOnly(end-5:end));
        return;
    end
    yymmdd = NaN;
    return;
end

% numeric
if isnumeric(val) && isscalar(val) && ~isnan(val)
    v = double(val);
    if v >= 1e5 && v <= 991231
        % already YYMMDD
        yymmdd = v;
        return;
    elseif v >= 19000101 && v <= 20991231
        % YYYYMMDD -> YYMMDD
        s = sprintf('%.0f', v);
        yymmdd = str2double(s(3:8));
        return;
    elseif v > 20000 && v < 80000
        % likely Excel serial date (days)
        try
            dt = datetime(v, 'ConvertFrom', 'excel');
            yymmdd = str2double(datestr(dt, 'yymmdd'));
            return;
        catch
        end
    end
end

yymmdd = NaN;
end

function Phase = assignPhaseFromCutoffs(sessName, phases)
sess = normalizeYYMMDD(sessName);
if isnan(sess)
    Phase = NaN;   % Unknown phase
    return;
end
p1 = phases.P1; p2 = phases.P2; p3 = phases.P3;
if ~isnan(p3) && sess >= p3
    Phase = 3;
elseif ~isnan(p2) && sess >= p2
    Phase = 2;
elseif ~isnan(p1)
    Phase = 1;
else
    Phase = NaN;   % Unknown
end
end

function Duration_Min = computeDurationMinutes(T)
% Timestamps from extractor are ROS nanoseconds; handle both ns and s.
getcol = @(name) getColIfExists(T, name);

startCol = getcol('StartTime');
endtrialCol = getcol('EndTrialTime');
endCol = getcol('EndTime');

if isempty(startCol), Duration_Min = NaN; return; end
t0 = min(startCol(~isnan(startCol)));
if ~isempty(endtrialCol) && any(~isnan(endtrialCol))
    t1 = max(endtrialCol(~isnan(endtrialCol)));
elseif ~isempty(endCol) && any(~isnan(endCol))
    t1 = max(endCol(~isnan(endCol)));
else
    Duration_Min = NaN; return;
end

dt = double(t1) - double(t0);
% Heuristic: if values are in nanoseconds, dt will be huge
if dt > 1e7  % > 10 million seconds threshold is too big; better: check raw scale
    % Better: detect unit by magnitude of absolute timestamps
    % If absolute timestamps > 1e10, they are likely ns
    if (abs(double(t1)) > 1e10) || (abs(double(t0)) > 1e10)
        Duration_Min = dt / 1e9 / 60;
    else
        Duration_Min = dt / 60; % assume seconds
    end
else
    % if dt small, assume seconds
    Duration_Min = dt / 60;
end
end

function col = getColIfExists(T, varName)
if any(strcmp(varName, T.Properties.VariableNames))
    col = T.(varName);
else
    col = [];
end
end

function n = countEqual(T, varName, targetVal)
if any(strcmp(varName, T.Properties.VariableNames))
    v = T.(varName);
    n = sum(v == targetVal);
else
    n = 0;
end
end

function mask = cueMask(T, varName, trueTokens)
% Return logical vector for cue == any(trueTokens)
% Missing column -> all-false vector.
trueTokens = lower(string(trueTokens));
if ~any(strcmp(varName, T.Properties.VariableNames))
    mask = false(height(T),1);
    return;
end
vals = lower(string(T.(varName)));
mask = false(size(vals));
for k = 1:numel(trueTokens)
    mask = mask | (vals == trueTokens(k));
end
end

function p = perfMasked(T, mask)
if ~any(mask)
    p = NaN;
    return;
end
p = mean(strcmp(T.Result(mask), 'SUCCESS'));
end
