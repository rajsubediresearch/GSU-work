% Select a file (CSV/XLSX/XLS), compute excess deaths and uncertainty bounds
% using:
%   Excess(t) = max(Observed - Predicted, 0)
%   LB(t)     = max(Observed - UpperPI, 0)
%   UB(t)     = max(Observed - LowerPI, 0)
%
% Compute totals ONLY for rows within YEAR_START:YEAR_END (default 2020–2023),
% and save:
%   1) Full output with new columns (same format as input extension): *_excess.*
%   2) Summary totals table labeled by window
%   3) (Optional) A filtered file containing ONLY rows within the window
%
% Expected typical column names:
%   time/year/date, data/observed, median/predicted, LB (LowerPI), UB (UpperPI)
%
% Author: (your name)
% Date: (today)

clear; clc;

% ----------------------------
% SETTINGS (edit here)
% ----------------------------
YEAR_START = 2020;
YEAR_END   = 2023;

% Save a second dataset containing only rows within YEAR_START:YEAR_END
SAVE_FILTERED_WINDOW_FILE = true;

% ----------------------------
% 1) Select file
% ----------------------------
[fname, fpath] = uigetfile( ...
    {'*.csv;*.xlsx;*.xls', 'Data files (*.csv, *.xlsx, *.xls)'; '*.*','All files'}, ...
    'Select Observed/Predicted file');

if isequal(fname,0)
    disp('No file selected. Exiting.');
    return;
end

infile = fullfile(fpath, fname);
[~, baseName, ext] = fileparts(infile);

% ----------------------------
% 2) Read file
% ----------------------------
T = readtable(infile, 'PreserveVariableNames', true);

% ----------------------------
% 3) Detect columns
% ----------------------------
rawNames  = T.Properties.VariableNames;
normNames = regexprep(lower(rawNames), '[^a-z0-9]', ''); % normalize for matching

getCol = @(cands) localGetCol(rawNames, normNames, cands);

timeCol  = getCol({'time','year','date','t'});
obsCol   = getCol({'data','observed','obs','y','count','deaths'});
predCol  = getCol({'median','predicted','pred','fit','expected','mean'});
lowerCol = getCol({'lb','lowerpi','lower','lpi','lowerbound','pilower','pi_lower'});
upperCol = getCol({'ub','upperpi','upper','upi','upperbound','piupper','pi_upper'});

% Validate required columns
missing = {};
if isempty(obsCol),   missing{end+1} = 'Observed';  end %#ok<SAGROW>
if isempty(predCol),  missing{end+1} = 'Predicted'; end %#ok<SAGROW>
if isempty(lowerCol), missing{end+1} = 'LowerPI';   end %#ok<SAGROW>
if isempty(upperCol), missing{end+1} = 'UpperPI';   end %#ok<SAGROW>

if ~isempty(missing)
    error(['Could not detect required columns: %s' newline ...
           'Your file columns are: %s' newline ...
           'Tip: rename columns to: time, data, median, LB, UB (or similar).'], ...
           strjoin(missing, ', '), strjoin(rawNames, ', '));
end

if isempty(timeCol)
    warning(['Time column not detected. Window totals (%d–%d) cannot be computed.' ...
             ' The script will still compute per-row Excess/LB/UB.'], YEAR_START, YEAR_END);
end

% Extract numeric arrays
Observed  = localToDouble(T.(obsCol));
Predicted = localToDouble(T.(predCol));
LowerPI   = localToDouble(T.(lowerCol));
UpperPI   = localToDouble(T.(upperCol));

% ----------------------------
% 4) Compute Excess(t), LB(t), UB(t)
% ----------------------------
Excess_t = max(Observed - Predicted, 0);

% According to your definitions:
% LB(t) = Observed - UpperPI if Observed > UpperPI else 0
LB_t = max(Observed - UpperPI, 0);

% UB(t) = Observed - LowerPI if Observed > LowerPI else 0
UB_t = max(Observed - LowerPI, 0);

% ----------------------------
% 5) Compute totals for YEAR_START:YEAR_END
% ----------------------------
TotalExcess = NaN; TotalLB = NaN; TotalUB = NaN;
idxWindow = [];

if ~isempty(timeCol)
    timeVar = T.(timeCol);

    % Convert time variable to years
    if isdatetime(timeVar)
        yr = year(timeVar);
    else
        yr = localToDouble(timeVar);
    end

    idxWindow = (yr >= YEAR_START) & (yr <= YEAR_END);

    TotalExcess = nansum(Excess_t(idxWindow));
    TotalLB     = nansum(LB_t(idxWindow));
    TotalUB     = nansum(UB_t(idxWindow));
end

% ----------------------------
% 6) Create output tables
% ----------------------------
Tout = T;
Tout.Excess     = Excess_t;
Tout.Excess_LB  = LB_t;
Tout.Excess_UB  = UB_t;

summaryT = table(YEAR_START, YEAR_END, TotalExcess, TotalLB, TotalUB);

% ----------------------------
% 7) Save outputs
% ----------------------------
outMain = fullfile(fpath, [baseName '_excess' ext]);

switch lower(ext)
    case '.csv'
        % Full output
        writetable(Tout, outMain);

        % Summary output
        outSummary = fullfile(fpath, sprintf('%s_excess_summary_%d_%d.csv', baseName, YEAR_START, YEAR_END));
        writetable(summaryT, outSummary);

    case {'.xlsx','.xls'}
        % Full output and summary as separate sheets in same workbook
        writetable(Tout, outMain, 'Sheet', 'excess_by_t');
        writetable(summaryT, outMain, 'Sheet', sprintf('totals_%d_%d', YEAR_START, YEAR_END));

    otherwise
        % Fallback to CSV
        warning('Unknown extension "%s". Saving outputs as CSV.', ext);

        outMain = fullfile(fpath, [baseName '_excess.csv']);
        writetable(Tout, outMain);

        outSummary = fullfile(fpath, sprintf('%s_excess_summary_%d_%d.csv', baseName, YEAR_START, YEAR_END));
        writetable(summaryT, outSummary);
end

% Optional: save filtered window dataset (CSV)
if SAVE_FILTERED_WINDOW_FILE && ~isempty(idxWindow)
    Tw = Tout(idxWindow, :);
    outWindow = fullfile(fpath, sprintf('%s_excess_%d_%d.csv', baseName, YEAR_START, YEAR_END));
    writetable(Tw, outWindow);
end

% ----------------------------
% 8) Print confirmation
% ----------------------------
disp('Done.');
disp(['Input file:  ' infile]);
disp(['Output file: ' outMain]);
if ~isempty(timeCol)
    fprintf('Totals (%d–%d):\n', YEAR_START, YEAR_END);
    disp(summaryT);
else
    disp('Totals not computed (no time/year/date column detected).');
end

% ============================
% Local helper functions
% ============================
function colName = localGetCol(rawNames, normNames, candidates)
% Return the FIRST matching column name from "candidates" based on normalized names.
    colName = '';
    for i = 1:numel(candidates)
        c = regexprep(lower(candidates{i}), '[^a-z0-9]', '');
        k = find(strcmp(normNames, c), 1, 'first');
        if ~isempty(k)
            colName = rawNames{k};
            return;
        end
    end
end

function x = localToDouble(v)
% Convert a table variable to double safely (numeric, logical, string, cell).
    if isnumeric(v)
        x = double(v);
    elseif islogical(v)
        x = double(v);
    elseif isdatetime(v)
        x = double(year(v)); % only used if user mistakenly passes datetime where year expected
    elseif iscell(v)
        x = cellfun(@(z) str2double(string(z)), v);
    else
        x = str2double(string(v));
    end
end
