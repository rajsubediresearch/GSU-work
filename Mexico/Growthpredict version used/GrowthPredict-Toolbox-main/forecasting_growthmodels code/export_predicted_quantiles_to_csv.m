% =============================================================================
%
% Convert a GrowthPredict QuantilesCalibration .mat output into a CSV with
% proper column headers, saved in the same \output folder.
% Raj Subedi | Chowell Lab | GSU | 2026-07-12
%
% What's actually in this file (from computeQuantiles.m in the
% GrowthPredict-Toolbox source, gchowell/GrowthPredict-Toolbox):
%
%   quantilescs is a [n_weeks x 23] matrix.
%   - Row i        = calibration week i (1 = first week of the calibration
%                    window, up to your calibrationperiod).
%   - Columns 1-23 = quantiles of the bootstrap-fitted curve distribution
%                    at that week, at these exact alpha levels:
%                    0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250,
%                    0.300, 0.350, 0.400, 0.450, 0.500 (median), 0.550,
%                    0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900,
%                    0.950, 0.975, 0.990
%                  So column 12 = median, column 2 = 2.5% bound,
%                  column 22 = 97.5% bound.
%
% Usage:
%   1. Set RUNNAME below to the part of the filename that comes AFTER
%      "QuantilesCalibration-", e.g. for:
%        QuantilesCalibration-growthModel-JALISCO_2025-08-18_2026-06-29-trimmed.txt-flag1-2-fixI0-1-method-3-dist-3-tstart-1-tend-1-calibrationperiod-40-forecastingperiod-3.mat
%      RUNNAME is everything after "QuantilesCalibration-" and before ".mat".
%   2. Run. Two CSVs are written into \output, both with a "Week" column
%      plus one column per quantile level (labeled by its alpha value,
%      e.g. "q0.500" for the median):
%        - QuantilesCalibration-<runname>.csv           (full name, as-is)
%        - <first 16 characters of the above name>.csv  (short name)
% =============================================================================

clc; clear;

% ── 0. USER INPUT — only this line should need to change ───────────────────

runname = 'Forecast-growthModel-JALISCO_2025-08-18_2026-07-06.txt-flag1-1-fixI0-1-method-3-dist-3-tstart-1-tend-1-calibrationperiod-40-forecastingperiod-4';

% ── 1. BUILD PATHS ───────────────────────────────────────────────────────

output_dir = fullfile(pwd, 'output');
matfile    = fullfile(output_dir, ['QuantilesCalibration-' runname '.mat']);

fullbase   = ['QuantilesCalibration-' runname];
shortbase  = fullbase(1:min(16, length(fullbase)));

csvfile_full  = fullfile(output_dir, [fullbase '.csv']);
csvfile_short = fullfile(output_dir, [shortbase '.csv']);

if ~isfile(matfile)
    error('Could not find .mat file:\n  %s', matfile);
end

fprintf('Loading: %s\n', matfile);
S = load(matfile);

if ~isfield(S, 'quantilescs')
    error(['Expected variable "quantilescs" not found in this .mat file.\n' ...
           'Variables actually present: %s'], strjoin(fieldnames(S), ', '));
end

quantilescs = S.quantilescs;
[nweeks, ncols] = size(quantilescs);

% The exact alpha levels used by computeQuantiles.m in GrowthPredict-Toolbox
alphaquantiles = [0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, ...
                   0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, ...
                   0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990];

% ── 2. BUILD HEADERS ─────────────────────────────────────────────────────

if ncols == numel(alphaquantiles)
    colnames = arrayfun(@(a) sprintf('q%.3f', a), alphaquantiles, 'UniformOutput', false);
else
    % Toolbox version produced a different number of quantile columns than
    % expected — fall back to generic labels rather than mislabeling them.
    fprintf('WARNING: expected %d quantile columns, found %d. Using generic column labels.\n', ...
        numel(alphaquantiles), ncols);
    colnames = arrayfun(@(c) sprintf('Col%d', c), 1:ncols, 'UniformOutput', false);
end

% ── 3. WRITE CSV ─────────────────────────────────────────────────────────

week = (1:nweeks)';
T = array2table([week, quantilescs]);
T.Properties.VariableNames = [{'Week'}, colnames];

writetable(T, csvfile_full);
fprintf('Saved (full name):  %s  [%d weeks x %d quantile levels]\n', csvfile_full, nweeks, ncols);

writetable(T, csvfile_short);
fprintf('Saved (short name): %s  [%d weeks x %d quantile levels]\n', csvfile_short, nweeks, ncols);

fprintf('Done.\n');