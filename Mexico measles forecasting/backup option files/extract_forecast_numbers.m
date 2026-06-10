% =============================================================================
% extract_forecast_numbers.m
% Extract forecast quantiles from GrowthPredict .mat output file
% Raj Subedi | Chowell Lab | GSU | 2026-06-08
%
% Usage: Run this script in MATLAB from the GrowthPredict-Toolbox directory
% =============================================================================

clc; clear;

% ── 1. LOAD THE .MAT FILE ────────────────────────────────────────────────────

matfile = 'D:\PhD coursework\GRA\Spring 2026\Mexico_Measles\DGE_data\outputs\08_matlab_inputs\GrowthPredict-Toolbox-main\forecasting_growthmodels code\output\Forecast-growthModel-cdmx_weekly-trimmed.txt-flag1-2-fixI0-1-method-3-dist-3-tstart-1-tend-1-calibrationperiod-26-forecastingperiod-6.mat';

fprintf('Loading: %s\n\n', matfile);
load(matfile);

% ── 2. INSPECT AVAILABLE VARIABLES ───────────────────────────────────────────

fprintf('Variables in .mat file:\n');
whos

% ── 3. EXTRACT FORECAST TRAJECTORIES ─────────────────────────────────────────
% GrowthPredict typically stores forecast curves in one of these variable names
% Try each until one works

if exist('forecast_model1', 'var')
    fc = forecast_model1;
    fprintf('\nFound: forecast_model1  [%d x %d]\n', size(fc,1), size(fc,2));
else
    fprintf('\nWARNING: Could not find forecast matrix.\n');
    fc = [];
end

% ── 4. EXTRACT OBSERVED DATA ──────────────────────────────────────────────────

if exist('data1', 'var')
    obs = data1;
elseif exist('data', 'var')
    obs = data;
elseif exist('Ydata', 'var')
    obs = Ydata;
else
    obs = [];
    fprintf('WARNING: Could not find observed data variable.\n');
end

% ── 5. COMPUTE QUANTILES PER FORECAST WEEK ───────────────────────────────────

if exist('data1', 'var')
    obs = data1(:,2);  % column 2 = case counts
elseif exist('data', 'var')
    obs = data(:,2);
else
    obs = [];
    fprintf('WARNING: Could not find observed data variable.\n');
end

% ── 5. USE PRE-COMPUTED QUANTILES IF AVAILABLE ───────────────────────────────
% quantilesf is 6x23 — rows = forecast weeks, cols = quantile levels
% Standard quantile columns in GrowthPredict:
% col 1=0.025, col 2=0.05, ..., col 12=0.5 (median), ..., col 23=0.975

if exist('quantilesf', 'var') && ~isempty(quantilesf)
    fprintf('\n================================================================\n');
    fprintf('  JALISCO — 6-WEEK AHEAD FORECAST (GRM, NegBin, fixI0=1)\n');
    fprintf('  Using pre-computed quantilesf [%dx%d]\n', size(quantilesf,1), size(quantilesf,2));
    fprintf('  Calibration: weeks 5-40 (n=37) | Forecast: weeks 41-46\n');
    fprintf('================================================================\n');
    fprintf('%-12s %-10s %-12s %-12s %-12s\n', ...
            'Fore.Week', 'Mean', 'Median', 'CI_low(2.5%)', 'CI_high(97.5%)');
    fprintf('%-12s %-10s %-12s %-12s %-12s\n', ...
            '---------', '----', '------', '------------', '--------------');

    forecastperiod = size(quantilesf, 1);  % should be 6
    n_quantiles    = size(quantilesf, 2);  % should be 23

    % Identify quantile columns
    % GrowthPredict uses 23 quantiles: 0.025,0.05,0.10,...,0.90,0.95,0.975
    % col 1 = 2.5%, col 12 = 50% (median), col 23 = 97.5%
    col_lo  = 1;   % 2.5%
    col_med = 12;  % 50%
    col_hi  = 23;  % 97.5%

    % Mean from bootstrap samples
    fc_forecast_only = forecast_model1(end-forecastperiod+1:end, :);

    forecast_totals = zeros(forecastperiod, 5);
    for w = 1:forecastperiod
        mn  = mean(fc_forecast_only(w, :));
        med = quantilesf(w, col_med);
        lo  = quantilesf(w, col_lo);
        hi  = quantilesf(w, col_hi);
        fprintf('Week +%-4d   %-10.1f %-12.1f %-12.1f %-12.1f\n', w, mn, med, lo, hi);
        forecast_totals(w,:) = [w, mn, med, lo, hi];
    end

    % WC period total
    total_mean = sum(forecast_totals(:,2));
    total_lo   = sum(forecast_totals(:,4));
    total_hi   = sum(forecast_totals(:,5));

    fprintf('\n----------------------------------------------------------------\n');
    fprintf('  TOTAL FORECAST — 6 weeks (covers World Cup period):\n');
    fprintf('  Mean:   %.0f cases\n', total_mean);
    fprintf('  95%% PI: %.0f – %.0f cases\n', total_lo, total_hi);
    fprintf('----------------------------------------------------------------\n');

    % Also print parameter estimates
    fprintf('\n  PARAMETER ESTIMATES (mean, 95%% CI):\n');
    fprintf('  r:  %.3f (%.3f – %.3f)\n', param_r(1), param_r(2), param_r(3));
    fprintf('  p:  %.3f (%.3f – %.3f)\n', param_p(1), param_p(2), param_p(3));
    fprintf('  a:  %.3f (%.3f – %.3f)\n', param_a(1), param_a(2), param_a(3));
    fprintf('  K:  %.1f (%.1f – %.1f)\n', param_K(1), param_K(2), param_K(3));
    fprintf('  AICc: %.2f\n', AICc);

    % Last 5 observed weeks
    if ~isempty(obs)
        fprintf('\n  Last 5 observed weeks:\n');
        n_obs = length(obs);
        for i = max(1,n_obs-4):n_obs
            fprintf('  Week %2d (t=%2d): %d cases\n', i+4, i, obs(i));
        end
    end

    % Save summary
    [mat_dir, mat_name, ~] = fileparts(matfile);
    outfile = fullfile(mat_dir, [mat_name '_summary.txt']);
    fid = fopen(outfile, 'w');
    fprintf(fid, '================================================================\n');
    fprintf(fid, '  JALISCO FORECAST SUMMARY\n');
    fprintf(fid, '  GRM | NegBin | fixI0=1 | Calibration=37wk | Forecast=6wk\n');
    fprintf(fid, '  Generated: %s\n', datestr(now));
    fprintf(fid, '================================================================\n\n');
    fprintf(fid, 'PARAMETERS:\n');
    fprintf(fid, '  r:    %.3f (95%% CI: %.3f-%.3f)\n', param_r(1), param_r(2), param_r(3));
    fprintf(fid, '  p:    %.3f (95%% CI: %.3f-%.3f)\n', param_p(1), param_p(2), param_p(3));
    fprintf(fid, '  a:    %.3f (95%% CI: %.3f-%.3f)\n', param_a(1), param_a(2), param_a(3));
    fprintf(fid, '  K:    %.1f (95%% CI: %.1f-%.1f)\n', param_K(1), param_K(2), param_K(3));
    fprintf(fid, '  AICc: %.2f\n\n', AICc);
    fprintf(fid, 'WEEKLY FORECASTS:\n');
    fprintf(fid, '%-12s %-10s %-12s %-12s %-12s\n', 'Week', 'Mean', 'Median', 'CI_2.5', 'CI_97.5');
    for w = 1:forecastperiod
        fprintf(fid, 'Week +%-4d   %-10.1f %-12.1f %-12.1f %-12.1f\n', ...
            w, forecast_totals(w,2), forecast_totals(w,3), ...
            forecast_totals(w,4), forecast_totals(w,5));
    end
    fprintf(fid, '\nTOTAL 6-week (WC period):\n');
    fprintf(fid, '  Mean:   %.0f\n  95%% PI: %.0f - %.0f\n', total_mean, total_lo, total_hi);
    fclose(fid);
    fprintf('\n✔ Summary saved to:\n    %s\n', outfile);

else
    fprintf('\nquantilesf not found — extracting directly from forecast_model1...\n');
end

fprintf('\n✅ Done.\n');
