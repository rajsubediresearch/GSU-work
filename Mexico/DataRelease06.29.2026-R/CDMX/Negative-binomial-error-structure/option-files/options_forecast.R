# =============================================================================
# =============================================================================
# method1 — estimation method (selects the objective function / likelihood
# actually optimized in objective_function(), growth_models.R)
# =============================================================================
#   0 = LSQ            Least squares (SSE). Assumes ~constant-variance/normal
#                       errors implicitly. No dispersion parameter (alpha, d
#                       are fixed at 0 in get_bounds()).
#   1 = MLE Poisson     Poisson negative log-likelihood: Var = Mean. No
#                       dispersion parameter (alpha, d fixed at 0).
#   3 = MLE NB1         Negative binomial, Var = Mean + alpha*Mean
#                       (dispersion scales linearly with the mean).
#                       alpha is a free parameter (d fixed at 1).
#   4 = MLE NB2         Negative binomial, Var = Mean + alpha*Mean^2
#                       (dispersion scales quadratically with the mean).
#                       alpha is a free parameter (d fixed at 1).
#   5 = MLE NB(d)       Negative binomial, Var = Mean + alpha*Mean^d
#                       (dispersion exponent d also estimated -- most
#                       flexible/general variance form). alpha AND d are
#                       both free parameters.
# =============================================================================
# dist1 — error structure used for SIMULATING bootstrap replicate data in
# AddErrorStructure() (growth_models.R). NOT an independent choice for MLE
# methods: for method1 in {1,3,4,5}, Run_Fit_GrowthModels.R /
# Run_Forecasting_GrowthModels.R / screen_all_models() (identifiability.R)
# all silently force dist1 to match method1 right before use, so the
# likelihood actually optimized and the noise actually simulated can't get
# out of sync. Setting dist1 yourself for these methods is harmless but
# redundant -- it's only a genuinely free, independent choice when
# method1 == 0 (LSQ), where the fit doesn't depend on a likelihood at all
# and dist1 purely controls how bootstrap noise is simulated on top of it.
# =============================================================================
#   0 = Normal          Var = factor1^2 (constant variance, from residual SD
#                       of the LSQ fit). Only meaningful with method1 == 0.
#   1 = Poisson         Var = Mean. Pairs with method1 == 1, or usable
#                       standalone with method1 == 0.
#   2 = NB(factor*mean) Var = factor1 * Mean, factor1 estimated empirically
#                       from a rolling mean-variance ratio of the data
#                       (getMeanVarianceRatio()). Only meaningful with
#                       method1 == 0 (falls back to dist1=1/Poisson if the
#                       empirical ratio comes out < 1, i.e. underdispersed).
#   3 = NB1             Var = Mean + alpha*Mean. Pairs with method1 == 3.
#   4 = NB2             Var = Mean + alpha*Mean^2. Pairs with method1 == 4.
#   5 = NB(d)           Var = Mean + alpha*Mean^d. Pairs with method1 == 5.
# options_forecast.R
# User configuration for forecasting
# Part of PhenoGrowthR: phenomenological epidemic growth-model fitting and forecasting
# =============================================================================

options_forecast <- function() {

  # =========================================================================
  # Dataset Properties
  # =========================================================================
  cadfilename1 <- "CIUDAD_DE_MEXICO_2025-08-04_2026-06-22-trimmed"   # Data file name (without .txt)
  caddisease   <- "Measles"                # Disease name
  datatype     <- "cases"                  # Data type

  # =========================================================================
  # Parameter Estimation
  # =========================================================================
  method1        <- 3     # see top for details
  dist1          <- 3     # see top for details
  numstartpoints <- 10    # Random starting points
  M              <- 300   # Bootstrap realizations

  # =========================================================================
  # Growth Model
  # =========================================================================
  # EXP = -1, GGM = 0, GLM = 1, GRM = 2, LM = 3, RICH = 4, GOM = 5
  # Accepts a single model or a vector — the toolbox loops over each one
  # automatically and saves separate, distinctly-named outputs per model.
  #   flag1 <- 1                          # just GLM
  flag1 <- c(1, 3, 4)                  # GLM, LM, RICH
  #   flag1 <- c(-1, 0, 1, 2, 3, 4, 5)      # every standard model
  # Model names are derived automatically — no need to set them separately.
  #flag1       <- 1        # GLM
  fixI0       <- 1        # Fix I0 to first data point

  # =========================================================================
  # Forecasting Parameters
  # =========================================================================
  getperformance   <- 1   # 0=real-time (no ground truth), 1=retrospective
  forecastingperiod <- 2  # n-week ahead forecast horizon

  # =========================================================================
  # Rolling Window Parameters
  # =========================================================================
  windowsize1 <- 29       # Calibration window size
  tstart1     <- 1        # Start index
  tend1       <- 1        # End index (last permissible start index)
  window_step <- 1        # Spacing between consecutive window starts (default 1 = consecutive weeks).
                          # E.g. tstart1=1, tend1=9, window_step=4 -> starts at 1, 5, 9.

  # =========================================================================
  # Reproducibility
  # =========================================================================
  seed <- 123    # Random seed for multistart + bootstrap; NA = no fixed seed

  # =========================================================================
  # Bootstrap Performance Tuning
  # =========================================================================
  numstartpoints_boot <- 2      # Multistart attempts per bootstrap replicate (1 warm start + this many random starts)
  boot_maxit          <- 2000   # optim() iteration cap for bootstrap fits (main fit keeps 20000)
  boot_factr          <- 1e7    # optim() convergence tolerance for bootstrap fits (R's own default; main fit keeps 1e-6)
  boot_rtol           <- 1e-4   # ODE solver relative tolerance for bootstrap fits (main fit keeps 1e-6)
  boot_atol           <- 1e-4   # ODE solver absolute tolerance for bootstrap fits (main fit keeps 1e-6)

  # =========================================================================
  # Advanced: Uncertainty Method & Identifiability (requires FME package)
  # =========================================================================
  # See options_fit.R for full explanation. "bootstrap" (default) or "mcmc".
  uncertainty_method <- "mcmc"
  mcmc_niter         <- 5000
  mcmc_burnin        <- NULL
  n_chains           <- 4       # independent MCMC chains for R-hat (1 = old single-chain behavior)
  rhat_threshold     <- 1.1

  run_identifiability_check <- FALSE
  identifiability_threshold <- 20

  # =========================================================================
  # Parallel Processing
  # =========================================================================
  use_parallel <- TRUE    # TRUE = parallel bootstrap; FALSE = sequential (safe/debug mode)
  n_cores      <- parallel::detectCores() - 1  # Workers to use; override e.g. n_cores <- 10

  return(list(
    cadfilename1    = cadfilename1,
    caddisease      = caddisease,
    datatype        = datatype,
    method1         = method1,
    dist1           = dist1,
    numstartpoints  = numstartpoints,
    M               = M,
    flag1           = flag1,
    fixI0           = fixI0,
    getperformance  = getperformance,
    forecastingperiod = forecastingperiod,
    windowsize1     = windowsize1,
    tstart1         = tstart1,
    tend1           = tend1,
    window_step     = window_step,
    seed            = seed,
    numstartpoints_boot = numstartpoints_boot,
    boot_maxit      = boot_maxit,
    boot_factr      = boot_factr,
    boot_rtol       = boot_rtol,
    boot_atol       = boot_atol,
    uncertainty_method = uncertainty_method,
    mcmc_niter      = mcmc_niter,
    mcmc_burnin     = mcmc_burnin,
    n_chains        = n_chains,
    rhat_threshold  = rhat_threshold,
    run_identifiability_check = run_identifiability_check,
    identifiability_threshold = identifiability_threshold,
    use_parallel    = use_parallel,
    n_cores         = n_cores
  ))
}
