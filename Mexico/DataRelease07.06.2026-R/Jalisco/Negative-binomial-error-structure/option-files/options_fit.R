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
# options_fit.R
# User configuration for model fitting
# Part of PhenoGrowthR: phenomenological epidemic growth-model fitting and forecasting
# =============================================================================

options_fit <- function() {

  # =========================================================================
  # Dataset Properties
  # =========================================================================
  cadfilename1 <- "JALISCO_2025-08-18_2026-06-29-trimmed"   # Data file name (without .txt)
  caddisease   <- "Measles"                # Disease name
  datatype     <- "cases"                  # Data type (cases, deaths, etc.)

  # =========================================================================
  # Parameter Estimation
  # =========================================================================
  method1        <- 3     # 0=LSQ, 1=MLE Poisson, 3=MLE NB(mean+alpha*mean),
                          # 4=MLE NB(mean+alpha*mean^2), 5=MLE NB(mean+alpha*mean^d)
  dist1          <- 3     # Error structure: 0=Normal, 1=Poisson, 2=NB(factor*mean),
                          # 3=NB(mean+alpha*mean), 4=NB(mean+alpha*mean^2), 5=NB(mean+alpha*mean^d)
  numstartpoints <- 10    # Number of random starting points for optimization
  M              <- 300   # Number of bootstrap realizations

  # =========================================================================
  # Growth Model
  # =========================================================================
  # Model flags:  EXP = -1, GGM = 0, GLM = 1, GRM = 2, LM = 3, RICH = 4, GOM = 5
  # Accepts a single model or a vector — the toolbox loops over each one
  # automatically and saves separate, distinctly-named outputs per model, so
  # you don't need to edit this file and re-run for each model.
  #   flag1 <- 1                          # just GLM
  flag1 <- c(1, 3, 4)                  # GLM, LM, RICH
  #   flag1 <- c(-1, 0, 1, 2, 3, 4, 5)      # every standard model
  # Model names are derived automatically — no need to set them separately.
  #flag1       <- 1        # GLM
  fixI0       <- 1        # Fix I0 to first data point value (>0 means fix to this value)

  # =========================================================================
  # Rolling Window Parameters
  # =========================================================================
  windowsize1 <- 40       # Calibration window size (number of data points)
  tstart1     <- 1        # Start index for rolling window
  tend1       <- 1        # End index (last permissible start index) for rolling window
  window_step <- 1        # Spacing between consecutive window starts. E.g. tstart1=1,
                          # tend1=9, window_step=4 gives windows starting at 1, 5, 9
                          # (i.e. weeks 1-25, 5-29, 9-33 for windowsize1=25).
                          # Default 1 = original consecutive-week behavior.

  # =========================================================================
  # Reproducibility
  # =========================================================================
  seed <- 123    # Random seed for multistart + bootstrap; NA = no fixed seed

  # =========================================================================
  # Bootstrap Performance Tuning
  # =========================================================================
  # These apply ONLY to the 300 bootstrap replicates, not the main point
  # estimate (which still uses the full `numstartpoints` above, and the
  # optim()/ODE defaults baked into fit_model.R / growth_models.R). Bootstrap
  # replicates don't need per-fit precision beyond what the resampling noise
  # already introduces, so looser settings here trade a little per-replicate
  # precision for substantially faster, more evenly-loaded bootstrapping.
  numstartpoints_boot <- 2      # Multistart attempts per bootstrap replicate (1 warm start + this many random starts)
  boot_maxit          <- 2000   # optim() iteration cap for bootstrap fits (main fit keeps 20000)
  boot_factr          <- 1e7    # optim() convergence tolerance for bootstrap fits (R's own default; main fit keeps 1e-6)
  boot_rtol           <- 1e-4   # ODE solver relative tolerance for bootstrap fits (main fit keeps 1e-6)
  boot_atol           <- 1e-4   # ODE solver absolute tolerance for bootstrap fits (main fit keeps 1e-6)

  # =========================================================================
  # Advanced: Uncertainty Method & Identifiability (requires FME package)
  # =========================================================================
  # uncertainty_method: "bootstrap" (default, no extra dependencies) or "mcmc"
  # (samples the real parameter posterior via FME::modMCMC instead of
  # resample-and-refit; usually faster AND more principled for hard-to-fit
  # windows, since it reveals genuinely flat/wide posteriors for weakly
  # identifiable parameters rather than a bootstrap histogram that can look
  # artificially peaked). Install with install.packages("FME").
  uncertainty_method <- "mcmc"
  mcmc_niter         <- 5000   # total MCMC iterations per chain (only used if uncertainty_method="mcmc")
  mcmc_burnin        <- NULL   # NULL = niter/5

  # n_chains: number of independent MCMC chains run from dispersed starting
  # points (1 warm start from the point estimate + n_chains-1 random starts).
  # Enables the Gelman-Rubin R-hat convergence diagnostic (needs >=2 chains
  # and the coda package, which comes with FME). Chains dispatch in parallel
  # when use_parallel=TRUE below, reusing the same worker pool as bootstrap —
  # so more chains costs little extra wall-clock time on a multi-core machine.
  # Set n_chains <- 1 to skip this and match the old single-chain behavior.
  n_chains       <- 4
  rhat_threshold <- 1.1   # R-hat below this = converged (standard rule of thumb)

  # run_identifiability_check: if TRUE, screens the fitted model's free
  # parameters for practical identifiability (FME::sensFun + collin()) right
  # after AICc, before spending time on bootstrap/MCMC uncertainty. A
  # collinearity index above identifiability_threshold means this data
  # window cannot jointly constrain that parameter combination — worth
  # knowing before trusting the resulting histograms. See identifiability.R.
  run_identifiability_check <- FALSE
  identifiability_threshold <- 20

  # =========================================================================
  # Parallel Processing
  # =========================================================================
  use_parallel <- TRUE    # TRUE = parallel bootstrap; FALSE = sequential (safe/debug mode)
  n_cores      <- parallel::detectCores() - 1  # Workers to use; override e.g. n_cores <- 10

  return(list(
    cadfilename1   = cadfilename1,
    caddisease     = caddisease,
    datatype       = datatype,
    method1        = method1,
    dist1          = dist1,
    numstartpoints = numstartpoints,
    M              = M,
    flag1          = flag1,
    fixI0          = fixI0,
    windowsize1    = windowsize1,
    tstart1        = tstart1,
    tend1          = tend1,
    window_step    = window_step,
    seed           = seed,
    numstartpoints_boot = numstartpoints_boot,
    boot_maxit     = boot_maxit,
    boot_factr     = boot_factr,
    boot_rtol      = boot_rtol,
    boot_atol      = boot_atol,
    uncertainty_method = uncertainty_method,
    mcmc_niter     = mcmc_niter,
    mcmc_burnin    = mcmc_burnin,
    n_chains       = n_chains,
    rhat_threshold = rhat_threshold,
    run_identifiability_check = run_identifiability_check,
    identifiability_threshold = identifiability_threshold,
    use_parallel   = use_parallel,
    n_cores        = n_cores
  ))
}
