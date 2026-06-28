## Running mif2 and pmcmc for our model

## We observe daily incidence counts. A random starting point failed to yield estimates that gives a finite value of the log-likelihood, 
## which is needed to run pmcmc. 
## Even after giving good starting values, sometimes we saw some bad estimates, so we for each iteration, we made n_rep = 10 trials and
## chose the estimates with least absolute bias. In case all trial fails, we also give the true value as the starting point as an option.  

# install.packages("pomp")
library("pomp")

run_mif2_repeated <- function(
    obs_data,              # A data frame with 2 columns: "time" = days and "cases" = incidence counts
    N_init,                # Initial number of susceptibles
    start_params,          # Starting value of the MLE. A named vector: c(beta=, gamma=, rho=)
    true_values  = NULL,   # True values of simulated data. A named vector: c(beta=, gamma=, rho=)  
    Nmif         = 200,    # number of filtering iterations to perform (check mif2 function)
    Np_mif       = 2000,   # The number of particles to use. (check Np of mif2 function)
    Np_eval      = 2000,   # Number of particles for evaluating the likelihood using pfilter of pomp
    n_eval       = 10,     # Number of times repeating pfilter
    n_rep        = 10     # Number of times, we are repeating the MLE evaluations
) {
  
  # ── Csnippets ────────────────────────────────────────────────────────────────
  sir_step <- Csnippet("
    double rate[2];
    double dN[2];
    rate[0] = beta * I / N;
    rate[1] = gamma;
    reulermultinom(1, nearbyint(S), &rate[0], dt, &dN[0]);
    reulermultinom(1, nearbyint(I), &rate[1], dt, &dN[1]);
    S += -dN[0];
    I +=  dN[0] - dN[1];
    R +=  dN[1];
    C +=  dN[0];
  ")
  
  sir_init <- Csnippet("
    double M = nearbyint(rho * N);
    S = nearbyint(N);
    I = M;
    R = 0;
    C = 0;
  ")
  
  sir_dmeas <- Csnippet("
    lik = ((int) nearbyint(C) == cases) ? (give_log ? 0.0 : 1.0)
                                        : (give_log ? R_NegInf : 0.0);
  ")
  
  sir_rmeas <- Csnippet("
    cases = (int) nearbyint(C);
  ")
  
  sir_dprior <- Csnippet("
    lik = dgamma(beta,  0.1, 10.0, 1)
        + dgamma(gamma, 0.1, 10.0, 1)
        + dunif(rho, 0.0, 1.0, 1);
    if (!give_log) lik = exp(lik);
  ")
  
  # ── Validate start_params ────────────────────────────────────────────────────
  required <- c("beta", "gamma", "rho")
  if (!all(required %in% names(start_params)))
    stop("start_params must contain: beta, gamma, rho")
  
  all_params <- c(start_params[required], N = N_init)
  
  # ── Build base pomp object ───────────────────────────────────────────────────
  sir_model <- pomp(
    data       = obs_data,
    times      = "time",
    t0         = 0,
    rprocess   = euler(sir_step, delta.t = 1/24),
    rinit      = sir_init,
    dmeasure   = sir_dmeas,
    rmeasure   = sir_rmeas,
    dprior     = sir_dprior,
    accumvars  = "C",
    statenames = c("S", "I", "R", "C"),
    paramnames = c("beta", "gamma", "rho", "N"),
    params     = all_params,
    partrans   = parameter_trans(
      log   = c("beta", "gamma"),
      logit = "rho"
    )
  )
  
  rws <- rw_sd(beta = 0.02, gamma = 0.02, rho = 0.02)
  
  # ── Storage ──────────────────────────────────────────────────────────────────
  all_estimates <- vector("list", n_rep)
  all_ll        <- rep(NA_real_, n_rep)
  all_ll_se     <- rep(NA_real_, n_rep)
  all_mif_out   <- vector("list", n_rep)
  all_times     <- rep(NA_real_, n_rep)      # elapsed seconds per replicate
  
  cat(sprintf("\n=== Running %d MIF2 replicates ===\n", n_rep))
  
  # ── n_rep replicates ─────────────────────────────────────────────────────────
  for (i in seq_len(n_rep)) {
    cat(sprintf("\n--- Replicate %d/%d ---\n", i, n_rep))
    tryCatch({
      t_start <- proc.time()["elapsed"]
      
      mif_out <- mif2(
        sir_model,
        Nmif                = Nmif,
        Np                  = Np_mif,
        rw.sd               = rws,
        cooling.fraction.50 = 0.5
      )
      pf_list  <- replicate(n_eval, pfilter(mif_out, Np = Np_eval), simplify = FALSE)
      ll_vals  <- sapply(pf_list, logLik)
      ll_est   <- logmeanexp(ll_vals, se = TRUE)
      
      t_end <- proc.time()["elapsed"]
      
      all_mif_out[[i]]   <- mif_out
      all_estimates[[i]] <- coef(mif_out)[required]
      all_ll[[i]]        <- ll_est[1]
      all_ll_se[[i]]     <- ll_est[2]
      all_times[[i]]     <- as.numeric(t_end - t_start)
      
      cat(sprintf("  LL = %.4f (SE = %.4f)\n", ll_est[1], ll_est[2]))
      cat(sprintf("  beta = %.4f, gamma = %.4f, rho = %.4f\n",
                  all_estimates[[i]]["beta"],
                  all_estimates[[i]]["gamma"],
                  all_estimates[[i]]["rho"]))
      cat(sprintf("  Time = %.2f sec\n", all_times[[i]]))
      
    }, error = function(e) {
      cat(sprintf("  Replicate %d failed: %s\n", i, e$message))
    })
  }
  
  # ── Check non-NA ─────────────────────────────────────────────────────────────
  non_na_idx <- which(!is.na(all_ll))
  n_non_na   <- length(non_na_idx)
  cat(sprintf("\n=== %d/%d replicates returned non-NA likelihood ===\n",
              n_non_na, n_rep))
  
  # ── Extra run with start_params if all NA ────────────────────────────────────
  extra_est    <- NULL
  extra_ll     <- NA_real_
  extra_ll_se  <- NA_real_
  extra_run    <- NULL
  extra_time   <- NA_real_               # time for extra run
  
  if (n_non_na == 0) {
    cat("\nAll replicates gave NA. Running extra replicate with start_params...\n")
    tryCatch({
      t_start <- proc.time()["elapsed"]
      
      sir_model_start <- pomp(sir_model, params = all_params)
      mif_extra <- mif2(
        sir_model_start,
        Nmif                = Nmif,
        Np                  = Np_mif,
        rw.sd               = rws,
        cooling.fraction.50 = 0.5
      )
      pf_list     <- replicate(n_eval, pfilter(mif_extra, Np = Np_eval), simplify = FALSE)
      ll_vals     <- sapply(pf_list, logLik)
      ll_est      <- logmeanexp(ll_vals, se = TRUE)
      
      t_end <- proc.time()["elapsed"]
      
      extra_run   <- mif_extra
      extra_ll    <- ll_est[1]
      extra_ll_se <- ll_est[2]
      extra_est   <- coef(mif_extra)[required]
      extra_time  <- as.numeric(t_end - t_start)
      
      cat(sprintf("  Extra run LL = %.4f (SE = %.4f)\n", extra_ll, extra_ll_se))
      cat(sprintf("  beta = %.4f, gamma = %.4f, rho = %.4f\n",
                  extra_est["beta"], extra_est["gamma"], extra_est["rho"]))
      cat(sprintf("  Time = %.2f sec\n", extra_time))
      
    }, error = function(e) {
      cat(sprintf("  Extra run failed: %s\n", e$message))
    })
  }
  
  # ── Select best ───────────────────────────────────────────────────────────────
  best_idx  <- NULL
  best_est  <- NULL
  best_ll   <- NA_real_
  best_time <- NA_real_
  
  if (n_non_na > 0) {
    if (!is.null(true_values)) {
      dists    <- sapply(non_na_idx, function(i) {
        est <- all_estimates[[i]][names(true_values)]
        sum(abs(true_values - est), na.rm = TRUE)
      })
      best_idx <- non_na_idx[which.min(dists)]
      cat(sprintf("\nBest replicate (min param distance): Replicate %d\n", best_idx))
    } else {
      best_idx <- non_na_idx[which.max(all_ll[non_na_idx])]
      cat(sprintf("\nBest replicate (max likelihood): Replicate %d\n", best_idx))
    }
    best_est  <- all_estimates[[best_idx]]
    best_ll   <- all_ll[[best_idx]]
    best_time <- all_times[[best_idx]]
    
  } else if (!is.null(extra_run)) {
    cat("\nBest: extra run with start_params\n")
    best_est  <- extra_est
    best_ll   <- extra_ll
    best_time <- extra_time
  }
  
  # ── Summary table ─────────────────────────────────────────────────────────────
  cat("\n=== SUMMARY TABLE ===\n")
  cat(sprintf("%-5s  %-12s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
              "Rep", "LogLik", "SE", "beta", "gamma", "rho", "Time(sec)"))
  cat(strrep("-", 75), "\n")
  
  for (i in seq_len(n_rep)) {
    est <- all_estimates[[i]]
    cat(sprintf("%-5d  %-12s  %-10s  %-10s  %-10s  %-10s  %-10s%s\n",
                i,
                ifelse(is.na(all_ll[i]),    "NA", sprintf("%.4f", all_ll[i])),
                ifelse(is.na(all_ll_se[i]), "NA", sprintf("%.4f", all_ll_se[i])),
                ifelse(is.null(est),        "NA", sprintf("%.4f", est["beta"])),
                ifelse(is.null(est),        "NA", sprintf("%.4f", est["gamma"])),
                ifelse(is.null(est),        "NA", sprintf("%.4f", est["rho"])),
                ifelse(is.na(all_times[i]), "NA", sprintf("%.2f",  all_times[i])),
                ifelse(!is.null(best_idx) && i == best_idx, "  <-- BEST", "")))
  }
  
  if (!is.null(extra_run)) {
    cat(sprintf("%-5s  %-12s  %-10s  %-10s  %-10s  %-10s  %-10s  <-- BEST (extra)\n",
                "EXT",
                ifelse(is.na(extra_ll),    "NA", sprintf("%.4f", extra_ll)),
                ifelse(is.na(extra_ll_se), "NA", sprintf("%.4f", extra_ll_se)),
                sprintf("%.4f", extra_est["beta"]),
                sprintf("%.4f", extra_est["gamma"]),
                sprintf("%.4f", extra_est["rho"]),
                ifelse(is.na(extra_time),  "NA", sprintf("%.2f",  extra_time))))
  }
  
  cat("\n=== BEST MLE ===\n")
  if (!is.null(best_est)) {
    cat(sprintf("  beta      = %.6f\n", best_est["beta"]))
    cat(sprintf("  gamma     = %.6f\n", best_est["gamma"]))
    cat(sprintf("  rho       = %.6f\n", best_est["rho"]))
    cat(sprintf("  LogLik    = %.4f\n", best_ll))
    cat(sprintf("  Time(sec) = %.2f\n", best_time))
  } else {
    cat("  No valid estimate found.\n")
  }
  
  cat("\n=== TIMING SUMMARY ===\n")
  valid_times <- c(all_times, if (!is.na(extra_time)) extra_time)
  valid_times <- valid_times[!is.na(valid_times)]
  cat(sprintf("  Total time  : %.2f sec (%.2f min)\n",
              sum(valid_times), sum(valid_times) / 60))
  cat(sprintf("  Mean per rep: %.2f sec\n", mean(valid_times)))
  cat(sprintf("  Min / Max   : %.2f / %.2f sec\n",
              min(valid_times), max(valid_times)))
  
  # ── Return ────────────────────────────────────────────────────────────────────
  invisible(list(
    data = obs_data,
    N = N_init,
    all_estimates = all_estimates,
    all_ll        = all_ll,
    all_ll_se     = all_ll_se,
    all_times     = all_times,       # vector of n_rep elapsed times (sec)
    extra_est     = extra_est,
    extra_ll      = extra_ll,
    extra_time    = extra_time,      # extra run time (or NA)
    best_est      = best_est,
    best_ll       = best_ll,
    best_time     = best_time,       # time of best replicate
    best_idx      = best_idx,
    all_mif_out   = all_mif_out
  ))
}


##########################################################################
## Running pmcmc




run_pmcmc_sir <- function(
    obs_data,        # Data frame with two columns: days and incidence count
    N_init,          # population size (integer);
    mle_obj,         # MLE estimate, with best_obj as the MLE; see run_mif2_repeated
    true_values,     # named vector: c(beta=, gamma=, rho=)  [mandatory]
    n_iter   = 1e5,  # number of MCMC iterations per chain
    Np       = 2000, # number of particles for the particle filter
    rw_sd    = 0.007,  # random-walk SD (scalar or per-chain vector)
    n_chains = 2,    # number of independent chains
    burn_in  = 0.5,  # fraction of iterations discarded as burn-in
    seed     = 1234  # RNG seed for reproducibility
) {
  
  mle_est         <- mle_obj$best_est   # named vector: beta, gamma, rho
  required_params <- c("beta", "gamma", "rho")
  
  if (!all(required_params %in% names(mle_est)))
    stop("best_est must contain: beta, gamma, rho")
  
  cat(sprintf("\n=== PMCMC for SIR Model (N = %d) ===\n", N_init))
  cat("MLE starting values:\n")
  cat(sprintf("  beta  = %.6f\n", mle_est["beta"]))
  cat(sprintf("  gamma = %.6f\n", mle_est["gamma"]))
  cat(sprintf("  rho   = %.6f\n", mle_est["rho"]))
  cat("True values:\n")
  cat(sprintf("  beta  = %.6f\n", true_values["beta"]))
  cat(sprintf("  gamma = %.6f\n", true_values["gamma"]))
  cat(sprintf("  rho   = %.6f\n", true_values["rho"]))
  cat(sprintf(
    "\nSettings: n_iter = %d, Np = %d, n_chains = %d, burn_in = %.0f%%\n",
    n_iter, Np, n_chains, burn_in * 100
  ))
  
  # ── Csnippets ───────────────────────────────────────────────────────────────
  sir_step <- Csnippet("
    double rate[2];
    double dN[2];
    rate[0] = beta * I / N;
    rate[1] = gamma;
    reulermultinom(1, nearbyint(S), &rate[0], dt, &dN[0]);
    reulermultinom(1, nearbyint(I), &rate[1], dt, &dN[1]);
    S += -dN[0];
    I +=  dN[0] - dN[1];
    R +=  dN[1];
    C +=  dN[0];
  ")
  
  sir_init <- Csnippet("
    double M = nearbyint(rho * N);
    S = nearbyint(N);
    I = M;
    R = 0;
    C = 0;
  ")
  
  sir_dmeas <- Csnippet("
    lik = ((int) nearbyint(C) == cases) ? (give_log ? 0.0 : 1.0)
                                        : (give_log ? R_NegInf : 0.0);
  ")
  
  sir_rmeas <- Csnippet("
    cases = (int) nearbyint(C);
  ")
  
  sir_dprior <- Csnippet("
    lik = dgamma(beta,  0.1, 10.0, 1)
        + dgamma(gamma, 0.1, 10.0, 1)
        + dunif(rho, 0.0, 1.0, 1);
    if (!give_log) lik = exp(lik);
  ")
  
  # ── Base pomp object ──────────────────────────────────────────────────────
  sir_model <- pomp(
    data       = obs_data,
    times      = "time",
    t0         = 0,
    rprocess   = euler(sir_step, delta.t = 1/24),
    rinit      = sir_init,
    dmeasure   = sir_dmeas,
    rmeasure   = sir_rmeas,
    dprior     = sir_dprior,
    accumvars  = "C",
    statenames = c("S", "I", "R", "C"),
    paramnames = c("beta", "gamma", "rho", "N"),
    params     = c(mle_est[required_params], N = N_init),
    partrans   = parameter_trans(
      log   = c("beta", "gamma"),
      logit = "rho"
    )
  )
  
  # ── Storage ─────────────────────────────────────────────────────────────────
  chain_draws       <- vector("list", n_chains)
  chain_draws_full  <- vector("list", n_chains)
  chain_times       <- rep(NA_real_,  n_chains)   # total wall-clock per chain (sum of all attempts)
  chain_accept      <- rep(NA_real_,  n_chains)
  chain_ll          <- vector("list", n_chains)
  chain_ll_full     <- vector("list", n_chains)
  chain_seeds       <- seed + seq_len(n_chains)
  chain_fail_reason <- vector("list", n_chains)   # per-chain: character vector of attempt failure messages
  chain_start_label <- rep(NA_character_, n_chains)  # which candidate succeeded ("jittered"/"MLE"/"true values")
  chain_start_params <- vector("list", n_chains)     # the actual param values used for the successful run
  chain_attempt_log  <- vector("list", n_chains)     # full per-attempt timing/label/status, for diagnostics
  
  # ── Helper: attempt pmcmc() directly with a candidate starting point ───────
  # No pilot pfilter step anymore — we try the real pmcmc() call and let it
  # succeed or fail on its own. Returns timing + outcome either way.
  .try_pmcmc <- function(model, params, label, rw_sd_ch) {
    t0 <- proc.time()["elapsed"]
    result <- tryCatch({
      res <- pomp::pmcmc(
        pomp(model, params = c(params[required_params], N = N_init)),
        Nmcmc    = n_iter,
        Np       = Np,
        proposal = mvn_diag_rw(
          setNames(rep(rw_sd_ch, length(required_params)), required_params)
        )
      )
      list(ok = TRUE, pmcmc_out = res, msg = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, pmcmc_out = NULL, msg = e$message)
    })
    t1 <- proc.time()["elapsed"]
    result$label   <- label
    result$elapsed <- as.numeric(t1 - t0)
    result
  }
  
  cat(sprintf("\n=== Running %d PMCMC chain(s) ===\n", n_chains))
  
  # ── Run chains ──────────────────────────────────────────────────────────────
  for (ch in seq_len(n_chains)) {
    
    cat(sprintf("\n--- Chain %d/%d ---\n", ch, n_chains))
    set.seed(chain_seeds[ch])
    rw_sd_ch <- if (length(rw_sd) >= ch) rw_sd[ch] else rw_sd[1]
    
    # 1. Jittered start (small perturbation around MLE)
    jittered          <- mle_est
    jittered["beta"]  <- abs(mle_est["beta"]  + rnorm(1, 0, rw_sd_ch * 0.05))
    jittered["gamma"] <- abs(mle_est["gamma"] + rnorm(1, 0, rw_sd_ch * 0.05))
    jittered["rho"]   <- pmin(pmax(
      mle_est["rho"] + rnorm(1, 0, rw_sd_ch * 0.05), 1e-6), 1 - 1e-6)
    
    candidates <- list(
      list(label = "jittered",    params = jittered),
      list(label = "MLE",         params = mle_est),
      list(label = "true values", params = true_values)
    )
    
    attempts <- list()
    success  <- NULL
    
    for (cand in candidates) {
      cat(sprintf("  Trying pmcmc() with %s starting values (beta=%.6f, gamma=%.6f, rho=%.6f)...\n",
                  cand$label, cand$params["beta"], cand$params["gamma"], cand$params["rho"]))
      
      attempt <- .try_pmcmc(sir_model, cand$params, cand$label, rw_sd_ch)
      attempts[[length(attempts) + 1]] <- attempt
      
      if (attempt$ok) {
        cat(sprintf("  pmcmc() SUCCEEDED with %s starting values (%.2f sec).\n",
                    cand$label, attempt$elapsed))
        success <- attempt
        chain_start_label[ch]  <- cand$label
        chain_start_params[[ch]] <- cand$params
        break
      } else {
        cat(sprintf("  pmcmc() FAILED with %s starting values (%.2f sec): %s\n",
                    cand$label, attempt$elapsed, attempt$msg))
      }
    }
    
    chain_attempt_log[[ch]] <- attempts
    total_elapsed <- sum(vapply(attempts, function(a) a$elapsed, numeric(1)))
    chain_times[ch] <- total_elapsed
    
    if (is.null(success)) {
      cat(sprintf("  Chain %d: all 3 starting points failed; skipping.\n", ch))
      chain_fail_reason[[ch]] <- vapply(attempts, function(a)
        sprintf("%s (%.2f sec): %s", a$label, a$elapsed, a$msg), character(1))
      next
    }
    
    pmcmc_out <- success$pmcmc_out
    
    # ── Extract draws ──────────────────────────────────────────────────────
    raw_traces <- as.data.frame(pmcmc_out@traces)
    n_rows     <- nrow(raw_traces)
    n_burnin   <- floor(n_iter * burn_in)
    
    draws_full <- raw_traces[, required_params, drop = FALSE]
    draws_post <- draws_full[(n_burnin + 1):n_rows, , drop = FALSE]
    ll_full    <- raw_traces[, "loglik"]
    ll_post    <- ll_full[(n_burnin + 1):n_rows]
    
    accept_rate <- mean(diff(draws_full[, "beta"]) != 0)
    
    chain_draws[[ch]]      <- draws_post
    chain_draws_full[[ch]] <- draws_full
    chain_ll[[ch]]         <- ll_post
    chain_ll_full[[ch]]    <- ll_full
    chain_accept[ch]       <- accept_rate
    
    cat(sprintf("  Started from      : %s\n", chain_start_label[ch]))
    cat(sprintf("  Acceptance rate   = %.3f\n", accept_rate))
    cat(sprintf("  Post-burnin draws : %d\n",   nrow(draws_post)))
    cat(sprintf("  Time = %.2f sec (%.2f min)\n", total_elapsed, total_elapsed / 60))
  }
  
  # ── Pool chains ─────────────────────────────────────────────────────────────
  valid_chains <- which(!sapply(chain_draws, is.null))
  n_valid      <- length(valid_chains)
  cat(sprintf("\n=== %d/%d chains completed successfully ===\n", n_valid, n_chains))
  
  total_time_sec <- sum(chain_times, na.rm = TRUE)
  
  # Common meta block used in every return path
  base_meta <- list(
    N_init             = N_init,
    mle_start          = mle_est,
    true_values        = true_values,
    rw_sd              = rw_sd,
    n_iter             = n_iter,
    Np                 = Np,
    n_chains           = n_chains,
    burn_in            = burn_in,
    seed               = seed,
    chain_seeds        = chain_seeds,
    chain_times        = chain_times,        # always preserved now, even on total failure
    chain_accept       = chain_accept,
    chain_start_label  = chain_start_label,   # "jittered" / "MLE" / "true values" / NA
    chain_start_params = chain_start_params,
    chain_fail_reason  = chain_fail_reason,   # per-chain vector of attempt error messages
    chain_attempt_log  = chain_attempt_log,   # full timing/label/ok detail per attempt
    total_time_sec     = total_time_sec,
    valid_chains       = valid_chains
  )
  
  if (n_valid == 0) {
    cat("No valid chains. Both chains failed.\n")
    return(list(
      status       = "FAILED",
      message      = "all chains failed for all starting-point attempts",
      summary      = NULL,
      pooled_draws = NULL,
      chain_draws  = NULL,
      chain_draws_full = NULL,
      chain_ll     = NULL,
      chain_ll_full    = NULL,
      meta = base_meta
    ))
  }
  
  pooled_draws  <- do.call(rbind, chain_draws[valid_chains])
  n_post_burnin <- floor(n_iter * burn_in)
  n_post_kept   <- n_iter - n_post_burnin
  n_post_total  <- nrow(pooled_draws)
  
  # ── Posterior summaries ───────────────────────────────────────────────────
  cat("\n=== POSTERIOR SUMMARY ===\n")
  
  summary_df <- do.call(rbind, lapply(required_params, function(par) {
    x       <- pooled_draws[[par]]
    mn      <- mean(x)
    sd_     <- sd(x)
    q       <- quantile(x, c(0.025, 0.975))
    covered <- true_values[par] >= q[1] && true_values[par] <= q[2]
    
    ess_val     <- tryCatch(mcmcse::ess(x), error = function(e) NA_real_)
    ess_per_sec <- if (is.finite(ess_val) && total_time_sec > 0)
      ess_val / total_time_sec else NA_real_
    
    data.frame(
      parameter = par,
      post_mean = mn,
      post_sd   = sd_,
      q2.5      = q[1],
      q97.5     = q[2],
      true_val  = true_values[par],
      covered   = covered,
      ESS       = ess_val,
      ESS_per_s = ess_per_sec,
      stringsAsFactors = FALSE
    )
  }))
  
  cat(sprintf("%-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s  %-10s  %-12s\n",
              "Param", "Post.Mean", "Post.SD", "Q2.5", "Q97.5",
              "TrueVal", "Covered", "ESS", "ESS/s"))
  cat(strrep("-", 95), "\n")
  for (k in seq_len(nrow(summary_df))) {
    r <- summary_df[k, ]
    cat(sprintf(
      "%-8s  %-10.6f  %-10.6f  %-10.6f  %-10.6f  %-10.6f  %-8s  %-10.2f  %-12.4f\n",
      r$parameter, r$post_mean, r$post_sd,
      r$q2.5, r$q97.5, r$true_val,
      ifelse(r$covered, "YES", "NO"),
      r$ESS, r$ESS_per_s
    ))
  }
  
  cat("\n=== TIMING SUMMARY ===\n")
  cat(sprintf("  Total time    : %.2f sec (%.2f min)\n",
              total_time_sec, total_time_sec / 60))
  cat(sprintf("  Valid chains  : %d / %d\n", n_valid, n_chains))
  for (ch in seq_len(n_chains)) {
    if (ch %in% valid_chains) {
      cat(sprintf("  Chain %d       : %.2f sec  |  accept = %.3f  |  start = %s\n",
                  ch, chain_times[ch], chain_accept[ch], chain_start_label[ch]))
    } else {
      cat(sprintf("  Chain %d       : FAILED  (%.2f sec across attempts)\n",
                  ch, chain_times[ch]))
    }
  }
  
  mean_accept <- mean(chain_accept[valid_chains], na.rm = TRUE)
  if (mean_accept < 0.05 || mean_accept > 0.50)
    cat(sprintf(
      "\n[!] Mean acceptance rate = %.3f (ideal range 0.05-0.30). Consider tuning rw_sd.\n",
      mean_accept
    ))
  
  # ── Return ──────────────────────────────────────────────────────────────────
  list(
    status       = if (n_valid == n_chains) "OK" else "PARTIAL",
    message      = if (n_valid == n_chains) "all chains succeeded"
    else sprintf("%d/%d chains succeeded", n_valid, n_chains),
    summary          = summary_df,
    pooled_draws     = pooled_draws,
    chain_draws      = chain_draws[valid_chains],
    chain_draws_full = chain_draws_full[valid_chains],
    chain_ll         = chain_ll[valid_chains],
    chain_ll_full    = chain_ll_full[valid_chains],
    meta = c(base_meta, list(
      n_post_burnin = n_post_burnin,
      n_post_kept   = n_post_kept,
      n_post_total  = n_post_total
    ))
  )
}

############## Trial fitting

set.seed(1234)
obs_data = data_processing(n = 1e4, rho = 0.05, beta = 2, gamma = 0.5,T.max_analysis = 10)
mif_result = run_mif2_repeated(obs_data = obs_data,
                         N_init = 10000,
                         start_params = c(beta = 1.8, gamma = 0.4, rho = 0.04),         
                         true_values  = c(beta = 2, gamma = 0.5, rho = 0.05),  
                         Nmif         = 400,
                         Np_mif       = 3000,
                         Np_eval      = 2000,
                         n_eval       = 10,
                         n_rep        = 2) ## Repeated 10 times for our run

## Starting something random like, start_params = c(beta= 0.8, gamma=0.4, rho=0.01), will not give good estimates.


####### Running pmcmc

pmcmc_result = run_pmcmc_sir(
    obs_data = obs_data,        # Data frame with two columns: days and incidence count
    N_init = 1e4,          # population size (integer);
    mle_obj = mif_result,         # MLE estimate, with best_obj as the MLE; see run_mif2_repeated
    true_values = c(beta = 2, gamma = 0.5, rho = 0.05),     # named vector: c(beta=, gamma=, rho=)  
    n_iter   = 10,  # For trial run, true run, 1e5
    Np       = 2000, # number of particles for the particle filter
    rw_sd    = 0.007,  # random-walk SD (scalar or per-chain vector)
    n_chains = 2,    # number of independent chains
    burn_in  = 0.5,  # fraction of iterations discarded as burn-in
    seed     = 1234 )
