require("tibble", "dplyr")
# ----------------------- simulation settings -----------------------
m_tests    <- 10000    # number of tests
pi0_true   <- 0.60     # true null fraction (fixed)
rho_eq     <- 0.30     # equicorrelation
df_t       <- 360      # degrees of freedom for t
delta_alt  <- 1.50     # noncentral shift for true positives (weak)
seed_sim   <- 123      # random seed

bayes_lambda <- 0.05                 # posterior-null threshold
pi0_priors   <- c(0.30, 0.50, 0.70)  # prior skepticism levels to sweep
tau_prior    <- 1.00                 # prior scale over s.e. (for z-equivalents)

bh_q <- 0.05                         # BH FDR level

out_csv <- "head_to_head_table.csv"

# ----------------------- helpers -----------------------
p_two_sided_t <- function(t, df) 2 * (1 - pt(abs(t), df = df))

storey_pi0 <- function(p, lambdas = c(0.5, 0.8, 0.9)) {
  m <- length(p)
  ests <- sapply(lambdas, function(l) {
    num <- sum(p > l); denom <- (1 - l) * m
    val <- if (denom > 0) num / denom else 1
    min(max(val, 0), 1)
  })
  median(ests)
}

bh_rejections <- function(p, q = 0.05) {
  m <- length(p); o <- order(p); ps <- p[o]
  thr <- q * (1:m) / m
  k <- max(which(ps <= thr), 0)
  rej <- rep(FALSE, m); if (k > 0) rej[o[1:k]] <- TRUE
  tau_star <- if (k > 0) max(ps[1:k]) else NA_real_
  list(reject = rej, tau_star = tau_star)
}

fdr_estimate <- function(p, R, tau, lambdas = c(0.5, 0.8, 0.9)) {
  if (is.na(R) || R <= 0) return(0)     # by convention: if no discoveries, FDR = 0
  if (is.na(tau)) return(NA_real_)
  m <- length(p); pi0_hat <- storey_pi0(p, lambdas = lambdas)
  fdr <- (pi0_hat * m * tau) / R
  min(max(fdr, 0), 1)
}

posterior_null_prob_z <- function(z, pi0, tau) {
  c <- (z * tau) / sqrt(1 + tau^2)
  p_plus <- pnorm(c)
  odds_10 <- ((1 - pi0) / pi0) * (p_plus / pmax(1 - p_plus, .Machine$double.eps))
  1 / (1 + odds_10)
}

bayes_decisions_from_z <- function(z, pi0, tau, lambda_B = 0.05) {
  pnull <- posterior_null_prob_z(z, pi0 = pi0, tau = tau)
  rej <- pnull <= lambda_B
  R <- sum(rej)
  bayes_fdr <- if (R > 0) mean(pnull[rej]) else 0
  list(reject = rej, pnull = pnull, bayes_fdr = bayes_fdr)
}

t3_rule <- function(t) {
  rej <- abs(t) >= 3
  tau_fixed <- 2 * (1 - pnorm(3))  # Normal reference for |t| >= 3
  list(reject = rej, tau_star = tau_fixed)
}

simulate_correlated_t <- function(m, pi0_true, rho, df, delta, seed = 123) {
  set.seed(seed)
  m1 <- round(m * (1 - pi0_true))
  H1_idx <- if (m1 > 0) sample.int(m, m1) else integer(0)
  Z0 <- rnorm(1); Zi <- rnorm(m)
  G  <- sqrt(rho) * Z0 + sqrt(1 - rho) * Zi
  mu <- rep(0, m); if (m1 > 0) mu[H1_idx] <- delta
  Wi <- rchisq(m, df = df)
  T  <- (mu + G) / sqrt(Wi / df)
  list(t = T, H1 = (1:m) %in% H1_idx)
}

head_to_head_table <- function(t, df, pi0_prior, tau_prior, bayes_lambda, bh_q) {
  m <- length(t)
  # p-values (t distribution under null)
  p_ts <- p_two_sided_t(t, df = df)
  
  # Signed z equivalents for Bayesian step
  z_abs <- qnorm(1 - p_ts / 2)
  z_signed <- sign(t) * z_abs
  
  # Bayesian
  B <- bayes_decisions_from_z(z_signed, pi0 = pi0_prior, tau = tau_prior, lambda_B = bayes_lambda)
  rej_B <- B$reject; R_B <- sum(rej_B); fdr_B <- B$bayes_fdr
  
  # 3-sigma |t|
  T3 <- t3_rule(t); rej_T3 <- T3$reject; R_T3 <- sum(rej_T3)
  fdr_T3 <- fdr_estimate(p_ts, R_T3, T3$tau_star)
  
  # BH
  BH <- bh_rejections(p_ts, q = bh_q); rej_BH <- BH$reject; R_BH <- sum(rej_BH)
  fdr_BH <- fdr_estimate(p_ts, R_BH, BH$tau_star)
  
  # Overlaps with Bayesian
  overlap_T3 <- sum(rej_T3 & rej_B)
  overlap_BH <- sum(rej_BH & rej_B)
  
  tibble::tibble(
    Prior_pi0 = sprintf("%.2f", pi0_prior),
    Procedure = c("Bayesian", "3-sigma |t|", "BH"),
   # ThresholdRule = c(
   #   sprintf("P(null|data) <= %.3f", bayes_lambda),
   #  "|t| >= 3",
   #   sprintf("q = %.2f", bh_q)
   # ),
    Tests_m = m,
    Discoveries_R = c(R_B, R_T3, R_BH),
    Estimated_FDR = round(c(fdr_B, fdr_T3, fdr_BH), 4),
    Overlap_with_Bayesian = c(NA_integer_, overlap_T3, overlap_BH)
  )
}

# ----------------------- simulate once, run all priors -----------------------
sim <- simulate_correlated_t(
  m   = m_tests,
  pi0_true = pi0_true,
  rho = rho_eq,
  df  = df_t,
  delta = delta_alt,
  seed = seed_sim
)

tables <- lapply(pi0_priors, function(p0) {
  head_to_head_table(
    t = sim$t,
    df = df_t,
    pi0_prior = p0,
    tau_prior = tau_prior,
    bayes_lambda = bayes_lambda,
    bh_q = bh_q
  )
})
tab_all <- dplyr::bind_rows(tables)

message("Head-to-head (simulation, delta_alt = 1.50; pi0_true = 0.60).")
print(tab_all, n = nrow(tab_all))
tryCatch(write.csv(tab_all, out_csv, row.names = FALSE),
         error = function(e) warning("Could not write output CSV: ", conditionMessage(e)))
message("Wrote combined table to: ", out_csv)

# Quick summary of the simulated mixture
R_true <- sum(sim$H1)
message(sprintf("\nSimulation summary: m = %d, true positives = %d (%.1f%%), rho = %.2f, df = %d, delta = %.2f",
                m_tests, R_true, 100*(1 - pi0_true), rho_eq, df_t, delta_alt))
