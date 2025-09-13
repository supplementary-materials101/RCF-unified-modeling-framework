require("tidyverse", "ggplot2")

# ---- Read the data ----
library(tidyverse)
df <- read_csv("USA_INVESTMENT_MONTHLY_VW_CAP_1963_2024.csv",
               col_types = cols())

## Y: capped, value-weighted monthly returns of investment-theme portfolios on the U.S. market

## X: monthly series for the five Fama–French benchmark factors: MKT − RF, SMB, HML, RMW, and CMA.

## Anomaly performance metric: model-adjusted intercept (alpha/abnormal return)

# ---- Define sample windows ----
original_start <- 196307
original_end   <- 199407
replication_start  <- 199408
replication_end    <- 202412

# ---- Original-sample variables (7/1963 to 7/1994) ----
y_original <- df %>%
  filter(date >= original_start, date <= original_end) %>%
  pull(ret)

X_original <- df %>%
  filter(date >= original_start, date <= original_end) %>%
  select(Mkt_RF, SMB, HML, RMW, CMA)

# ---- Replication-sample variables (8/1994 to 12/2024) ----
y_replication <- df %>%
  filter(date >= replication_start, date <= replication_end) %>%
  pull(ret)

X_replication <- df %>%
  filter(date >= replication_start, date <= replication_end) %>%
  select(Mkt_RF, SMB, HML, RMW, CMA)

# ---- Posterior Inference ----
log_integrate <- function(f, lower, upper, rel.tol = 1e-12) {
  ## numerically integrate exp(log f) on a log scale
  g <- function(x) exp(f(x))        # recover the original scale for integrate
  stats::integrate(g, lower, upper, rel.tol = rel.tol)$value
}

posterior_null_prob <- function(y, X, pi_null = 0.5, tau = 1) {
  ## 1 — OLS estimates
  y <- as.vector(y);   X <- as.matrix(X)
  Xd <- cbind(1, X)                    # add intercept
  n  <- nrow(Xd);     p  <- ncol(Xd) - 1
  bh <- solve(crossprod(Xd), crossprod(Xd, y))
  resid <- y - Xd %*% bh
  df  <- n - p - 1
  s2  <- sum(resid^2) / df
  Vb  <- s2 * solve(crossprod(Xd))
  b0  <- as.numeric(bh[1]);      se0 <- sqrt(Vb[1,1])
  nu  <- df                      # Student-t df
  
  ## 2 — log-densities
  log_tdens  <- function(b) dt((b - b0)/se0, df = nu, log = TRUE) - log(se0)
  log_prior  <- function(b) log(2) + dnorm(b, sd = tau, log = TRUE)
  
  ## 3 — evidence on each half-line (log-sum-exp inside integrator)
  log_int_fun <- function(side) {
    function(b) log_tdens(b) + log_prior(b)          # log integrand
  }
  E0 <- pi_null * log_integrate(log_int_fun("neg"), -Inf, 0)
  E1 <- (1 - pi_null) * log_integrate(log_int_fun("pos"), 0,  Inf)
  
  ## 4 — posterior null probability
  post_M_null <- E0 / (E0 + E1)
  return(post_M_null)
}

## posterior null probability - original sample
posterior_null_prob(y=y_original, 
                    X=X_original, 
                    pi_null = 0.5,
                    tau=0.0041)
## [1] 5.796415e-16

## posterior null probability - replication sample
posterior_null_prob(y=y_replication, 
                    X=X_replication, 
                    pi_null = 0.7,
                    tau=0.0041)
## [1] 0.08838085

# ---- Sensitivity check for τ ---- 
## 1. Grid of τ values
tau_grid <- seq(from = 0.001, to = 0.005, by = 0.0001)

## 2. Vector of π_null values to examine
pi_null_grid <- c(0.30, 0.50, 0.70)

## 3. Evaluate Pr(M0 | data) for every (π_null, τ) pair
sensitivity_tau <- do.call(
  rbind,
  lapply(pi_null_grid, function(pp) {
    post_vals <- sapply(
      tau_grid,
      function(tt)
        posterior_null_prob(
          y   = y_replication,
          X   = X_replication,
          pi_null = pp,
          tau = tt
        )
    )
    data.frame(
      tau     = tau_grid,
      pi_null     = pp,
      post_M_null = post_vals
    )
  })
)

## 4. Inspect the combined data-frame
head(sensitivity_tau)
##      tau pi_null    post_M_null
## 1 0.0010 0.3 0.03154233
## 2 0.0011 0.3 0.02900435
## ...

## 'sensitivity_tau' now contains one row per τ–π_null combination
## with columns: tau, pi_null, post_M_null

## ---------------------------------------------------------------------------
## Function :  plot_sensitivity_tau
## Purpose  :  Overlay sensitivity curves of Pr(M0 | data) vs τ for
##             three prior weights (π_null = 0.3, 0.5, 0.7)
## ---------------------------------------------------------------------------
plot_sensitivity_tau <- function(df) {
  stopifnot(all(c("tau", "pi_null", "post_M_null") %in% names(df)))
  library(ggplot2)
  
  ## axis breaks (ensure 2.5 % and 1 % appear)
  pretty_y  <- scales::pretty_breaks(n = 6)(df$post_M_null)
  
  y_breaks <- sort(unique(c(scales::pretty_breaks(n = 6)(df$post_M_null),
                            0.05, 0.025, 0.01)))
  
  ggplot(df, aes(x = tau, y = post_M_null, colour = as.factor(pi_null))) +
    geom_line(linewidth = 0.8) +
    ## reference lines ------------------------------------------------------
  geom_hline(yintercept = 0.05,
             linetype   = "dashed",
             colour     = "black") +
    geom_hline(yintercept = 0.025,
               linetype   = "dashed",
               colour     = "black") +
    geom_hline(yintercept = 0.01,
               linetype   = "dashed",
               colour     = "black") +
    ## axes -----------------------------------------------------------------
  scale_y_continuous(
    breaks = y_breaks,
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    ## legend & labels ------------------------------------------------------
  scale_colour_manual(
    values = c("0.3" = "grey70", "0.5" = "grey45", "0.7" = "black"),
    name   = expression(pi[null])
  ) +
    labs(
      x = expression(tau),
      y = expression("Pr(" * M[null] * " | data)")
    ) +
    ## theme ----------------------------------------------------------------
  theme_bw(base_size = 12) +
    theme(
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major.y = element_line(colour = "grey85"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold"),
      legend.position    = "right"
    )
}

Figure1 <- plot_sensitivity_tau(sensitivity_tau)
print(Figure1)
