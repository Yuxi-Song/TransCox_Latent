# ------------------ simulation ------------------
simulate_lcph <- function(n,
                          p = 2,
                          beta_mat = matrix(c( -0.5,  2,
                                               0.5, 4), nrow = p, byrow = TRUE),
                          G = 2,
                          # class-2 vs class-1 membership
                          xi0 =  2.0,                  # intercept
                          xi1 = matrix(c(-1.0, 0.5), nrow = p),  # slopes
                          # true Weibull baselines per class (scale λ, shape k)
                          lambda = c(0.8, 1.2), # length-G Weibull scale (>0)
                          kshape = c(2.0, 3.0), # length-G Weibull shape (>0)
                          censor_lo = 0.0, # Lower bound of the censoring distribution
                          censor_hi = 1.2, # Upper bound of the censoring distribution
                          seed = 6) {
  set.seed(seed)
  X <- cbind(rnorm(n), rbinom(n, 1, 0.5))
  colnames(X) <- paste0("X", 1:ncol(X))

  # class probs (G=2, class 1 is reference)
  eta2 <- as.numeric(xi0 + X %*% xi1)       # length n
  pi2  <- plogis(eta2);  pi1 <- 1 - pi2
  Z    <- 1L + rbinom(n, 1, pi2)            # 1 or 2

  # event times: t = λ * (-log U / exp(η))^(1/k)
  U   <- runif(n)
  eta <- vapply(seq_len(n), function(i) drop(X[i, ] %*% beta_mat[, Z[i]]), 0.0)
  lam_i <- lambda[Z]; k_i <- kshape[Z]
  T_event <- lam_i * (-log(U) / exp(eta))^(1 / k_i)

  # independent administrative censoring
  C <- runif(n, censor_lo, censor_hi)
  time  <- pmin(T_event, C)
  event <- as.integer(T_event <= C)

  out <- data.frame(id = seq_len(n), time, event, class_true = Z, X)
  attr(out, "beta_true") <- beta_mat
  attr(out, "xi0_true")  <- xi0
  attr(out, "xi1_true")  <- xi1
  attr(out, "lambda")    <- lambda
  attr(out, "kshape")    <- kshape
  out
}
