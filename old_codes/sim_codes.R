library(joint.Cox)
library(parallel)

nrep <- 100

results <- mclapply(1:nrep, function(seed) {
  set.seed(seed)
  # --- simulate data ---
  dat <- simulate_lcph_spline(
    n = 100,
    beta_mat = matrix(c(-2.1, 4.1), ncol = 1),
    alpha_mat = matrix(c(0.2,0.2,0.05,0.2,0.35), nrow = 1),
    G = 1
  )

  dat1 <- dat[dat$class_true == 1, ]
  t.event <- dat1$time
  event   <- dat1$event
  Z       <- dat1[, c("X1", "X2")]

  # --- fit model ---
  fit <- tryCatch(
    splineCox.reg(
      t.event, event, Z,
      kappa = seq(10, 1e17, length = 30),
      LCV.plot = FALSE
    ),
    error = function(e) NULL
  )

  # --- extract coefficients ---
  if (!is.null(fit)) {
    return(fit$beta)
  } else {
    return(rep(NA, ncol(Z)))
  }

}, mc.cores = detectCores() - 1)

# --- Combine results ---
beta_mat <- do.call(rbind, results)
beta_mat <- beta_mat[complete.cases(beta_mat), , drop = FALSE]

colMeans(beta_mat)
