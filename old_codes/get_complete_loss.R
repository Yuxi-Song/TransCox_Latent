# dat: data.frame(time, event)
# X: covariate matrix
# beta_mat: p x G
# zeta_list: list of spline coefficients for each class (ζ_g)
# knots_list, degree_list: lists of spline basis settings
# xi0, xi1: softmax membership parameters
# kappa, R_list: optional penalty on 2nd derivative (each R_g is quadratic penalty matrix)
loglik_lcph_spline <- function(time, event, X, beta_mat,
                               zeta_list, knots_list, degree,
                               xi0 = min(time), xi1 = max(time),
                               kappa = c(seq(10,1e+17,length=30)), R_list = NULL,
                               return_post = FALSE) {
  X   <- as.matrix(X)
  xi1 <- as.matrix(xi1)
  n   <- nrow(X)
  G   <- ncol(beta_mat)
  t   <- dat$time
  d   <- dat$event

  ## 1) Softmax class-membership probabilities
  linpred <- sapply(1:G, function(g) xi0[g] + X %*% xi1[, g])
  linpred <- linpred - apply(linpred, 1, max)
  prior   <- exp(linpred)
  prior   <- prior / rowSums(prior)

  ## 2) Per-class log-likelihood
  logL_g <- matrix(NA_real_, n, G)
  penalty <- 0

  for (g in 1:G) {
    knots <- knots_list[[g]]
    deg   <- degree
    zeta  <- zeta_list[[g]]

    M <- splines2::mSpline(t,
                 knots = knots[-c(1, length(knots))],
                 degree = deg,
                 Boundary.knots = range(knots),
                 intercept = TRUE)
    I <- splines2::iSpline(t,
                 knots = knots[-c(1, length(knots))],
                 degree = deg,
                 Boundary.knots = range(knots),
                 intercept = TRUE)

    h0 <- as.numeric(M %*% zeta)             # baseline hazard at t_i
    H0 <- as.numeric(I %*% zeta)             # cumulative baseline hazard
    eta <- as.numeric(X %*% beta_mat[, g])   # linear predictor

    logL_g[, g] <- d * (log(pmax(h0, 1e-12)) + eta) - H0 * exp(eta)

    # penalty term if R_list supplied
    if (!is.null(R_list)) {
      penalty <- penalty + crossprod(zeta, R_list[[g]] %*% zeta)
    }
  }

  ## 3) Combine mixture
  log_mix <- log(prior) + logL_g
  m <- apply(log_mix, 1, max)
  Li <- m + log(rowSums(exp(log_mix - m)))
  ll <- sum(Li) - kappa * penalty

  ## 4) Return
  if (!return_post) return(ll)

  post <- exp(log_mix - Li) # posterior P(c_i=g|...)
  colnames(post) <- paste0("class", 1:G)
  list(loglik = ll, post = post, per_obs_loglik = Li, penalty = kappa * penalty)
}
