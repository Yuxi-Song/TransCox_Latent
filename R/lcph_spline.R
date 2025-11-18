# ====================== Packages ======================
# suppressPackageStartupMessages({
#   library(splines2)   # mSpline / iSpline
# })

init_lcph <- function(time, event, X, G = 2, target_p2 = 0.5, use_kmeans = TRUE) {
  X <- as.matrix(X)
  p <- ncol(X)

  # start with desired prevalence
  xi0 <- log(target_p2 / (1 - target_p2))
  xi1 <- matrix(0, nrow = p, ncol = G - 1)

  # provisional classes
  if (use_kmeans) {
    cl <- kmeans(cbind(scale(log(pmax(time, 1e-6))), event),
                 centers = G, nstart = 10)$cluster
    # logistic for xi (class-2 vs class-1)
    logit_fit <- glm(I(cl == 2) ~ X, family = binomial())
    xi0 <- 0.5 * unname(coef(logit_fit)[1])               # shrink a bit
    xi1 <- 0.5 * matrix(unname(coef(logit_fit)[-1]), ncol = 1)
  }

  # beta from Cox
  beta0 <- matrix(0, nrow = p, ncol = G)
  if (use_kmeans) {
    for (g in 1:G) {
      idg <- (cl == g)
      fitg <- survival::coxph(survival::Surv(time[idg], event[idg]) ~ X[idg, , drop = FALSE])
      beta0[, g] <- unname(coef(fitg))
    }
  } else {
    fit_all <- survival::coxph(survival::Surv(time, event) ~ X)
    beta_all <- unname(coef(fit_all))
    beta0 <- cbind(beta_all, beta_all)
  }
  # small jitter/shrink
  beta0 <- 0.9 * beta0 + matrix(rnorm(length(beta0), 0, 0.05), nrow = p)

  list(beta0 = beta0, xi0 = xi0, xi1 = xi1)
}

# ------------------ helpers ------------------
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# build M/I bases and roughness penalty Ω ≈ ∫ (f''(t))^2 dt
make_bases <- function(time, K = 5, degree = 3) {
  tmin <- min(time); tmax <- max(time)
  n_int <- max(0, K - (degree + 1))
  knots <- if (n_int > 0)
    as.numeric(quantile(time, probs = seq_len(n_int)/(n_int+1)))
  else NULL

  M <- splines2::mSpline(time, degree = degree, intercept = TRUE,
                         knots = knots, Boundary.knots = c(tmin, tmax))
  I <- splines2::iSpline(time, degree = degree, intercept = TRUE,
                         knots = knots, Boundary.knots = c(tmin, tmax))

  grid <- seq(tmin, tmax, length.out = 401)
  D2   <- splines2::mSpline(grid, degree = degree, intercept = TRUE,
                  knots = knots, Boundary.knots = c(tmin, tmax), derivs = 2)
  dt   <- (tmax - tmin)/(length(grid) - 1)
  Omega <- t(D2) %*% D2 * dt
  attr(M, "knots") <- knots; attr(M, "Boundary.knots") <- c(tmin, tmax)
  list(M = as.matrix(M), I = as.matrix(I), Omega = as.matrix(Omega),
       knots = knots, bknots = c(tmin, tmax))
}

# ------------------ E step ------------------
e_step <- function(X, time, event, beta, theta, xi0, xi1, M, I) {
  n <- nrow(X); G <- ncol(beta)
  w <- matrix(NA_real_, n, G)

  # π_ig (class 1 reference; only g=2 has params)
  eta_mat <- cbind(0, drop(X %*% xi1) + as.numeric(xi0))   # n x 2
  logpi   <- cbind(-log1p(exp(eta_mat[,2])),                 # log π1
                   eta_mat[,2] - log1p(exp(eta_mat[,2])))    # log π2

  for (g in 1:G) {
    bg <- exp(theta[, g])             # positive spline coefs
    h0 <- drop(M %*% bg);  h0 <- pmax(h0, 1e-12)
    H0 <- drop(I %*% bg)
    eta <- drop(X %*% beta[, g])
    w[, g] <- logpi[, g] + event * (log(h0) + eta) - exp(eta) * H0
  }
  w <- exp(w - apply(w, 1, logsumexp)) # normalize rows
  w
}

# ------------------ M step (per class) ------------------
mstep_class <- function(X, time, event, w_g, beta_g, theta_g, M, I, Omega, pen = 1e-2) {
  p <- ncol(X); K <- length(theta_g)

  obj_grad <- function(par) {
    beta <- par[1:p]
    b    <- par[(p+1):(p+K)]
    bg   <- exp(b)

    Mb <- drop(M %*% bg); Mb <- pmax(Mb, 1e-12)
    Ib <- drop(I %*% bg)
    eta <- drop(X %*% beta)
    mu  <- exp(eta) * Ib

    ll <- sum(w_g * (event * (log(Mb) + eta) - mu))
    pen_term <- -0.5 * pen * drop(t(bg) %*% Omega %*% bg)

    score_beta <- colSums(w_g * (event - mu) * X)

    grad_bg_linear <- colSums(w_g * (event * (M / Mb) - (exp(eta) * I))) -
      pen * drop(Omega %*% bg)
    score_b <- bg * grad_bg_linear

    list(obj = -(ll + pen_term), grad = -c(score_beta, score_b))
  }

  par0 <- c(beta_g, theta_g)
  opt <- optim(par0,
               fn = function(par) obj_grad(par)$obj,
               gr = function(par) obj_grad(par)$grad,
               method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-7))
  list(beta = opt$par[1:p], theta = opt$par[(p+1):(p+K)], value = opt$value)
}

# ------------------ M step (membership, G=2) ------------------
mstep_xi <- function(X, w) {
  p <- ncol(X)
  obj_grad <- function(par) {
    xi0 <- par[1]; xi1 <- par[1 + 1:p]
    eta2 <- xi0 + drop(X %*% xi1)
    logpi1 <- -log1p(exp(eta2))
    logpi2 <-  eta2 - log1p(exp(eta2))

    ll <- sum(w[,1] * logpi1 + w[,2] * logpi2)
    pi2 <- plogis(eta2)
    diff2 <- w[,2] - pi2
    grad <- -c(sum(diff2), colSums(diff2 * X))
    list(obj = -ll, grad = grad)
  }
  par0 <- c(0, rep(0, p))
  opt <- optim(par0,
               fn = function(par) obj_grad(par)$obj,
               gr = function(par) obj_grad(par)$grad,
               method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-8))
  list(xi0 = opt$par[1], xi1 = matrix(opt$par[-1], ncol = 1))
}

# ------------------ EM wrapper ------------------
fit_lcph_spline_em <- function(time, event, X,
                               K = 100, degree = 3,
                               pen = 1e-2,
                               maxit = 10000, tol = 1e-10,
                               verbose = TRUE, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X); G <- 2

  bs <- make_bases(time, K = K, degree = degree)
  M <- bs$M; I <- bs$I; Omega <- bs$Omega

  # light initialization
  w  <- matrix(1/G, n, G)
  beta  <- matrix(0, nrow = p, ncol = G)
  theta <- matrix(log(0.1), nrow = K, ncol = G)   # log-coefs
  # xi0   <- 0.1; xi1 <- matrix(c(0.1, 0.1), nrow = p, ncol = 1)
  Xstd <- scale(X)
  # crude prevalence from a median time split
  z0 <- as.integer(time > median(time) & event == 1)   # 0/1 indicator
  pi2 <- mean(z0)                                      # ~ class-2 prevalence
  xi0 <- qlogis(pmin(pmax(pi2, 1e-3), 1-1e-3))         # clamp for stability
  xi1 <- matrix(-0.1, nrow = ncol(Xstd), ncol = 1)        # start with no X effect
  # inits <- init_lcph(time, event, X, G = 2, target_p2 = 0.5, use_kmeans = TRUE)
  # beta0 <- inits$beta0
  # xi0   <- inits$xi0
  # xi1   <- inits$xi1

  loglik_obs <- function() {
    # log L = sum_i log( sum_g π_ig f_gi )
    eta2 <- xi0 + drop(X %*% xi1)
    logpi <- cbind(-log1p(exp(eta2)), eta2 - log1p(exp(eta2)))
    out <- numeric(n)
    for (i in 1:n) {
      lg <- numeric(G)
      for (g in 1:G) {
        bg <- exp(theta[, g])
        h0 <- max(1e-12, sum(M[i, ] * bg))
        H0 <- sum(I[i, ] * bg)
        eta <- sum(X[i, ] * beta[, g])
        lg[g] <- logpi[i, g] + event[i] * (log(h0) + eta) - exp(eta) * H0
      }
      out[i] <- logsumexp(lg)
    }
    sum(out)
  }

  ll_old <- -Inf
  for (iter in 1:maxit) {
    # E
    w <- e_step(X, time, event, beta, theta, xi0, xi1, M, I)

    # M: membership
    xi_fit <- mstep_xi(X, w)
    xi0 <- xi_fit$xi0; xi1 <- xi_fit$xi1

    # M: per-class (beta, theta)
    for (g in 1:G) {
      upd <- mstep_class(X, time, event, w[, g], beta[, g], theta[, g],
                         M, I, Omega, pen = pen)
      beta[, g]  <- upd$beta
      theta[, g] <- upd$theta
    }

    ll_new <- loglik_obs()
    if (verbose) cat(sprintf("EM iter %03d  logLik = %.6f\n", iter, ll_new))
    if (abs(ll_new - ll_old) < tol * (1 + abs(ll_old))) break
    ll_old <- ll_new
  }

  list(beta = beta, theta = theta, xi0 = xi0, xi1 = xi1,
       w = w, M = M, I = I, Omega = Omega,
       knots = attr(M, "knots"), bknots = attr(M, "Boundary.knots"),
       logLik = ll_old, iters = iter, converged = (iter < maxit))
}

# ================== Example ==================
# set.seed(42)
# dat <- simulate_lcph(n = 500)
# X   <- as.matrix(dat[, c("X1","X2")])
#
# fit <- fit_lcph_spline_em(time = dat$time,
#                           event = dat$event,
#                           X = X,
#                           K = 5, degree = 3,
#                           pen = 1e-2, maxit = 10000, tol = 1e-8,
#                           verbose = TRUE)
#
# cat("\nEstimated beta (columns = classes):\n"); print(round(fit$beta, 3))
# cat("\nEstimated xi0, xi1 (class-2 vs class-1):\n"); print(round(c(fit$xi0, t(fit$xi1)), 3))
# cat(sprintf("\nConverged: %s in %d iterations (logLik = %.3f)\n",
#             fit$converged, fit$iters, fit$logLik))
#
# # ----- hard class calls & accuracy -----
# Zhat <- max.col(fit$w)          # MAP class
# acc  <- mean(Zhat == dat$class_true)
# cat(sprintf("\nClass accuracy = %.1f%%\n", 100*acc))
#
# # ----- overlay true vs estimated baseline hazards -----
# tgrid <- seq(min(dat$time), max(dat$time), length.out = 300)
#
# # estimated h0_g(t)
# Mgrid <- mSpline(tgrid, degree = 3, intercept = TRUE,
#                  knots = fit$knots, Boundary.knots = fit$bknots)
# hhat1 <- drop(Mgrid %*% exp(fit$theta[,1]))
# hhat2 <- drop(Mgrid %*% exp(fit$theta[,2]))
#
# # true Weibull h0_g(t) used in simulate_lcph()
# lam <- attr(dat, "lambda"); ksh <- attr(dat, "kshape")
# htrue1 <- (ksh[1]/lam[1]) * (tgrid/lam[1])^(ksh[1]-1)
# htrue2 <- (ksh[2]/lam[2]) * (tgrid/lam[2])^(ksh[2]-1)
#
# par(mfrow = c(1,2))
# plot(tgrid, htrue1, type="l", lwd=2, xlab="time", ylab="hazard",
#      main="Class 1: true vs estimated")
# lines(tgrid, hhat1, lwd=2, lty=2)
# legend("topleft", c("true Weibull", "estimated M-spline"), lwd=2, lty=c(1,2), bty="n")
#
# plot(tgrid, htrue2, type="l", lwd=2, xlab="time", ylab="hazard",
#      main="Class 2: true vs estimated")
# lines(tgrid, hhat2, lwd=2, lty=2)
# legend("topleft", c("true Weibull", "estimated M-spline"), lwd=2, lty=c(1,2), bty="n")
