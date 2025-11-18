# Generate 1000 seeds
sim_num <- 500
seeds <- sample.int(1e9, sim_num)
# beta <- list()
# beta[[1]] <- matrix(0, nrow = sim_num, ncol = 2)
# beta[[2]] <- matrix(0, nrow = sim_num, ncol = 2)
#
# theta <- list()
# theta[[1]] <- matrix(0, nrow = sim_num, ncol = 5)
# theta[[2]] <- matrix(0, nrow = sim_num, ncol = 5)
# converge <- vector("numeric", length = sim_num)
# xi0 <- vector("numeric", length = sim_num)
# xi1 <- matrix(0, nrow = sim_num, ncol = 2)
# acc <- vector("numeric", length = sim_num)

# for(j in seq_along(seeds)) {
#   dat <- simulate_lcph(n = 2000, seed = seeds[j]) # Simulate 2000 obs
#   X   <- as.matrix(dat[, c("X1","X2")])
#   fit <- fit_lcph_spline_em(time = dat$time,
#                             event = dat$event,
#                             X = X,
#                             K = 5, degree = 3,
#                             pen = 1e-2, maxit = 10000, tol = 1e-8,
#                             verbose = TRUE)
#   beta[[1]][j, ] <- fit$beta[, 1]
#   beta[[2]][j, ] <- fit$beta[, 2]
#   theta[[1]][j, ] <- fit$theta[, 1]
#   theta[[1]][j, ] <- fit$theta[, 2]
#   converge[j] <- fit$converged
#   xi0[j] <- fit$xi0
#   xi1[j, ] <- fit$xi1
#   Zhat <- max.col(fit$w)          # MAP class
#   acc[j]  <- mean(Zhat == dat$class_true)
# }


library(doParallel)
# devtools::load_all()

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

results <- foreach(j = seq_along(seeds)) %dopar% {

  dat <- simulate_lcph(n = 2000, seed = seeds[j])
  X   <- as.matrix(dat[, c("X1","X2")])

  fit <- fit_lcph_spline_em(
    time = dat$time,
    event = dat$event,
    X = X,
    K = 5, degree = 3,
    pen = 1e-2,
    maxit = 10000, tol = 1e-8,
    verbose = FALSE
  )

  Zhat <- max.col(fit$w)

  list(
    beta1   = fit$beta[,1],
    beta2   = fit$beta[,2],
    theta1  = fit$theta[,1],
    theta2  = fit$theta[,2],
    xi0     = fit$xi0,
    xi1     = fit$xi1,
    acc     = mean(Zhat == dat$class_true),
    conv    = fit$converged
  )
}

# allocate
beta = list(
  matrix(NA, length(seeds), ncol(fit$beta)),
  matrix(NA, length(seeds), ncol(fit$beta))
)
theta = list(
  matrix(NA, length(seeds), ncol(fit$theta)),
  matrix(NA, length(seeds), ncol(fit$theta))
)
xi0 = numeric(length(seeds))
xi1 = matrix(NA, length(seeds), ncol(fit$xi1))
acc = numeric(length(seeds))
converge = logical(length(seeds))

# fill them
for (j in seq_along(seeds)) {
  beta[[1]][j,]  <- results[[j]]$beta1
  beta[[2]][j,]  <- results[[j]]$beta2
  theta[[1]][j,] <- results[[j]]$theta1
  theta[[2]][j,] <- results[[j]]$theta2
  xi0[j]         <- results[[j]]$xi0
  xi1[j,]        <- results[[j]]$xi1
  acc[j]         <- results[[j]]$acc
  converge[j]    <- results[[j]]$conv
}

stopCluster(cl)

beta_est_1_res <- colMeans(sweep(beta[[1]], 2, c(-0.5, 0.5), "-"))
beta_est_2_res <- colMeans(sweep(beta[[2]], 2, c(2, 4), "-"))

acc_mean_est <- mean(acc)
xi0_est <- mean(xi0 - 2)
xi1_est <- colMeans(sweep(xi1, 2, c(-1.0, 0.5), "-"))

theta_est_1 <- colMeans(theta[[1]])
theta_est_2 <- colMeans(theta[[2]])

# ----- overlay true vs estimated baseline hazards -----
tgrid <- seq(min(dat$time), max(dat$time), length.out = 300)

# estimated h0_g(t)
Mgrid <- mSpline(tgrid, degree = 3, intercept = TRUE,
                 knots = fit$knots, Boundary.knots = fit$bknots)
hhat1 <- drop(Mgrid %*% exp(theta_est_1))
hhat2 <- drop(Mgrid %*% exp(theta_est_2))

# true Weibull h0_g(t) used in simulate_lcph()
lam <- attr(dat, "lambda"); ksh <- attr(dat, "kshape")
htrue1 <- (ksh[1]/lam[1]) * (tgrid/lam[1])^(ksh[1]-1)
htrue2 <- (ksh[2]/lam[2]) * (tgrid/lam[2])^(ksh[2]-1)

par(mfrow = c(1,2))
plot(tgrid, htrue1, type="l", lwd=2, xlab="time", ylab="hazard",
     main="Class 1: true vs estimated")
lines(tgrid, hhat1, lwd=2, lty=2)
legend("topleft", c("true Weibull", "estimated M-spline"), lwd=2, lty=c(1,2), bty="n")

plot(tgrid, htrue2, type="l", lwd=2, xlab="time", ylab="hazard",
     main="Class 2: true vs estimated")
lines(tgrid, hhat2, lwd=2, lty=2)
legend("topleft", c("true Weibull", "estimated M-spline"), lwd=2, lty=c(1,2), bty="n")

