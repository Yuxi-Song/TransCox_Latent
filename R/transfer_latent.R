## ==========================================================
## Transfer learning for Cox PH with fixed Bernstein basis
##  - Known source params (beta_s, zeta_s)
##  - True deviations (eta_true, zeta_true)
##  - Simulate TARGET data from beta_t_true, zeta_t_true
##  - Estimate eta, zeta via penalized MLE (ridge)
##  - Compare estimated vs true beta & baseline hazard
## ==========================================================

## ----- 4) Transfer objective: estimate (eta, zeta) on target -----
negloglik_transfer <- function(par, X, time, status,
                               beta_s, zeta_s, m, tau,
                               lambda_eta = 0, lambda_zeta = 0){
  p <- ncol(X)
  eta   <- par[1:p]
  zeta <- par[(p+1):(p + (m+1))]
  beta_t <- beta_s + eta
  zeta_t <- zeta_s + zeta

  # vectorized hazard & integrand for integrate()
  M <- splines2::mSpline(time, degree = degree, intercept = TRUE,
                         knots = knots, Boundary.knots = c(tmin, tmax))
  I <- splines2::iSpline(time, degree = degree, intercept = TRUE,
                         knots = knots, Boundary.knots = c(tmin, tmax))
  h0 <- drop(M %*% zeta_t); h0 <- pmax(h0, 1e-12)
  H0 <- drop(I %*% zeta_t)
  # eta <- drop(X %*% beta_t)

  xb <- drop(X %*% beta_t)

  ll <- sum(-status * (log(h0) + xb) + exp(xb) * H0)

  # Ridge penalties (set to 0 for pure MLE)
  pen <- lambda_eta * sum(eta^2) + lambda_zeta * sum(zeta^2)
  ll + pen
}

# Solve (start from "no deviation")
init <- rep(0, p + (m+1))
fit <- optim(init, negloglik_transfer,
             X = X_t, time = time_t, status = status_t,
             beta_s = beta_s, zeta_s = zeta_s, m = m, tau = tau,
             lambda_eta = 0.01, lambda_zeta = 0.01,        # small ridge helps stability
             method = "BFGS", control = list(maxit = 400, reltol = 1e-8))

eta_hat   <- fit$par[1:p]
zeta_hat <- fit$par[(p+1):(p + (m+1))]
beta_t_hat <- beta_s + eta_hat
zeta_t_hat <- zeta_s + zeta_hat

## ----- 5) Numeric comparison: β and baseline hazard -----
cat("\nTrue beta_t:      ", round(beta_t_true, 3), "\n")
cat("Estimated beta_t: ", round(beta_t_hat, 3), "\n")
cat("Eta (true vs est):\n")
print(rbind(true = round(eta_true, 3), est = round(eta_hat, 3)))

# Hazard comparison on a display grid
tgrid <- seq(0, tau, length.out = 300)
haz_true <- haz_from_alpha(tgrid, zeta_t_true, tau, m)
haz_hat  <- haz_from_alpha(tgrid, zeta_t_hat,  tau, m)

rmse_haz  <- sqrt(mean( (haz_hat - haz_true)^2 ))
mae_beta  <- mean(abs(beta_t_hat - beta_t_true))

cat("\nRMSE (baseline hazard on grid) =", round(rmse_haz, 5), "\n")
cat("MAE  (beta)                     =", round(mae_beta, 5), "\n")

## ----- 6) Plots: estimated vs true (hazard & cumulative) -----
Lambda_true <- cumsum(haz_true) * (tau/length(tgrid))
Lambda_hat  <- cumsum(haz_hat)  * (tau/length(tgrid))

par(mfrow = c(1,2))
plot(tgrid, haz_true, type="l", col="red", lwd=2,
     main="Baseline hazard: True vs Estimated",
     xlab="Time", ylab=expression(lambda[0](t)))
lines(tgrid, haz_hat, col="blue", lwd=2)
legend("topleft", c("True","Estimated"), col=c("red","blue"), lty=1, bty="n")

plot(tgrid, Lambda_true, type="l", col="red", lwd=2,
     main="Cumulative baseline hazard",
     xlab="Time", ylab=expression(Lambda[0](t)))
lines(tgrid, Lambda_hat, col="blue", lwd=2)
legend("topleft", c("True","Estimated"), col=c("red","blue"), lty=1, bty="n")
