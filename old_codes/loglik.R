comp_log_fun <- function(
    time,
    event,
    X,
    G = 2, # Number of latent classes
    beta, # dim = (p + 1) * g
    xi, # dim = (p + 1) * (G - 1), class parameter
    zeta,
    K
) {
  l = 0
  xi1 <- min(time);
  xi3 <- max(time);
  Omega=c(192,-132,24,12,0,
         -132,96,-24,-12,12,
         24,-24,24,-24,24,
         12,-12,-24,96,-132,
         0,12,24,-132,192)
  Omega=matrix(Omega,5,5)/( (xi3-xi1)/2 )^5

  p = ncol(X)
  X_tilda <- cbind(1, X)
  # for(i in 1:length(time)) {
  #   for (class in 1 : G) {
  #     g=exp(pmin(zeta[, class],500))
  #     l=-K*t(g)%*%Omega%*%g
  #     bas_h=as.vector(joint.Cox::M.spline(time,xi1=xi1,xi3=xi3)%*%g)
  #     bas_H=as.vector(joint.Cox::I.spline(time,xi1=xi1,xi3=xi3)%*%g)
  #     Xi_tilda <- X_tilda[i, ]
  #     bXi = as.numeric(as.matrix(Xi_tilda) %*% beta[, class])
  #     xi0 <- rep(0, p+1)
  #     xi <- cbind(xi0, xi)
  #     pig = exp(as.matrix(Xi_tilda) %*% as.matrix(xi[, class])) / sum(exp(as.matrix(Xi_tilda) %*% as.matrix(xi)))
  #     l = l + log(pig * ((bas_h[i] * exp(bXi)) ^ (event[i])) * exp(-exp(bXi) * bas_H[i]))
  #     }
  #   }
  L = 0
  for (i in 1:length(time)) {
    lik_i = 0
    for (class in 1:G) {
      g = exp(pmin(zeta[, class], 500))
      bas_h = as.vector(M.spline(time, xi1=xi1, xi3=xi3) %*% g)
      bas_H = as.vector(I.spline(time, xi1=xi1, xi3=xi3) %*% g)
      Xi_tilda <- X_tilda[i, ]
      xi0 <- rep(0, p+1)
      xi <- cbind(xi0, xi)
      bXi = as.numeric(as.matrix(Xi_tilda) %*% beta[, class])
      eta_all <- as.matrix(Xi_tilda) %*% xi  # gives 1 x G
      eta_all <- eta_all - apply(eta_all, 1, max)
      pig <- exp(eta_all[class]) / sum(exp(eta_all))
      lik_i = lik_i + pig * ((bas_h[i] * exp(bXi)) ^ event[i]) *
        exp(-exp(bXi) * bas_H[i])
    }
    L = L + log(lik_i)
  }

  pen = 0
  for (class in 1:G) {
    g = exp(pmin(zeta[, class], 500))
    pen = pen - K * t(g) %*% Omega %*% g
  }
  L = L + pen
  return(L)
}
