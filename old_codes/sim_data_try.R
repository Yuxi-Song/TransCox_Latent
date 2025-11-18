# Simulate latent-class proportional hazards data
# with class-specific Weibull baselines (PH parameterization)
# H0_g(t) = (t / lambda_g)^k_g,  h0_g(t) = (k_g / lambda_g) * (t / lambda_g)^(k_g - 1)
# => Inversion:  t_i = lambda_g * [ -log(U_i) * exp(-eta_i) ]^(1/k_g)

simulate_lcph <- function(
    n,                           # number of observations
    G = 2,                       # number of latent classes
    p = 2,                       # number of covariates (intercept not included)
    # ---- class-membership model (multinomial logit) ----
    xi0 = 0,                  # length-G intercepts; xi0[1] = 0 (reference) recommended
    xi1 = matrix(c(0.2, -0.2), ncol = G - 1),                  # p x G slopes; xi1[,1] = 0 (reference) recommended
    # ---- survival model (class-specific PH with Weibull baseline) ----
    beta_mat = matrix(c(2, 2, -2, -2), ncol = 2, byrow = TRUE),             # p x G class-specific log-HR
    lambda_vec = c(0.102, 0.109),           # length-G Weibull scale (>0)
    k_vec = c(4.719, 4.432),                # length-G Weibull shape (>0)
    # ---- censoring ----
    lunif = 0,
    censor_type = c("uniform","exponential")[1],
    seed = 123,
    return_probs = TRUE
){
  set.seed(seed)

  ## 1) Covariates: user can pass p, here we simulate p columns (example: 1 normal, 1 binary, then normals)
  if (p < 1) stop("p must be >= 1")
  X <- as.matrix(cbind(
    rnorm(n, 3, 2),
    rbinom(n, 1, 0.5),
    if (p > 2) replicate(p - 2, rnorm(n)) else NULL
  ))
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  p <- ncol(X)

  # ## 2) Class-membership parameters (soft identifiability defaults)
  # if (is.null(xi0)) {
  #   xi0 <- rep(0, G)               # set all zero; reference achieved by subtracting max in softmax
  #   # You may set xi0[1] <- 0 explicitly if you later constrain columns.
  # }
  # if (is.null(xi1)) {
  #   xi1 <- matrix(0, nrow = p + 1, ncol = G)  # default no covariate effects on class membership
  #   # Optionally: xi1[,1] <- 0 to mark class 1 as reference
  # }
  # stopifnot(length(xi0) == G, all(dim(xi1) == c(p, G)))

  ## 3) Softmax for class probabilities, then sample class Z
  xi0 <- c(0, xi0)
  xi1 <- cbind(0, xi1)
  linpred <- sapply(1:G, function(g) as.numeric(xi0[g] + X %*% (xi1[, g])))
  # stabilise softmax
  linpred_centered <- linpred
  # - apply(linpred, 1, max)
  exp_lp <- exp(linpred_centered)
  pi_mat <- exp_lp / rowSums(exp_lp)
  Z <- matrix(0, nrow=n, ncol = G)
  for (i in 1:n){
    Z[i, ] = t(rmultinom(1,1,pi_mat[i,]))
  }
  Z <- apply(Z, 1, which.max)

  # ## 4) Survival parameters
  # if (is.null(beta_mat))  beta_mat  <- matrix(0.2, nrow = p, ncol = G)   # mild effects by default
  # if (is.null(lambda_vec)) lambda_vec <- rep(2, G)                        # scales
  # if (is.null(k_vec))      k_vec      <- rep(1.5, G)                      # shapes
  # stopifnot(all(dim(beta_mat) == c(p, G)),
  #           length(lambda_vec) == G, length(k_vec) == G,
  #           all(lambda_vec > 0), all(k_vec > 0))

  ## 5) Simulate event times via closed-form inversion
  eta <- sapply(seq_len(n), function(i) crossprod(as.matrix(X[i, ]), beta_mat[, Z[i]]))
  eta <- as.numeric(eta)
  lambda_i <- lambda_vec[Z]
  k_i <-  k_vec[Z]
  scale_i  <- lambda_i * exp(-eta / k_i) # PH adjustment of scale
  T_event  <- rweibull(n, shape = k_i, scale = scale_i)
  ## 6) Censoring
  if (censor_type == "uniform") {
    # c_max <- as.numeric(quantile(T_event, probs = 1 - target_censor_rate))
    # if (!is.finite(c_max) || c_max <= 0) c_max <- max(T_event)
    C <- runif(n, lunif, lunif + 1)
  } else {
    lam <- 1 / (quantile(T_event, 0.7) + 1e-6)
    C <- rexp(n, rate = lam)
  }
  time  <- pmin(T_event, C)
  event <- as.integer(T_event <= C)

  ## 7) Output
  out <- data.frame(
    id = seq_len(n),
    time = time,
    event = event,
    class_true = Z,
    X
  )
  if (return_probs) {
    colnames(pi_mat) <- paste0("pi_class", seq_len(G))
    out <- cbind(out, pi_mat)
  }

  ## Attributes for reproducibility
  attr(out, "beta_mat")  <- beta_mat
  attr(out, "lambda_vec") <- lambda_vec
  attr(out, "k_vec") <- k_vec
  #attr(out, "xi0") <- xi0
  attr(out, "xi1") <- xi1
  attr(out, "G")   <- G
  attr(out, "p")   <- p
  out
}


# set.seed(156)
# dat <- simulate_lcph(
#   200,
#   G = 2,
#   beta_mat = matrix(c(2, -2, 0.1, 0.3, 1, 2), ncol = 2),
#   lambda_vec = c(1, 2),
#   k_vec = c(0.4, 1.1),
#   xi0 = c(0, -0.5),
#   xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2, byrow = FALSE))

param_target <- list(
  n = 200,
  G = 2,
  beta_mat = matrix(c(0.5, -0.5, 0.1, 0.3, 1, 2), ncol = 2),
  lambda_vec = c(1, 2),
  k_vec = c(0.4, 1.1),
  xi0 = c(0, -0.5),
  xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2)
)

# param_source1 <- list(
#   G = 2,
#   beta_mat = matrix(c(0.5, -0.5, -0.1, -2, 1, 2), ncol = 2),
#   lambda_vec = c(1, 5),
#   k_vec = c(0.4, 0.5),
#   xi0 = c(0, -0.5),
#   xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2)
# )

param_source1 <- list(
  G = 2,
  beta_mat = matrix(c(0.5, -0.5, 0, -2, 1, 2), ncol = 2),
  lambda_vec = c(1, 5),
  k_vec = c(0.4, 0.5),
  xi0 = c(0, -0.5),
  xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2)
)

param_source2 <- list(
  G = 2,
  beta_mat = matrix(c(0.5, -0.5, 0.1, 0.3, -1, 2), ncol = 2),
  lambda_vec = c(1, 2),
  k_vec = c(-0.4, 1.1),
  xi0 = c(0, -0.5),
  xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2)
)

param_source3 <- list(
  G = 2,
  beta_mat = matrix(c(0.5, -0.5, -0.1, -0.3, -1, -2), ncol = 2),
  lambda_vec = c(10, 20),
  k_vec = c(1.1, 0.4),
  xi0 = c(0, -0.5),
  xi1 = matrix(c(0, -0.5, 0, 1), ncol = 2)
)

param <- list(
  source1 = param_source1,
  source2 = param_source1,
  source3 = param_source1,
  target  = param_target
)

GetSimData <- function(
    nt = 200,
    ns = c(600, 700, 800),
    num_s = 3,
    param_setting = param) {
  if(length(ns) != num_s) {
    stop("Number of source data num_s must align with the length of vector ns")
  }
  if (length(param_setting) < num_s + 1) {
    stop("param_setting must contain all source and target parameter lists")
  }

  # ---- simulate each source separately ----
  data_list <- list()

  for (i in seq_len(num_s)) {
    par_i <- param_setting[[paste0("source", i)]]
    par_i$n <- ns[i]
    data_list[[paste0("source", i)]] <- do.call(simulate_lcph, par_i)
  }

  # ---- simulate target ----
  par_t <- param_setting$target
  par_t$n <- nt
  data_list$target <- do.call(simulate_lcph, par_t)

  return(data_list)
}
dat <- GetSimData()

