#' get_C_index: Calculate Harrel's C-index.
#' @description Calculate the Harrel's concordance index(C-index) based on the estimated beta
#' @param cov The covariate matrix.
#' @param delta a vector of the censoring indicator corresponding to time points of dH. 0 indicating censoring, while 1 indicating non-censoring.
#' @param beta The estimated beta vector.
#' @param time Time vector
#' @author Yuxi Song. Email: lucillesong@ucla.edu
#' @examples
#' n <- 100
#' time <- rexp(n, rate = 0.1)  # Survival times
#' status <- sample(0:1, n, replace = TRUE)  # Censoring indicator
#' X <- matrix(rnorm(n * 2), n, 2)  # Two covariates
#' cox_model <- coxph(Surv(time, status) ~ X[,1] + X[,2])
#' beta <- summary(cox_model)$coefficients[, 1]
#' C <- get_C_index(X, beta, time, status)
#' @export

get_C_index <- function(cov, beta, time, delta) {
  risk_scores <- as.matrix(cov) %*% beta
  # Compute C-Index
  C_index <- survcomp::concordance.index(risk_scores, time, delta)$c.index
  return(C_index)
}
