pre_whiten <- function(variance_matrix, outcome) {
  eigen_decomposition <- eigen(variance_matrix)

  if (min(eigen_decomposition$values) <= 0) {
    stop("Variance matrix has zero or negative eigenvalues.", call. = FALSE)
  }

  inverse_square_root <-
    eigen_decomposition$vectors %*%
    diag(1 / sqrt(eigen_decomposition$values)) %*%
    t(eigen_decomposition$vectors)

  as.numeric(inverse_square_root %*% outcome)
}

constrain_variance_params <- function(params, eigenvalues, margin = 0.01) {
  constraint <- 0
  for (param_index in 2:length(params)) {
    constraint <- constraint + params[param_index] * eigenvalues^(param_index - 1)
  }
  constraint <- min(constraint) + params[1]

  if (constraint <= 0) {
    params <- c(params[1] - constraint + margin, params[-1])
  }

  params
}

make_mle_objective <- function(adjacency) {
  eigen_decomposition <- eigen(adjacency)
  eigenvalues <- eigen_decomposition$values
  eigenvectors <- eigen_decomposition$vectors

  function(params, y) {
    sig0 <- params[1]
    sig1 <- params[2]
    sig2 <- params[3]
    z <- t(eigenvectors) %*% y

    variance_values <- sig1 * eigenvalues + sig2 * eigenvalues^2 + sig0
    min_variance <- min(variance_values)
    if (min_variance <= 0) {
      sig0 <- sig0 - min_variance + 0.01
      variance_values <- sig1 * eigenvalues + sig2 * eigenvalues^2 + sig0
    }

    sum(log(abs(variance_values))) + t(z) %*% diag(1 / variance_values) %*% z
  }
}

estimate_mle_prewhiten_vectors <- function(adjacency, x, y, initial_values = NULL) {
  popn <- nrow(adjacency)
  eigenvalues <- eigen(adjacency)$values
  objective <- make_mle_objective(adjacency)

  if (is.null(initial_values)) {
    initial_values <- runif(3, c(0.01, 0.01, 0.01), max = c(1, 0.5, 0.5))
    initial_values <- constrain_variance_params(initial_values, eigenvalues)
  }

  fit_x <- optim(
    par = initial_values,
    fn = objective,
    method = "L-BFGS-B",
    y = x,
    lower = c(0.001, 0.001, 0.01),
    upper = c(Inf, Inf, Inf)
  )

  fit_y <- optim(
    par = initial_values,
    fn = objective,
    method = "L-BFGS-B",
    y = y,
    lower = c(0.001, 0.001, 0.01),
    upper = c(Inf, Inf, Inf)
  )

  coef_x <- constrain_variance_params(fit_x$par, eigenvalues)
  coef_y <- constrain_variance_params(fit_y$par, eigenvalues)

  cov_x <- coef_x[1] * diag(popn) + coef_x[2] * adjacency + coef_x[3] * adjacency %*% adjacency
  cov_y <- coef_y[1] * diag(popn) + coef_y[2] * adjacency + coef_y[3] * adjacency %*% adjacency

  list(
    covariate_pre = pre_whiten(cov_x, x),
    outcome_pre = pre_whiten(cov_y, y),
    est_coef_x = coef_x,
    est_coef_y = coef_y
  )
}

estimate_mle_prewhiten <- function(adjacency, simulation, initial_values = NULL) {
  estimate_mle_prewhiten_vectors(
    adjacency = adjacency,
    x = simulation$covariate[[2]],
    y = simulation$outcome[[2]],
    initial_values = initial_values
  )
}
