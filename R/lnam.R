estimate_lnam <- function(adjacency, x, y) {
  sna::lnam(y, x = x, W1 = adjacency)
}

estimate_lnam_prewhiten_vectors <- function(adjacency, x, y) {
  popn <- nrow(adjacency)
  identity_matrix <- diag(popn)

  fit <- estimate_lnam(adjacency, x, y)
  rho_est <- stats::coef(fit)[2]
  beta_est <- stats::coef(fit)[1]
  t_value <- fit$beta / fit$beta.se
  lnam_p <- 2 * stats::pt(abs(t_value), fit$df.residual, lower.tail = FALSE)

  fit_x <- sna::lnam(x, W1 = adjacency)
  fit_y <- sna::lnam(y, W1 = adjacency)
  rho_est_x <- as.numeric(fit_x$rho1)
  rho_est_y <- as.numeric(fit_y$rho1)
  estimated_variance_x <-
    solve(identity_matrix - rho_est_x * adjacency) %*%
    solve(identity_matrix - rho_est_x * adjacency) *
    as.numeric(fit_x$sigmasq)
  estimated_variance_y <-
    solve(identity_matrix - rho_est_y * adjacency) %*%
    solve(identity_matrix - rho_est_y * adjacency) *
    as.numeric(fit_y$sigmasq)

  list(
    covariate_pre = pre_whiten(estimated_variance_x, x),
    outcome_pre = pre_whiten(estimated_variance_y, y),
    beta = beta_est,
    rho = rho_est,
    lnam_p = lnam_p
  )
}

estimate_lnam_prewhiten <- function(adjacency, simulation) {
  estimate_lnam_prewhiten_vectors(
    adjacency = adjacency,
    x = simulation$covariate[[2]],
    y = simulation$outcome[[2]]
  )
}
