estimate_lnam <- function(adjacency, x, y) {
  sna::lnam(y, x = x, W1 = adjacency)
}

estimate_lnam_prewhiten_vectors <- function(adjacency, x, y) {
  popn <- nrow(adjacency)
  identity_matrix <- diag(popn)

  fit <- estimate_lnam(adjacency, x, y)
  rho_est <- stats::coef(fit)[2]
  beta_est <- stats::coef(fit)[1]
  sigsq_est <- as.numeric(fit$sigmasq)
  t_value <- fit$beta / fit$beta.se
  lnam_p <- 2 * stats::pt(abs(t_value), fit$df.residual, lower.tail = FALSE)

  estimated_variance <-
    solve(identity_matrix - rho_est * adjacency) %*%
    solve(identity_matrix - rho_est * adjacency) *
    sigsq_est

  list(
    covariate_pre = pre_whiten(estimated_variance, x),
    outcome_pre = pre_whiten(estimated_variance, y),
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
