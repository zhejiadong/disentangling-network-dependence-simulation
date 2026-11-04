# Adjusted Gaussian finite-covariance MLE with full observed-Hessian inference.
# input: symmetric adjacency A, numeric outcome Y, full-rank mean design D,
# and beta_column (the coefficient index used for the scalar inference fields).
# Include an intercept column in D when required; whitening transforms it too.
# Fit lm.fit(D_whitened, Y_whitened) without adding another intercept.
# d is the covariance polynomial degree.
prepare_spectral <- function(input) {
  eigen_A <- eigen(input$A, symmetric = TRUE)

  list(
    lambda = eigen_A$values,
    U = eigen_A$vectors,
    Y_rot = as.numeric(crossprod(eigen_A$vectors, input$Y)),
    D_rot = crossprod(eigen_A$vectors, input$D)
  )
}

fit_joint_finite <- function(
  input,
  spectral,
  d,
  n_random_starts = 12L,
  optimizer_maxit = 3000L
) {
  n <- length(input$Y)
  p <- ncol(input$D)
  q <- d + 1L
  lambda <- spectral$lambda
  U <- spectral$U
  Y_rot <- spectral$Y_rot
  D_rot <- spectral$D_rot
  powers <- 0:d

  variance_from_theta <- function(theta) {
    out <- rep(0, n)
    for (m in powers) {
      out <- out + theta[m + 1L] * lambda^m
    }
    out
  }

  profile_fit <- function(log_theta) {
    theta <- exp(log_theta)
    variance_values <- variance_from_theta(theta)

    if (
      any(!is.finite(variance_values)) ||
      min(variance_values) <= 1e-8
    ) {
      return(list(value = 1e100))
    }

    weights <- 1 / variance_values
    information <- crossprod(
      D_rot,
      D_rot * weights
    )
    score <- crossprod(
      D_rot,
      Y_rot * weights
    )

    alpha <- tryCatch(
      solve(information, score),
      error = function(e) NULL
    )

    if (is.null(alpha) || any(!is.finite(alpha))) {
      return(list(value = 1e100))
    }

    residual_rot <- as.numeric(
      Y_rot - D_rot %*% alpha
    )

    value <- 0.5 * (
      n * log(2 * pi) +
      sum(log(variance_values)) +
      sum(residual_rot^2 / variance_values)
    )

    list(
      value = value,
      alpha = as.numeric(alpha),
      theta = theta,
      variance_values = variance_values,
      information = information,
      residual_rot = residual_rot
    )
  }

  objective <- function(log_theta) {
    profile_fit(log_theta)$value
  }

  ordinary <- stats::lm.fit(input$D, input$Y)
  residual_variance <- sum(ordinary$residuals^2) / n

  base_theta <- numeric(q)
  base_theta[1] <- max(residual_variance, 1e-4)

  for (m in seq_len(d)) {
    eigen_scale <- max(abs(lambda)^m)
    base_theta[m + 1L] <- max(
      residual_variance * 0.01 / max(eigen_scale, 1),
      1e-8
    )
  }

  base_start <- log(base_theta)

  set.seed(20260820 + d)
  random_starts <- if (n_random_starts > 0L) {
    sweep(
      matrix(
        stats::rnorm(
          n_random_starts * q,
          0,
          0.8
        ),
        ncol = q
      ),
      2L,
      base_start,
      "+"
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = q)
  }

  starts <- rbind(base_start, random_starts)

  fits <- lapply(
    seq_len(nrow(starts)),
    function(i) {
      stats::optim(
        par = starts[i, ],
        fn = objective,
        method = "BFGS",
        control = list(maxit = optimizer_maxit)
      )
    }
  )

  fit_rows <- lapply(
    seq_along(fits),
    function(i) {
      fit <- fits[[i]]
      details <- profile_fit(fit$par)
      theta <- if (is.null(details$theta)) {
        rep(NA_real_, q)
      } else {
        details$theta
      }

      row <- data.frame(
        start = i,
        objective = fit$value,
        convergence = fit$convergence,
        min_variance_eigenvalue = if (
          is.null(details$variance_values)
        ) {
          NA_real_
        } else {
          min(details$variance_values)
        },
        stringsAsFactors = FALSE
      )

      for (m in powers) {
        row[[paste0("sigma", m)]] <- theta[m + 1L]
      }
      row
    }
  )

  fit_table <- do.call(rbind, fit_rows)

  valid <- (
    fit_table$convergence == 0 &
    is.finite(fit_table$objective) &
    fit_table$objective < 1e99
  )

  if (!any(valid)) {
    stop(
      paste("All finite d =", d, "profile starts failed."),
      call. = FALSE
    )
  }

  candidates <- which(valid)
  best_index <- candidates[
    which.min(fit_table$objective[candidates])
  ]
  best_optim <- fits[[best_index]]
  best <- profile_fit(best_optim$par)

  full_nll <- function(parameters) {
    alpha <- parameters[seq_len(p)]
    log_theta <- parameters[p + seq_len(q)]
    theta <- exp(log_theta)
    variance_values <- variance_from_theta(theta)

    if (
      any(!is.finite(variance_values)) ||
      min(variance_values) <= 1e-8
    ) {
      return(1e100)
    }

    residual_rot <- as.numeric(
      Y_rot - D_rot %*% alpha
    )

    0.5 * (
      n * log(2 * pi) +
      sum(log(variance_values)) +
      sum(residual_rot^2 / variance_values)
    )
  }

  full_start <- c(best$alpha, best_optim$par)

  full_optim <- stats::optim(
    par = full_start,
    fn = full_nll,
    method = "BFGS",
    control = list(maxit = optimizer_maxit)
  )

  if (
    full_optim$convergence != 0 ||
    !is.finite(full_optim$value)
  ) {
    stop(
      paste("Finite d =", d, "full optimization failed."),
      call. = FALSE
    )
  }

  log_theta_hat <- full_optim$par[p + seq_len(q)]
  theta_hat <- exp(log_theta_hat)
  variance_values <- variance_from_theta(theta_hat)
  weights <- 1 / variance_values

  information <- crossprod(
    D_rot,
    D_rot * weights
  )
  score <- crossprod(
    D_rot,
    Y_rot * weights
  )
  alpha_hat <- as.numeric(solve(information, score))

  conditional_vcov <- solve(information)
  conditional_se <- sqrt(
    conditional_vcov[
      input$beta_column,
      input$beta_column
    ]
  )

  joint_parameters <- c(alpha_hat, log_theta_hat)

  hessian <- numDeriv::hessian(
    func = full_nll,
    x = joint_parameters
  )
  hessian <- (hessian + t(hessian)) / 2
  hessian_values <- eigen(
    hessian,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  full_vcov <- tryCatch(
    solve(hessian),
    error = function(e) NULL
  )

  full_se <- if (
    min(hessian_values) > 1e-8 &&
    !is.null(full_vcov) &&
    all(is.finite(full_vcov)) &&
    full_vcov[
      input$beta_column,
      input$beta_column
    ] > 0
  ) {
    sqrt(
      full_vcov[
        input$beta_column,
        input$beta_column
      ]
    )
  } else {
    NA_real_
  }

  selected_se <- if (is.finite(full_se)) {
    full_se
  } else {
    conditional_se
  }

  se_type <- if (is.finite(full_se)) {
    "full_observed_hessian"
  } else {
    "conditional_gls_fallback"
  }

  beta_hat <- alpha_hat[input$beta_column]
  statistic <- beta_hat / selected_se
  p_value <- 2 * stats::pnorm(
    abs(statistic),
    lower.tail = FALSE
  )

  residual_rot <- as.numeric(
    Y_rot - D_rot %*% alpha_hat
  )

  whitened_residual <- as.numeric(
    U %*% (
      residual_rot / sqrt(variance_values)
    )
  )

  whitened_coefficients <- as.numeric(
    stats::lm.fit(
      D_rot / sqrt(variance_values),
      Y_rot / sqrt(variance_values)
    )$coefficients
  )

  equivalence_error <- max(
    abs(whitened_coefficients - alpha_hat)
  )

  list(
    d = d,
    beta = beta_hat,
    se = selected_se,
    se_type = se_type,
    conditional_se = conditional_se,
    full_se = full_se,
    statistic = statistic,
    p_value = p_value,
    ci_lower = beta_hat - 1.96 * selected_se,
    ci_upper = beta_hat + 1.96 * selected_se,
    coefficients = setNames(alpha_hat, colnames(input$D)),
    theta = setNames(
      theta_hat,
      paste0("sigma", powers)
    ),
    objective = full_nll(joint_parameters),
    profile_convergence = best_optim$convergence,
    full_convergence = full_optim$convergence,
    min_variance_eigenvalue = min(variance_values),
    hessian_min_eigenvalue = min(hessian_values),
    equivalence_error = equivalence_error,
    # Symmetric whitening in the original observation coordinates.
    Y_whitened = as.numeric(U %*% (Y_rot / sqrt(variance_values))),
    D_whitened = U %*% (D_rot / sqrt(variance_values)),
    # Retain the full joint covariance, not a redesigned uncertainty estimate.
    full_vcov = full_vcov,
    conditional_vcov = conditional_vcov,
    coefficient_se = setNames(
      if (min(hessian_values) > 1e-8 && !is.null(full_vcov) &&
          all(is.finite(full_vcov))) {
        sqrt(diag(full_vcov)[seq_len(p)])
      } else {
        rep(NA_real_, p)
      }, colnames(input$D)
    ),
    whitened_residual = whitened_residual,
    fit_table = fit_table
  )
}
