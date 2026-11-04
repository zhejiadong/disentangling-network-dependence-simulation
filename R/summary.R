summarize_association <- function(covariate, outcome) {
  if (!is.matrix(covariate) || !is.matrix(outcome) ||
      !identical(dim(covariate), dim(outcome)) || nrow(covariate) == 0L) {
    stop("Association matrices must have matching, nonempty replication rows.")
  }
  if (any(!is.finite(covariate)) || any(!is.finite(outcome))) {
    stop("Association matrices contain missing or nonfinite values.")
  }
  n_sim <- nrow(covariate)
  cor_coef <- numeric(n_sim)
  distance_cor <- numeric(n_sim)
  distance_p <- numeric(n_sim)
  lm_p <- numeric(n_sim)
  lm_beta <- numeric(n_sim)
  z_beta <- numeric(n_sim)

  for (sim_index in seq_len(n_sim)) {
    cor_coef[sim_index] <- stats::cor(covariate[sim_index, ], outcome[sim_index, ])
    distance_cor[sim_index] <- energy::dcor(covariate[sim_index, ], outcome[sim_index, ])
    distance_test <- energy::dcorT.test(covariate[sim_index, ], outcome[sim_index, ])
    distance_p[sim_index] <- distance_test$p.value
    lm_fit <- stats::lm(outcome[sim_index, ] ~ covariate[sim_index, ])
    lm_summary <- summary(lm_fit)
    lm_beta[sim_index] <- lm_summary$coefficients[2, "Estimate"]
    lm_p[sim_index] <- lm_summary$coefficients[2, 4]
    z_beta[sim_index] <-
      lm_beta[sim_index] /
      (lm_summary$coefficients[2, "Std. Error"] * stats::sd(covariate[sim_index, ]))
  }

  if (length(lm_p) != n_sim || length(distance_p) != n_sim ||
      any(!is.finite(lm_p)) || any(!is.finite(distance_p))) {
    stop("Association p-values must be finite and complete for every replication.")
  }
  list(
    cor_coef = cor_coef,
    distance_cor = distance_cor,
    distance_p = distance_p,
    lm_p = lm_p,
    lm_beta = lm_beta,
    z_beta = z_beta,
    lm_rej = mean(lm_p < 0.05),
    distance_rej = mean(distance_p < 0.05)
  )
}

matrix_from_simulations <- function(simulations, scenario_index, variable, time_index = 2) {
  matrix <- do.call(
    rbind,
    lapply(simulations[[scenario_index]], function(simulation) {
      simulation[[variable]][[time_index]]
    })
  )
  if (is.null(nrow(matrix)) || nrow(matrix) != length(simulations[[scenario_index]])) {
    stop("Missing replication rows in matrix_from_simulations.")
  }
  matrix
}

matrix_from_adjustments <- function(adjustments, scenario_index, variable) {
  matrix <- do.call(
    rbind,
    lapply(adjustments[[scenario_index]], function(adjustment) {
      adjustment[[variable]]
    })
  )
  if (is.null(nrow(matrix)) || nrow(matrix) != length(adjustments[[scenario_index]])) {
    stop("Missing replication rows in matrix_from_adjustments.")
  }
  matrix
}

summarize_rejection_rates <- function(simulations, mle_adjustments, lnam_adjustments) {
  scenario_count <- length(simulations)
  if (length(mle_adjustments) != scenario_count || length(lnam_adjustments) != scenario_count) {
    stop("Adjustment scenario counts must match simulations.")
  }
  for (i in seq_len(scenario_count)) {
    n <- length(simulations[[i]])
    if (n == 0L || length(mle_adjustments[[i]]) != n || length(lnam_adjustments[[i]]) != n) {
      stop("Every summary must have the same nonzero replication rows in scenario ", i)
    }
  }
  output <- data.frame(
    Scenario = names(simulations),
    `Linear model, raw` = NA_real_,
    `Distance correlation, raw` = NA_real_,
    `LNAM coefficient test` = NA_real_,
    `Linear model, pre-whitened` = NA_real_,
    `Distance correlation, pre-whitened` = NA_real_,
    `Linear model, LNAM pre-whitened` = NA_real_,
    `Distance correlation, LNAM pre-whitened` = NA_real_,
    check.names = FALSE
  )

  for (scenario_index in seq_len(scenario_count)) {
    raw_summary <- summarize_association(
      covariate = matrix_from_simulations(simulations, scenario_index, "covariate", 2),
      outcome = matrix_from_simulations(simulations, scenario_index, "outcome", 2)
    )

    mle_summary <- summarize_association(
      covariate = matrix_from_adjustments(mle_adjustments, scenario_index, "covariate_pre"),
      outcome = matrix_from_adjustments(mle_adjustments, scenario_index, "outcome_pre")
    )

    lnam_summary <- summarize_association(
      covariate = matrix_from_adjustments(lnam_adjustments, scenario_index, "covariate_pre"),
      outcome = matrix_from_adjustments(lnam_adjustments, scenario_index, "outcome_pre")
    )

    lnam_p_values <- unlist(lapply(lnam_adjustments[[scenario_index]], function(adjustment) {
      adjustment$lnam_p
    }))

    if (length(lnam_p_values) != length(simulations[[scenario_index]]) ||
        any(!is.finite(lnam_p_values))) {
      stop("LNAM p-values must be finite and complete for every replication in scenario ", scenario_index)
    }
    output[scenario_index, "Linear model, raw"] <- raw_summary$lm_rej
    output[scenario_index, "Distance correlation, raw"] <- raw_summary$distance_rej
    output[scenario_index, "LNAM coefficient test"] <- mean(lnam_p_values < 0.05)
    output[scenario_index, "Linear model, pre-whitened"] <- mle_summary$lm_rej
    output[scenario_index, "Distance correlation, pre-whitened"] <- mle_summary$distance_rej
    output[scenario_index, "Linear model, LNAM pre-whitened"] <- lnam_summary$lm_rej
    output[scenario_index, "Distance correlation, LNAM pre-whitened"] <- lnam_summary$distance_rej
  }

  output
}
