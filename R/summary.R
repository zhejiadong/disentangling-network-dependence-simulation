summarize_association <- function(covariate, outcome) {
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

  list(
    cor_coef = cor_coef,
    distance_cor = distance_cor,
    distance_p = distance_p,
    lm_p = lm_p,
    lm_beta = lm_beta,
    z_beta = z_beta,
    lm_rej = mean(lm_p < 0.05, na.rm = TRUE),
    distance_rej = mean(distance_p < 0.05, na.rm = TRUE)
  )
}

matrix_from_simulations <- function(simulations, scenario_index, variable, time_index = 2) {
  do.call(
    rbind,
    lapply(simulations[[scenario_index]], function(simulation) {
      simulation[[variable]][[time_index]]
    })
  )
}

matrix_from_adjustments <- function(adjustments, scenario_index, variable) {
  do.call(
    rbind,
    lapply(adjustments[[scenario_index]], function(adjustment) {
      adjustment[[variable]]
    })
  )
}

summarize_rejection_rates <- function(simulations, mle_adjustments, lnam_adjustments) {
  scenario_count <- length(simulations)
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

    output[scenario_index, "Linear model, raw"] <- raw_summary$lm_rej
    output[scenario_index, "Distance correlation, raw"] <- raw_summary$distance_rej
    output[scenario_index, "LNAM coefficient test"] <- mean(stats::na.omit(lnam_p_values) < 0.05)
    output[scenario_index, "Linear model, pre-whitened"] <- mle_summary$lm_rej
    output[scenario_index, "Distance correlation, pre-whitened"] <- mle_summary$distance_rej
    output[scenario_index, "Linear model, LNAM pre-whitened"] <- lnam_summary$lm_rej
    output[scenario_index, "Distance correlation, LNAM pre-whitened"] <- lnam_summary$distance_rej
  }

  output
}
