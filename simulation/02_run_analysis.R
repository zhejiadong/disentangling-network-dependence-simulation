# Run from simulation/ after 01_generate_data.R.
source("../R/prewhiten.R")
source("../R/lnam.R")
simulation_data <- readRDS("data/simulations.rds")
network <- simulation_data$network
adjacency <- network$adjacency
simulations <- simulation_data$simulations
rng_after_generation <- simulation_data$rng_after_generation
.Random.seed <- rng_after_generation

# Fit the marginal finite-covariance models for X and Y.
# Each saved fit includes covariance coefficients, optimizer output and whitened data.
mle_adjustments <- setNames(vector("list", length(simulations)), names(simulations))
for (scenario in names(simulations)) {
  mle_adjustments[[scenario]] <- vector("list", length(simulations[[scenario]]))
  for (replication in seq_along(simulations[[scenario]])) {
    dataset <- simulations[[scenario]][[replication]]
    x <- dataset$covariate[[2]]
    y <- dataset$outcome[[2]]
    mle_adjustments[[scenario]][[replication]] <-
      estimate_mle_prewhiten_vectors(adjacency, x, y)
  }
  cat(scenario, "finite pre-whitening complete\n")
}

# Fit the original simulation LNAM comparator and retain every fitted object.
lnam_adjustments <- setNames(vector("list", length(simulations)), names(simulations))
for (scenario in names(simulations)) {
  lnam_adjustments[[scenario]] <- vector("list", length(simulations[[scenario]]))
  for (replication in seq_along(simulations[[scenario]])) {
    dataset <- simulations[[scenario]][[replication]]
    x <- dataset$covariate[[2]]
    y <- dataset$outcome[[2]]
    lnam_adjustments[[scenario]][[replication]] <-
      estimate_lnam_prewhiten_vectors(adjacency, x, y)
  }
  cat(scenario, "LNAM pre-whitening complete\n")
}

# Inspect one fitted dataset without retaining a dense covariance for every repeat.
# Change these two indices to examine another scenario or replication.
inspect_scenario <- 1
inspect_replication <- 1
inspect_data <- simulations[[inspect_scenario]][[inspect_replication]]
mle_fit <- mle_adjustments[[inspect_scenario]][[inspect_replication]]
lnam_fit <- lnam_adjustments[[inspect_scenario]][[inspect_replication]]
coef_x <- mle_fit$est_coef_x
coef_y <- mle_fit$est_coef_y
covariance_x <- coef_x[1] * diag(nrow(adjacency)) +
  coef_x[2] * adjacency + coef_x[3] * adjacency %*% adjacency
covariance_y <- coef_y[1] * diag(nrow(adjacency)) +
  coef_y[2] * adjacency + coef_y[3] * adjacency %*% adjacency
# NAM retains separate X/Y whitened vectors, not a joint covariance object.
x_lnam_whitened <- lnam_fit$covariate_pre
y_lnam_whitened <- lnam_fit$outcome_pre
x_whitened <- mle_fit$covariate_pre
y_whitened <- mle_fit$outcome_pre

analysis_results <- list(mle_adjustments = mle_adjustments,
                         lnam_adjustments = lnam_adjustments)
dir.create("results", showWarnings = FALSE)
saveRDS(analysis_results, "results/analysis.rds")
