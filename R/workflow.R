run_parallel <- function(items, worker, cores = 1) {
  if (cores <= 1) {
    lapply(items, worker)
  } else {
    parallel::mclapply(items, worker, mc.cores = cores)
  }
}

run_simulation_comparison <- function(
    n_sim = 5,
    popn = 50,
    edges = 80,
    seed = 321,
    cores = 1) {
  network <- make_erdos_renyi_network(
    popn = popn,
    edges = edges,
    seed = seed
  )

  simulations <- simulate_scenarios(
    adjacency = network$adjacency,
    n_sim = n_sim,
    seed = seed
  )

  scenario_indices <- seq_along(simulations)

  mle_adjustments <- run_parallel(scenario_indices, function(scenario_index) {
    lapply(simulations[[scenario_index]], function(simulation) {
      estimate_mle_prewhiten(network$adjacency, simulation)
    })
  }, cores = cores)
  names(mle_adjustments) <- names(simulations)

  lnam_adjustments <- run_parallel(scenario_indices, function(scenario_index) {
    lapply(simulations[[scenario_index]], function(simulation) {
      estimate_lnam_prewhiten(network$adjacency, simulation)
    })
  }, cores = cores)
  names(lnam_adjustments) <- names(simulations)

  summarize_rejection_rates(
    simulations = simulations,
    mle_adjustments = mle_adjustments,
    lnam_adjustments = lnam_adjustments
  )
}
