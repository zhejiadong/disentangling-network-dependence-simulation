make_erdos_renyi_network <- function(popn = 500, edges = 500, seed = 321) {
  set.seed(seed)
  graph <- igraph::sample_gnm(popn, edges)
  adjacency <- as.matrix(igraph::as_adjacency_matrix(graph))

  list(
    graph = graph,
    adjacency = adjacency,
    popn = popn,
    edges = edges,
    seed = seed
  )
}

simulate_process <- function(
    adjacency,
    rho = 1,
    k = 0.3,
    alpha = 0.3,
    epsilon = 0.1,
    max_time = 1,
    lnam = FALSE) {
  popn <- nrow(adjacency)
  identity_matrix <- diag(popn)
  outcome <- matrix(0, ncol = popn, nrow = max_time + 1)
  outcome[1, ] <- rnorm(popn)

  if (lnam) {
    outcome[2, ] <- qr.solve(identity_matrix - rho * adjacency, outcome[1, ])
  } else {
    for (time_index in 2:(max_time + 1)) {
      outcome[time_index, ] <-
        k * adjacency %*% outcome[time_index - 1, ] +
        alpha * outcome[time_index - 1, ] +
        rnorm(popn, 0, epsilon)
    }
  }

  outcomes <- lapply(1:nrow(outcome), function(row_index) {
    as.numeric(outcome[row_index, ])
  })

  if (lnam) {
    names(outcomes) <- c("time0", "timeequl")
  } else {
    names(outcomes) <- paste0("time", 0:max_time)
  }

  outcomes
}

simulate_scenarios <- function(
    adjacency,
    n_sim = 500,
    rho = c(0.29, 0.27, 0.25),
    k = c(0.9, 0.8, 0.7),
    alpha = c(0.1, 0.2, 0.3),
    epsilon = 0.1,
    seed = 321) {
  set.seed(seed)
  scenario_count <- length(rho) + length(k)
  simulations <- vector("list", scenario_count)

  for (scenario_index in seq_len(scenario_count)) {
    simulations[[scenario_index]] <- vector("list", n_sim)
  }

  for (sim_index in seq_len(n_sim)) {
    for (scenario_index in seq_along(rho)) {
      simulations[[scenario_index]][[sim_index]] <- list(
        covariate = simulate_process(
          adjacency,
          rho = rho[scenario_index],
          epsilon = epsilon,
          lnam = TRUE
        ),
        outcome = simulate_process(
          adjacency,
          rho = rho[scenario_index],
          epsilon = epsilon,
          lnam = TRUE
        )
      )
    }

    for (scenario_index in seq_along(k)) {
      list_index <- scenario_index + length(rho)
      simulations[[list_index]][[sim_index]] <- list(
        covariate = simulate_process(
          adjacency,
          k = k[scenario_index],
          alpha = alpha[scenario_index],
          epsilon = epsilon,
          lnam = FALSE
        ),
        outcome = simulate_process(
          adjacency,
          k = k[scenario_index],
          alpha = alpha[scenario_index],
          epsilon = epsilon,
          lnam = FALSE
        )
      )
    }
  }

  names(simulations) <- c(
    paste0("lnam_rho_", rho),
    paste0("direct_k_", k, "_alpha_", alpha)
  )

  simulations
}
