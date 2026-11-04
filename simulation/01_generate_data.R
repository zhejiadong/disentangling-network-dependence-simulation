# Run from simulation/. Edit these settings before generating the data.
source("../R/simulation.R")
n_sim <- 500
popn <- 500
edges <- 500
seed <- 321
rho <- c(0.29, 0.27, 0.25)
k <- c(0.9, 0.8, 0.7)
alpha <- c(0.1, 0.2, 0.3)
epsilon <- 0.1

# One fixed network, followed by independent X and Y in each scenario.
network <- make_erdos_renyi_network(popn, edges, seed)
graph <- network$graph
adjacency <- network$adjacency
simulations <- simulate_scenarios(adjacency, n_sim, rho, k, alpha, epsilon, seed)

# Keep the RNG state so fitting in a new session uses the same starting values.
rng_after_generation <- .Random.seed
simulation_data <- list(network = network, simulations = simulations,
                        rng_after_generation = rng_after_generation)
dir.create("data", showWarnings = FALSE)
saveRDS(simulation_data, "data/simulations.rds")

# For example, inspect simulations[[1]][[1]] for the first dataset.
cat("Saved", n_sim, "replications for each of", length(simulations), "scenarios.\n")
