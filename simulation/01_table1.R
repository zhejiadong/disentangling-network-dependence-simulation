# Table 1: run from simulation/. Full paper settings; edit only for a smoke test.
n_sim <- 500L
popn <- 500L
edges <- 500L
seed <- 321L
source("../R/simulation.R")
source("../R/prewhiten.R")
source("../R/lnam.R")
source("../R/summary.R")
source("../fhs_application/reporting.R")

network <- make_erdos_renyi_network(popn, edges, seed)
simulations <- simulate_scenarios(network$adjacency, n_sim=n_sim, seed=seed)
mle_adjustments <- lnam_adjustments <- vector("list",length(simulations))
for (scenario_index in seq_along(simulations)) {
  message("Fitting scenario ",names(simulations)[scenario_index])
  mle_adjustments[[scenario_index]] <- lapply(simulations[[scenario_index]],function(simulation) {
    estimate_mle_prewhiten(network$adjacency,simulation)
  })
  lnam_adjustments[[scenario_index]] <- lapply(simulations[[scenario_index]],function(simulation) {
    estimate_lnam_prewhiten(network$adjacency,simulation)
  })
}
rates <- summarize_rejection_rates(simulations,mle_adjustments,lnam_adjustments)
# The accepted source generates strongest first, equilibrium first. Paper order is reversed.
table1 <- rates[c(6,5,4,3,2,1),]
rownames(table1) <- NULL
stopifnot(all(is.finite(as.matrix(table1[,-1]))))
write_table(table1,"table1")
saveRDS(list(settings=list(n_sim=n_sim,n=popn,edges=edges,seed=seed),
  network=network, simulations=simulations, proposed=mle_adjustments, nam=lnam_adjustments),
  "results/simulation_private.rds")
writeLines(capture.output(sessionInfo()),"results/session_info.txt")
print(table1,row.names=FALSE)
