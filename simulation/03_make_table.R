# Run from simulation/; reuse saved data and fits without refitting.
source("../R/summary.R")
simulation_data <- readRDS("data/simulations.rds")
analysis_results <- readRDS("results/analysis.rds")
simulations <- simulation_data$simulations
mle_adjustments <- analysis_results$mle_adjustments
lnam_adjustments <- analysis_results$lnam_adjustments

rejection_rates <- summarize_rejection_rates(simulations, mle_adjustments, lnam_adjustments)
# Generation is strongest/equilibrium first; Table 1 uses the reverse order.
table1 <- rejection_rates[c(6, 5, 4, 3, 2, 1), ]
rownames(table1) <- NULL
stopifnot(all(is.finite(as.matrix(table1[, -1]))))
write.csv(table1, "results/table1.csv", row.names = FALSE)
print(table1, row.names = FALSE)
