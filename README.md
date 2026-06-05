# Disentangling Network Dependence among Multiple Variables

This repository contains code for the paper [Disentangling network dependence among multiple variables](https://arxiv.org/abs/2506.20974).

## Repository Contents

- `R/simulation.R`: functions for network generation and data-generating process for the simulation studies.
- `R/prewhiten.R`: functions for the pre-whitening approach, constructing the proposed variance-covariance matrix, solving the mixed effect model through the optimization algorithm mentioned in the paper.
- `R/lnam.R`: functions to run the (linear) network autocorrelation model, fitted with `sna::lnam()`, and to extract the estimated covariance structure for pre-whitening.
- `R/summary.R`: functions for summarizing rejection rates across simulation replicates.
- `R/figure2.R`: code for Figure 2, illustrating how the variance-covariance structure changes under different transmission processes.
- `R/workflow.R`: a wrapper for running the simulation comparison from an R session.

## Dependencies

The main R packages are listed in `DESCRIPTION`. The analysis uses `igraph`, `sna`, `energy`, `MASS`, `expm`, and `optimx`.

## Simulated Data Structure

The function `simulate_scenarios()` returns a list of simulation scenarios. Each scenario contains a list of simulation replicates. Each replicate has the following structure:

```r
list(
  covariate = list(
    time0 = numeric(n),
    time1 = numeric(n)
  ),
  outcome = list(
    time0 = numeric(n),
    time1 = numeric(n)
  )
)
```

For LNAM-generated variables, the second element is named `timeequl` rather than `time1`, following the notation in the original simulation code. In both cases, the methods use the second element of each variable list as the observed network-dependent covariate and outcome.

The adjacency matrix is stored separately. It should be an `n` by `n` matrix with the same node order as the covariate and outcome vectors.

## Quick Check in R

To check that the code runs, start R from the repository root and run a small simulation with five replicates:

```r
source("R/simulation.R")
source("R/prewhiten.R")
source("R/lnam.R")
source("R/summary.R")
source("R/workflow.R")

small_table <- run_simulation_comparison(
  n_sim = 5,
  popn = 50,
  edges = 80,
  seed = 321,
  cores = 1
)

small_table
```

The output is a rejection-rate table. With only five replicates, the numbers are intended only as a code check, not as simulation evidence.

## Running the Proposed Pre-Whitening Method

With simulated data from this repository:

```r
source("R/simulation.R")
source("R/prewhiten.R")

network <- make_erdos_renyi_network(popn = 500, edges = 500, seed = 321)
simulations <- simulate_scenarios(network$adjacency, n_sim = 10)

fit <- estimate_mle_prewhiten(
  adjacency = network$adjacency,
  simulation = simulations[[1]][[1]]
)
```

For real data, use the vector-level function directly:

```r
source("R/prewhiten.R")

fit <- estimate_mle_prewhiten_vectors(
  adjacency = adjacency,
  x = exposure_vector,
  y = outcome_vector
)
```

Here `adjacency` is an `n` by `n` network adjacency matrix, and `exposure_vector` and `outcome_vector` are length-`n` numeric vectors ordered in the same way as the rows and columns of `adjacency`.

## Running the NAM method

The LNAM code uses the `lnam()` function from the `sna` package:

```r
source("R/prewhiten.R")
source("R/lnam.R")

fit <- estimate_lnam_prewhiten(
  adjacency = network$adjacency,
  simulation = simulations[[1]][[1]]
)
```

For real data:

```r
source("R/prewhiten.R")
source("R/lnam.R")

fit <- estimate_lnam_prewhiten_vectors(
  adjacency = adjacency,
  x = exposure_vector,
  y = outcome_vector
)
```

This fits `sna::lnam(y, x = x, W1 = adjacency)` and uses the estimated LNAM covariance structure for pre-whitening.

## Generating Figure 2

Figure 2 illustrates how the variance-covariance matrix changes under different underlying transmission processes. To generate the figure from R:

```r
source("R/simulation.R")
source("R/figure2.R")

make_figure2(
  output = "figure2_variance_covariance_structure.pdf",
  n_sim = 1000,
  seed = 124
)
```

The defaults follow the small synthetic network example in the original figure code: 10 nodes, 15 edges, and 1,000 simulation replicates. The output path can be changed through the `output` argument.
