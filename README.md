# Disentangling network dependence among multiple variables

R code for the JRSSA article ([paper](https://arxiv.org/abs/2506.20974)). The scripts generate Figures 1–2 and Tables 1–2, C1–C2. FHS data are not distributed.

## Run

Install the R dependencies once (`Matrix` is supplied with R):

```r
install.packages(c("igraph", "sna", "energy", "ape", "numDeriv"))
```

Run from the repository root:

```sh
(cd illustrations && Rscript 01_figures.R)
(cd simulation && Rscript 01_table1.R)
# Authorized FHS inputs must first be placed in fhs_application/data_raw/.
(cd fhs_application && Rscript 01_prepare_data.R)
(cd fhs_application && Rscript 02_analysis.R)
(cd fhs_application && Rscript 03_tables.R)
```

In RStudio, set the working directory to the script's folder and `source()` it. Settings appear at the top of each driver. The illustration uses 1,000 replicates and 10,000 updates; Table 1 uses 500 replicates per scenario, 500 nodes and 500 edges. Reduced settings are useful for code checks, not paper results.

## Article-to-code map

All filenames below are relative to the corresponding script folder.

| Item | Input | Code | Generated output |
|---|---|---|---|
| Figure 1 | Synthetic 10-node, 15-edge graph; seed 124 | `illustrations/01_figures.R` | `figures/figure1.pdf` (network, adjacency, geodesic, sample t=1, t=2, t=10,000) |
| Figure 2 | Same graph and t=1 process | `illustrations/01_figures.R` | `figures/figure2.pdf` (true vs sample t=1) |
| Table 1 | Independent X/Y on a synthetic 500-node, 500-edge graph | `simulation/01_table1.R`, `R/simulation.R`, `R/prewhiten.R`, `R/lnam.R`, `R/summary.R` | `results/table1.csv`, `results/table1.tex` |
| Table 2 | Authorized prepared FHS cohorts/graphs | `fhs_application/02_analysis.R`, `methods.R`, `03_tables.R` | `results/table2.csv`, `results/table2.tex` |
| Table C1 | Same two graphs | `fhs_application/03_tables.R`, `reporting.R` | `results/tableC1.csv`, `results/tableC1.tex` |
| Table C2 | Same aligned cohorts/graphs | `fhs_application/03_tables.R` | `results/tableC2.csv`, `results/tableC2.tex` |

CSVs preserve numeric results; LaTeX files use the article's table layouts. Figure scripts also save the individual panel PDFs. Figure 2 reuses Figure 1's sample covariance at t=1.

## Settings and reproduction

- Illustrations use a 10-node, 15-edge graph with seed 124. At t=1/2, kappa=0.9, alpha=0.1 and innovation SD=0.5; at t=10,000, kappa=0.234, alpha=0.09 and innovation SD=0.1. The theoretical t=1 covariance is `0.81 A^2 + 0.18 A + 0.26 I`.
- Table 1 uses direct kappa/alpha=(0.7,0.3), (0.8,0.2), (0.9,0.1), and equilibrium rho=0.25, 0.27, 0.29, with innovation SD=0.1. Proposed and NAM whitening fit separate X and Y covariance models. Incomplete/nonfinite p-values stop the summary rather than changing its denominator.
- FHS models use BMI as exposure and LVM/height as outcome, adjusted for age, sex and SBP. Proposed d=2 inference uses the joint observed Hessian and a 1.96 normal interval; OLS/NAM use residual-df t intervals. Table 2 displays residual Moran p-values, not coefficient p-values.

These are fresh seeded Monte Carlo runs. Historical intermediate RNG states and saved samples are unavailable, so exact printed-number replay is not promised; graph-library versions and base-R rendering may also differ. Session information is saved with outputs.

See [FHS inputs and construction](fhs_application/README.md) for accessions, fields and table definitions. Run protected FHS data only in the authorized environment. Raw/prepared data, IDs, graphs, fits, results and logs must remain private and are ignored by Git. Artificial-fixture checks are not a protected-data run: empirical reproduction still requires comparison of reconstructed IDs, edges and ordering with the frozen analysis inputs and regeneration of graph summaries.
