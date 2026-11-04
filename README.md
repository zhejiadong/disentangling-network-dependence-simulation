# Disentangling network dependence among multiple variables

R code for the [paper](https://arxiv.org/abs/2506.20974). Install dependencies once:

```r
install.packages(c("igraph", "sna", "energy", "ape", "numDeriv"))
```

## Simulation

From the project root:

```r
setwd("simulation")
source("01_generate_data.R")
source("02_run_analysis.R")
source("03_make_table.R")
setwd("..")
```

Edit settings at the top of `01_generate_data.R`: 500 replications, 500 nodes, 500 edges, seed 321; equilibrium rho = 0.29/0.27/0.25 and direct kappa/alpha = (0.9,0.1)/(0.8,0.2)/(0.7,0.3), innovation SD 0.1. Generation saves `data/simulations.rds`, including the RNG checkpoint. Analysis restores that checkpoint and saves inspectable fits in `results/analysis.rds`. Each stage can run in a fresh R session from `simulation/`. The report reads those files without refitting and writes `results/table1.csv` in paper order.

## FHS application

Place authorized inputs in `fhs_application/data_raw/`, then from the project root:

```r
setwd("fhs_application")
source("01_prepare_data.R")
source("02_analysis.R")
source("03_tables.R")
setwd("..")
```

Use exactly one c1 and one c2 TAB-delimited table for each accession below, plain text or gzip. Filenames must contain the accession and delimited c1/c2 token; subfolders are supported. Confirm versions, units and coding against the authorized dbGaP dictionaries.

| Accession | Required columns |
|---|---|
| `pht000031` | `shareid`, `B13`, `B14`, `B24`, `B439`, `B443`, `B448` |
| `pht003099` | `shareid`, `age2`, `sex` |
| `pht000836` | `shareid`, `sharealterid`, `ALTERTYPE`, `SPELLBEGIN`, `SPELLEND` |

B13 is weight (kg), B14 height (cm), B24 first-physician SBP, and B439/B443/B448 septal thickness/posterior wall thickness/LV internal dimension (mm). The scripts construct BMI and LVM/height and analyze Social and Relatives networks separately, adjusting for age, sex and SBP. The intercept/covariate-adjusted MLE is in `R/prewhiten_adjusted.R`; `02_analysis.R` calls it directly and saves the fits. `03_tables.R` reads saved results and writes numeric CSVs.

FHS data are not included.

## Illustrations

Independently, from the project root, `setwd("illustrations"); source("01_figures.R"); setwd("..")` generates the synthetic figures. Its settings remain at the top of that script.
