# Disentangling network dependence among multiple variables

R code for the [paper](https://arxiv.org/abs/2506.20974). Install dependencies once:

```r
install.packages(c("igraph", "sna", "energy", "ape", "numDeriv"))
```

## Simulation

Run from the repository root:

```r
setwd("simulation")
source("01_generate_data.R")
source("02_run_analysis.R")
source("03_make_table.R")
setwd("..")
```

Results are saved in `simulation/results/`; the summary table is `simulation/results/table1.csv`.

## FHS application

Run from the repository root:

```r
setwd("fhs_application")
source("01_prepare_data.R")
source("02_analysis.R")
source("03_tables.R")
setwd("..")
```

FHS data are not included.

## Illustrations

Run from the repository root:

```r
setwd("illustrations")
source("01_figures.R")
setwd("..")
```

Figures are saved in `illustrations/figures/`.
