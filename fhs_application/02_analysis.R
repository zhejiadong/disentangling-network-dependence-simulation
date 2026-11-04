# Continuous BMI -> continuous LVM/height; full mean designs on fixed cohorts.
n_random_starts <- 12L
optimizer_maxit <- 3000L
library(igraph)
source("methods.R", local = TRUE)

formulas <- list(
  no_height = lvm_height1 ~ bmi + age + female + sbp)
cohorts <- networks <- graphs <- adjacency <- outcomes <- list()
for (label in c("social", "relatives")) {
  networks[[label]] <- readRDS(file.path("data", label, "network_private.rds"))
  graphs[[label]] <- networks[[label]]$graph
  dat <- readRDS(file.path("data", label, "cohort_private.rds"))
  cohorts[[label]] <- dat[match(V(graphs[[label]])$name, as.character(dat$id)), ]
  stopifnot(!anyDuplicated(dat$id), setequal(as.character(dat$id), V(graphs[[label]])$name),
    is_simple(graphs[[label]]), !is_directed(graphs[[label]]), all(degree(graphs[[label]]) > 0),
    !anyNA(cohorts[[label]]))
  cohorts[[label]]$female <- as.integer(cohorts[[label]]$sex == 2L)
  adjacency[[label]] <- as.matrix(as_adjacency_matrix(graphs[[label]], sparse = FALSE, attr = NULL))
  storage.mode(adjacency[[label]]) <- "double"
  diag(adjacency[[label]]) <- 0
  outcomes[[label]] <- as.numeric(cohorts[[label]]$lvm_height1)
}
social_cohort <- cohorts$social
relatives_cohort <- cohorts$relatives
social_graph <- graphs$social
relatives_graph <- graphs$relatives
social_adjacency <- adjacency$social
relatives_adjacency <- adjacency$relatives

# Inspect every design, rotated design, fit and whitened variable by network/spec.
design_matrices <- inputs <- spectral_objects <- fits <- whitened <- list()
result_rows <- diagnostic_rows <- residual_moran_rows <- list()
dir.create("results", showWarnings = FALSE)
for (label in names(cohorts)) {
  design_matrices[[label]] <- inputs[[label]] <- spectral_objects[[label]] <-
    fits[[label]] <- whitened[[label]] <- list()
  for (spec in names(formulas)) {
    dat <- cohorts[[label]]
    D <- model.matrix(formulas[[spec]], data = dat)
    stopifnot(qr(D)$rank == ncol(D), all(is.finite(D)))
    design_matrices[[label]][[spec]] <- D
    input <- list(label = label, dat = dat, D = D, Y = outcomes[[label]],
      A = adjacency[[label]], beta_column = match("bmi", colnames(D)))
    inputs[[label]][[spec]] <- input
    spectral <- prepare_spectral(input)
    spectral_objects[[label]][[spec]] <- spectral

    ordinary <- lm(formulas[[spec]], data = dat)
    fits[[label]][[spec]] <- list(ordinary = ordinary)
    finite <- fit_joint_finite(input, spectral, d = 2L)
    fits[[label]][[spec]]$joint_finite_d2 <- finite
    whitened[[label]][[spec]] <- list(Y = finite$Y_whitened,
      D = finite$D_whitened, residual = finite$whitened_residual)
    finite_ok <- identical(finite$se_type, "full_observed_hessian") &&
      is.finite(finite$full_se) && finite$full_se > 0 &&
      is.finite(finite$min_variance_eigenvalue) && finite$min_variance_eigenvalue > 0 &&
      is.finite(finite$hessian_min_eigenvalue) && finite$hessian_min_eigenvalue > 1e-8 &&
      finite$profile_convergence == 0L && finite$full_convergence == 0L &&
      is.finite(finite$equivalence_error) && finite$equivalence_error < 1e-6
    diagnostic_rows[[paste(label, spec)]] <- data.frame(network = label, spec = spec,
      profile_convergence = finite$profile_convergence, full_convergence = finite$full_convergence,
      hessian_min_eigenvalue = finite$hessian_min_eigenvalue,
      min_variance_eigenvalue = finite$min_variance_eigenvalue,
      gls_equivalence_error = finite$equivalence_error, full_se = finite$full_se,
      finite_ok = finite_ok)
    print(diagnostic_rows[[paste(label, spec)]], row.names = FALSE)
    saveRDS(fits, "results/fits_private.rds")
    if (!finite_ok) stop("Finite fit needs review; inspect finite and fits. No fallback SE is reported.")

    lnam <- fit_joint_lnam(input)
    fits[[label]][[spec]]$joint_lnam <- lnam
    lnam_stable <- abs(lnam$rho) * max(abs(spectral$lambda)) < 1
    if (!setequal(names(lnam$coefficients), colnames(D)) ||
        !is.finite(lnam$se) || lnam$se <= 0 || !isTRUE(lnam_stable)) {
      stop("LNAM design, standard error or stability requires review; inspect lnam and fits.")
    }
    diagnostic_rows[[paste(label, spec)]]$lnam_rho <- lnam$rho
    diagnostic_rows[[paste(label, spec)]]$lnam_stable <- lnam_stable

    # OLS uses t intervals; finite GLS uses its original 1.96 normal interval;
    # sna LNAM uses its model-based SE and residual-df t interval.
    ordinary_table <- coef(summary(ordinary))
    ordinary_bmi <- ordinary_table["bmi", ]
    ordinary_critical <- qt(0.975, df = ordinary$df.residual)
    if (any(!is.finite(c(ordinary_bmi, finite$beta, finite$se, lnam$beta, lnam$se))) ||
        ordinary_bmi[2] <= 0) {
      stop("Non-finite coefficient/inference or nonpositive SE; inspect the fits before reporting.")
    }
    result_rows[[paste(label, spec)]] <- data.frame(
      spec = spec, network = label, method = c("ordinary", "joint_finite_d2", "joint_lnam"),
      estimate = c(ordinary_bmi[1], finite$beta, lnam$beta),
      std_error = c(ordinary_bmi[2], finite$se, lnam$se),
      statistic = c(ordinary_bmi[3], finite$statistic, lnam$statistic),
      p_value = c(ordinary_bmi[4], finite$p_value, lnam$p_value),
      ci_lower_95 = c(ordinary_bmi[1] - ordinary_critical * ordinary_bmi[2], finite$ci_lower, lnam$ci_lower),
      ci_upper_95 = c(ordinary_bmi[1] + ordinary_critical * ordinary_bmi[2], finite$ci_upper, lnam$ci_upper),
      se_type = c("ordinary_lm", finite$se_type, "sna_lnam_model_based"),
      scale = "continuous_raw_BMI_kg_m2_and_LVM_height1_g_per_m", n = nrow(dat), row.names = NULL)
    residuals_by_method <- list(ordinary = residuals(ordinary),
      joint_finite_d2 = finite$whitened_residual, joint_lnam = lnam$disturbances)
    for (method in names(residuals_by_method)) {
      moran <- ape::Moran.I(residuals_by_method[[method]], weight = input$A,
        scaled = FALSE, na.rm = FALSE, alternative = "greater")
      residual_moran_rows[[paste(label, spec, method)]] <- data.frame(
        network = label, spec = spec, method = method, moran_I = moran$observed,
        expected = moran$expected, sd = moran$sd, p_value_greater = moran$p.value)
    }
    saveRDS(fits, "results/fits_private.rds")
  }
}
all_results <- do.call(rbind, result_rows)
fit_diagnostics <- do.call(rbind, diagnostic_rows)
residual_moran <- do.call(rbind, residual_moran_rows)
write.csv(all_results, "results/combined_results_all_specs.csv", row.names = FALSE)
write.csv(fit_diagnostics, "results/fit_diagnostics.csv", row.names = FALSE)
write.csv(residual_moran, "results/residual_moran.csv", row.names = FALSE)

# Descriptives use each network's own nonisolated cohort, not a common network sample.
descriptive_rows <- raw_moran_rows <- list()
for (label in names(cohorts)) {
  dat <- cohorts[[label]]
  descriptive_rows[[label]] <- data.frame(network = label, n = nrow(dat),
    edges = ecount(graphs[[label]]), bmi_mean = mean(dat$bmi), bmi_sd = sd(dat$bmi),
    lvm_height1_mean = mean(dat$lvm_height1), lvm_height1_sd = sd(dat$lvm_height1),
    age_mean = mean(dat$age), age_sd = sd(dat$age),
    female_n = sum(dat$female), female_percent = 100 * mean(dat$female),
    sbp_mean = mean(dat$sbp), sbp_sd = sd(dat$sbp))
  for (variable in c("age", "female", "bmi", "lvm_height1", "sbp")) {
    moran <- ape::Moran.I(dat[[variable]], weight = adjacency[[label]],
      scaled = FALSE, na.rm = FALSE, alternative = "greater")
    raw_moran_rows[[paste(label, variable)]] <- data.frame(network = label,
      variable = variable, moran_I = moran$observed, expected = moran$expected,
      sd = moran$sd, standardized_I = (moran$observed - moran$expected) / moran$sd,
      analytic_p_greater = moran$p.value)
  }
}
descriptives <- do.call(rbind, descriptive_rows)
raw_moran <- do.call(rbind, raw_moran_rows)
write.csv(descriptives, "results/aggregate_descriptives.csv", row.names = FALSE)
write.csv(raw_moran, "results/raw_analytic_moran.csv", row.names = FALSE)
print(all_results, row.names = FALSE, digits = 5)
