# Table summaries use only computed fits/cohorts/graphs.
graph_summary <- function(graph, label) {
  stopifnot(!igraph::is_directed(graph), igraph::is_simple(graph))
  degree <- igraph::degree(graph)
  comp <- igraph::components(graph)
  distance <- igraph::distances(graph)
  # Unordered distinct pairs, within components only: no diagonal or Inf.
  distance <- distance[upper.tri(distance) & is.finite(distance)]
  data.frame(network=label, n=igraph::vcount(graph), edges=igraph::ecount(graph),
    degree_mean=mean(degree), degree_median=median(degree),
    degree_min=min(degree), degree_max=max(degree), components=comp$no,
    largest_component=max(comp$csize), largest_percent=100*max(comp$csize)/length(degree),
    distance_mean=if(length(distance)) mean(distance) else NA_real_,
    distance_median=if(length(distance)) median(distance) else NA_real_,
    distance_min=if(length(distance)) min(distance) else NA_real_,
    distance_max=if(length(distance)) max(distance) else NA_real_)
}

# Formatting is confined to LaTeX; CSVs retain the computed numeric columns.
write_table <- function(table, stem) {
  dir.create("results", showWarnings=FALSE)
  write.csv(table, file.path("results", paste0(stem,".csv")), row.names=FALSE)
  f <- function(x, digits=3) formatC(x, format="f", digits=digits)
  p <- function(x) ifelse(x < .001, "$<0.001$", f(x))
  row <- function(...) paste0(paste(..., collapse=" & "), " \\\\")
  networks <- c("relatives", "social")
  if (stem == "table2") {
    methods <- c("ordinary", "joint_lnam", "joint_finite_d2")
    index <- match(as.vector(t(outer(networks, methods, paste))), paste(table$network,table$method))
    stopifnot(!anyNA(index))
    t <- table[index,]
    data <- rep(c("$(\\mathbf{X},\\mathbf{Y})$", "$(\\mathbf{X},\\mathbf{Y})$",
                  "$(\\mathbf{X}^{*},\\mathbf{Y}^{*})$"),2)
    association <- rep(c("$\\beta_{\\mathrm{OLS}}$", "$\\beta_{\\mathrm{NAM}}$", "$\\beta_{\\mathrm{OLS}}$"),2)
    body <- vapply(seq_len(6), function(i) row(c(
      c("Relatives", "", "", "Social", "", "")[i], data[i], association[i],
      paste0(f(t$estimate[i])," (",f(t$std_error[i]),")"),
      paste0("[",f(t$ci_lower_95[i]),", ",f(t$ci_upper_95[i]),"]"),
      p(t$residual_moran_p_greater[i]))), character(1))
    header <- row(c("Network", "Data", "Association", "Estimate (SE)", "95\\% CI", "Residual Moran $p$"))
    spec <- "lllrlr"
  } else if (stem == "tableC1") {
    t <- table[match(networks,table$network),]
    values <- lapply(seq_len(2), function(i) c(
      as.character(t$n[i]), as.character(t$edges[i]), f(t$degree_mean[i],2),
      paste0(t$degree_median[i]," [",t$degree_min[i],", ",t$degree_max[i],"]"),
      as.character(t$components[i]), paste0(t$largest_component[i]," (",f(t$largest_percent[i],1),")"),
      paste0(f(t$distance_mean[i],2)," (",t$distance_median[i],")"),
      paste0(t$distance_min[i],"--",t$distance_max[i])))
    labels <- c("Participants (nodes)", "Undirected ties (edges)", "Degree, mean",
      "Degree, median [range]", "Connected components", "Largest component, $n$ (\\%)",
      "Distance among connected pairs, mean (median)", "Distance among connected pairs, range")
    body <- vapply(seq_along(labels),function(i) row(c(labels[i],values[[1]][i],values[[2]][i])),character(1))
    header <- row(c("Network statistic", "Relatives", "Social")); spec <- "lrr"
  } else if (stem == "tableC2") {
    variables <- c("age", "female", "bmi", "lvm_height1", "sbp")
    labels <- c("Age, years", "Female sex, $n$ (\\%)", "BMI, kg/m$^2$", "LVM/height, g/m", "Systolic blood pressure, mmHg")
    body <- vapply(seq_along(variables),function(i) {
      cells <- unlist(lapply(networks,function(network) {
        t <- table[table$network==network & table$variable==variables[i],]
        stopifnot(nrow(t)==1L)
        summary <- if (variables[i]=="female") paste0(t$female_count," (",f(t$female_percent,1),")") else
          paste0(f(t$mean,2)," (",f(t$sd,2),")")
        c(summary,p(t$analytic_p_greater))
      }))
      row(c(labels[i],cells))
    },character(1))
    counts <- vapply(networks,function(n) unique(table$n[table$network==n]),numeric(1))
    header <- c(row(c("",paste0("\\multicolumn{2}{c}{Relatives ($n=",counts[1],"$)}"),
      paste0("\\multicolumn{2}{c}{Social ($n=",counts[2],"$)}"))),
      row(c("Characteristic", "Summary", "$p$-value", "Summary", "$p$-value")))
    spec <- "lrrrr"
  } else if (stem == "table1") {
    scenario <- as.character(table$Scenario)
    direct <- grepl("^direct_k_",scenario); equilibrium <- grepl("^lnam_rho_",scenario)
    stopifnot(nrow(table)==6L, ncol(table)==8L, sum(direct)==3L, sum(equilibrium)==3L)
    parameter <- as.numeric(sub("_alpha_.*", "", sub("^(direct_k_|lnam_rho_)", "",scenario)))
    index <- c(which(direct)[order(parameter[direct])], which(equilibrium)[order(parameter[equilibrium])])
    t <- table[index,]
    body <- character()
    for (block in 0:1) {
      body <- c(body, row(paste0("\\multicolumn{9}{l}{", c("(a) Direct transmission process at $t=1$", "(b) Equilibrium transmission process")[block+1], "}")), "\\hline")
      for (j in 1:3) body <- c(body,row(c(if(j==1) "Strength" else "",c("Weak","Medium","Strong")[j],f(as.numeric(t[block*3+j,-1])))))
    }
    header <- c(row(c("", "", "\\multicolumn{3}{c|}{$(\\mathbf{X},\\mathbf{Y})$}",
      "\\multicolumn{2}{c|}{$(\\mathbf{X}^{*},\\mathbf{Y}^{*})$}",
      "\\multicolumn{2}{c}{$(\\mathbf{X}^{*}_{\\mathrm{NAM}},\\mathbf{Y}^{*}_{\\mathrm{NAM}})$}")),
      row(c("", "", "$\\beta_{\\mathrm{OLS}}$", "$\\mathrm{dCorr}$", "$\\beta_{\\mathrm{NAM}}$",
        "$\\beta_{\\mathrm{OLS}}$", "$\\mathrm{dCorr}$", "$\\beta_{\\mathrm{OLS}}$", "$\\mathrm{dCorr}$")))
    spec <- "ll|ccc|cc|cc"
  } else stop("Unknown article table: ",stem)
  writeLines(c(paste0("\\begin{tabular}{",spec,"}"),header,"\\hline",body,"\\end{tabular}"),
    file.path("results",paste0(stem,".tex")))
}
