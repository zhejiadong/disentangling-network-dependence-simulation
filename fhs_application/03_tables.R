# Work from fhs_application/ after 02_analysis.R. No real-network plots.
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
all_results <- read.csv("results/combined_results_all_specs.csv",stringsAsFactors=FALSE)
residual_moran <- read.csv("results/residual_moran.csv",stringsAsFactors=FALSE)
network_order <- c("relatives","social")
method_order <- c("ordinary","joint_lnam","joint_finite_d2")
key <- paste(all_results$network,all_results$method)
indices <- match(as.vector(t(outer(network_order,method_order,paste))),key)
stopifnot(!anyNA(indices), !anyDuplicated(key))
table2 <- all_results[indices,c("network","method","estimate","std_error","ci_lower_95","ci_upper_95","p_value","n")]
rkey <- paste(residual_moran$network,residual_moran$method)
table2$residual_moran_I <- residual_moran$moran_I[match(key[indices],rkey)]
table2$residual_moran_p_greater <- residual_moran$p_value_greater[match(key[indices],rkey)]
stopifnot(!anyNA(table2))

graphs <- cohorts <- list()
for (label in network_order) {
  graphs[[label]] <- readRDS(file.path("data",label,"network_private.rds"))$graph
  dat <- readRDS(file.path("data",label,"cohort_private.rds"))
  stopifnot(!anyDuplicated(dat$id),setequal(as.character(dat$id),igraph::V(graphs[[label]])$name))
  cohorts[[label]] <- dat[match(igraph::V(graphs[[label]])$name,as.character(dat$id)),]
  cohorts[[label]]$female <- as.integer(cohorts[[label]]$sex==2L)
}
tableC1 <- do.call(rbind,lapply(network_order,function(label) graph_summary(graphs[[label]],label)))
rows <- list()
for (label in network_order) {
  A <- as.matrix(igraph::as_adjacency_matrix(graphs[[label]],sparse=FALSE))
  for (variable in c("age","female","bmi","lvm_height1","sbp")) {
    x <- cohorts[[label]][[variable]]
    moran <- ape::Moran.I(x,A,scaled=FALSE,na.rm=FALSE,alternative="greater")
    rows[[paste(label,variable)]] <- data.frame(network=label,variable=variable,
      n=length(x),mean=mean(x),sd=sd(x),
      female_count=if(variable=="female") sum(x) else NA_real_,
      female_percent=if(variable=="female") 100*mean(x) else NA_real_,
      moran_I=moran$observed,analytic_p_greater=moran$p.value)
  }
}
tableC2 <- do.call(rbind,rows)
write.csv(table2, "results/table2.csv", row.names = FALSE)
write.csv(tableC1, "results/tableC1.csv", row.names = FALSE)
write.csv(tableC2, "results/tableC2.csv", row.names = FALSE)
print(table2,row.names=FALSE)
print(tableC1,row.names=FALSE)
print(tableC2,row.names=FALSE)
