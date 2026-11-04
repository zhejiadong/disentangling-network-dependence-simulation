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

# Save numeric results; table objects remain available in the R session.
write_table <- function(table, stem) {
  dir.create("results", showWarnings=FALSE)
  write.csv(table, file.path("results", paste0(stem,".csv")), row.names=FALSE)
}
