simulate_figure2_covariances <- function(
    popn = 10,
    edges = 15,
    n_sim = 1000,
    seed = 124,
    k = 0.9,
    alpha = 0.1,
    epsilon = 0.5,
    long_time = 10000) {
  set.seed(seed)
  graph <- igraph::sample_gnm(popn, edges)
  adjacency <- as.matrix(igraph::as_adjacency_matrix(graph))

  direct_t1 <- matrix(0, nrow = n_sim, ncol = popn)
  direct_t2 <- matrix(0, nrow = n_sim, ncol = popn)
  direct_long <- matrix(0, nrow = n_sim, ncol = popn)

  for (sim_index in seq_len(n_sim)) {
    direct_t1[sim_index, ] <- simulate_process(
      adjacency,
      k = k,
      alpha = alpha,
      epsilon = epsilon,
      max_time = 1
    )$time1

    direct_t2[sim_index, ] <- simulate_process(
      adjacency,
      k = k,
      alpha = alpha,
      epsilon = epsilon,
      max_time = 2
    )$time2

    direct_long[sim_index, ] <- simulate_direct_endpoint(
      adjacency,
      k = 0.234,
      alpha = 0.09,
      epsilon = 0.1,
      max_time = long_time
    )
  }

  list(
    graph = graph,
    adjacency = adjacency,
    geodesic = igraph::distances(graph),
    covariance_t1 = stats::cov(direct_t1),
    covariance_t2 = stats::cov(direct_t2),
    covariance_long = stats::cov(direct_long)
  )
}

simulate_direct_endpoint <- function(adjacency, k, alpha, epsilon, max_time) {
  popn <- nrow(adjacency)
  outcome <- stats::rnorm(popn)

  for (time_index in seq_len(max_time)) {
    outcome <- k * adjacency %*% outcome + alpha * outcome + stats::rnorm(popn, 0, epsilon)
  }

  as.numeric(outcome)
}

draw_matrix_panel <- function(matrix_value, title, palette = grDevices::heat.colors(20)) {
  matrix_value <- as.matrix(matrix_value)
  if (any(is.infinite(matrix_value))) {
    finite_max <- max(matrix_value[is.finite(matrix_value)])
    matrix_value[is.infinite(matrix_value)] <- finite_max + 1
  }

  graphics::image(
    t(matrix_value[nrow(matrix_value):1, ]),
    col = palette,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title
  )
  graphics::box()
}

make_figure2 <- function(
    output = "figure2_variance_covariance_structure.pdf",
    popn = 10,
    edges = 15,
    n_sim = 1000,
    seed = 124,
    long_time = 10000) {
  figure_data <- simulate_figure2_covariances(
    popn = popn,
    edges = edges,
    n_sim = n_sim,
    seed = seed,
    long_time = long_time
  )

  output_dir <- dirname(output)
  if (!identical(output_dir, ".") && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  grDevices::pdf(output, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))

  plot(
    figure_data$graph,
    vertex.size = 12,
    vertex.label.cex = 0.8,
    edge.width = 1.5,
    layout = igraph::layout_with_fr(figure_data$graph),
    main = "(a) Network"
  )

  draw_matrix_panel(figure_data$adjacency, "(b) Adjacency matrix", palette = c("white", "steelblue"))
  draw_matrix_panel(figure_data$geodesic, "(c) Geodesic distance", palette = rev(grDevices::heat.colors(20)))
  draw_matrix_panel(figure_data$covariance_t1, expression("(d) " * bold(widehat(V)(t == 1))))
  draw_matrix_panel(figure_data$covariance_t2, expression("(e) " * bold(widehat(V)(t == 2))))
  draw_matrix_panel(figure_data$covariance_long, expression("(f) " * bold(widehat(V)(t == infinity))))

  invisible(output)
}
