# Figure 1 (six panels) and Figure 2 (true vs sample t=1).
# Run from illustrations/. All networks here are synthetic.
popn <- 10L
edges <- 15L
graph_seed <- 124L
simulation_seed <- 124L
n_sim <- 1000L
long_time <- 10000L
k <- 0.9
alpha <- 0.1
noise_sd <- 0.5
long_k <- 0.234
long_alpha <- 0.09
long_noise_sd <- 0.1

set.seed(graph_seed)
graph <- igraph::sample_gnm(popn, edges, directed=FALSE, loops=FALSE)
A <- as.matrix(igraph::as_adjacency_matrix(graph,sparse=FALSE))
Geodesic <- igraph::distances(graph)
B <- k*A + alpha*diag(popn)
True_t1 <- B %*% t(B) + noise_sd^2*diag(popn)
# Explicit new Monte Carlo stream: historical intermediate RNG states are unavailable.
# Rows are independent replicates. Updating all rows together preserves the process law.
set.seed(simulation_seed)
Y0_t1 <- matrix(rnorm(n_sim*popn),nrow=n_sim)
Y_t1 <- Y0_t1 %*% t(B) + matrix(rnorm(n_sim*popn,sd=noise_sd),nrow=n_sim)
Y_t2 <- matrix(rnorm(n_sim*popn),nrow=n_sim)
for (time in 1:2) {
  Y_t2 <- Y_t2 %*% t(B) + matrix(rnorm(n_sim*popn,sd=noise_sd),nrow=n_sim)
}
B_long <- long_k*A + long_alpha*diag(popn)
stopifnot(max(abs(eigen(B_long,symmetric=TRUE,only.values=TRUE)$values)) < 1)
Y_tinfty <- matrix(rnorm(n_sim*popn),nrow=n_sim)
for (time in seq_len(long_time)) {
  Y_tinfty <- Y_tinfty %*% t(B_long) +
    matrix(rnorm(n_sim*popn,sd=long_noise_sd),nrow=n_sim)
}
Sample_t1 <- cov(Y_t1)
Sample_t2 <- cov(Y_t2)
Sample_tinfty <- cov(Y_tinfty)

# Shared numeric scale for true/sample t=1; other covariance panels use their own scale.
palette <- grDevices::hcl.colors(40,"YlOrRd",rev=TRUE)
t1_limits <- range(True_t1,Sample_t1)
draw_matrix <- function(M,title,limits=range(M[is.finite(M)])) {
  graphics::par(mar=c(3,3,3,1))
  finite <- M[is.finite(M)]
  if (diff(limits)==0) limits <- limits+c(-.5,.5)
  display <- M
  display[!is.finite(display)] <- NA_real_
  image(seq_len(ncol(M)),seq_len(nrow(M)),t(display[nrow(M):1,,drop=FALSE]),
    col=palette,zlim=limits,axes=FALSE,xlab="Node",ylab="Node",asp=1,main=title)
  axis(1,at=seq_len(ncol(M)),cex.axis=.7)
  axis(2,at=seq_len(nrow(M)),labels=nrow(M):1,las=1,cex.axis=.7)
  # Values make every scale explicit; nonconnected geodesics are shown as Inf.
  for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M))) {
    text(j,nrow(M)+1-i,if(is.finite(M[i,j])) format(round(M[i,j],2),trim=TRUE) else "Inf",cex=.48)
  }
  box()
}
set.seed(graph_seed)
network_layout <- igraph::layout_with_fr(graph)
draw_network <- function() {
  par(mar=c(3,3,3,1))
  plot(graph,layout=network_layout,vertex.size=17,vertex.label.cex=.8,
    edge.color="grey45",main="(a) Network")
}
dir.create("figures",showWarnings=FALSE)
pdf("figures/figure1.pdf",width=12,height=8,useDingbats=FALSE)
par(mfrow=c(2,3))
draw_network()
draw_matrix(A,"(b) Adjacency matrix")
draw_matrix(Geodesic,"(c) Geodesic distance")
draw_matrix(Sample_t1,"(d) Sample t = 1",t1_limits)
draw_matrix(Sample_t2,"(e) Sample t = 2")
draw_matrix(Sample_tinfty,"(f) Sample t = 10,000")
dev.off()
pdf("figures/figure2.pdf",width=8,height=4,useDingbats=FALSE)
par(mfrow=c(1,2))
draw_matrix(True_t1,"(a) True t = 1",t1_limits)
draw_matrix(Sample_t1,"(b) Sample t = 1",t1_limits)
dev.off()
# Manuscript asset names; all generated from the same objects as the composites.
pdf("figures/network.pdf",width=5,height=5,useDingbats=FALSE)
draw_network(); dev.off()
panels <- list(Adj=A,Geodesic=Geodesic,Sample_t1=Sample_t1,Sample_t2=Sample_t2,
  Sample_tinfty=Sample_tinfty,True_t1=True_t1)
for (name in names(panels)) {
  pdf(file.path("figures",paste0(name,".pdf")),width=5,height=5,useDingbats=FALSE)
  draw_matrix(panels[[name]],gsub("_"," ",name),
    if(name %in% c("True_t1","Sample_t1")) t1_limits else range(panels[[name]][is.finite(panels[[name]])]))
  dev.off()
}
dir.create("results",showWarnings=FALSE)
saveRDS(list(A=A,True_t1=True_t1,Sample_t1=Sample_t1,Sample_t2=Sample_t2,
  Sample_tinfty=Sample_tinfty),"results/illustration_matrices.rds")
writeLines(capture.output(sessionInfo()),"results/session_info.txt")
