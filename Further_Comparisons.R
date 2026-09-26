###############################################################################
# Script: Further_Comparisons.R
#
# Purpose: Simulation experiments for the article
# "A TWO-SAMPLE TEST BASED ON AVERAGED WILCOXON RANK SUMS OVER INTERPOINT DISTANCES",
# in particular, further power studies and benchmark comparisons for the averaged 
# interpoint Wilcoxon statistic (corresponding to Figures 10-14 of the supplement).
#
# Simulations include:
#   1) power comparisons on uniform distributions on spheres
#      (concentric spheres and shifted/intersecting spheres),
#   2) power and level comparisons in the uniform ellipse problem,
#   3) size-adjusted power comparisons in normal location and scale problems,
#   4) size-adjusted power comparisons for different choices of lp metrics,
#   5) sensitivity of power to the metric parameter p and to growing dimension.
#
# Analysis overview:
# - For each setting, samples X and Y are generated under the specified model.
# - Interpoint distance matrices are computed and the averaged interpoint
#   Wilcoxon statistic is evaluated using scenario-specific P1 values.
# - The performance of the averaged Wilcoxon test is compared to benchmark
#   procedures including depth-based, distance-based, kernel-based and
#   nearest-neighbor methods.
# - In the size-adjusted power experiments, null distributions are first
#   calibrated and empirical thresholds are estimated before power is computed
#   under alternatives.
# - Results are collected into data frames and plotted as multi-panel power
#   and level curves.
#
# Required packages:
# - ggplot2
# - dplyr
# - tidyr
# - patchwork
# - energy
# - SpatialNP
# - DepthProc
# - kernlab
# - Ecume
# - FNN
# - npmv
#
# Required source file:
# - Utilities.R
#   This script assumes that helper functions such as
#   `generate_quadruples()`, `generate_uniform_sphere()`,
#   `generate_sphere_samples()` and `generate_ellipse_samples()`
#   are available in the R session or sourced from Utilities.R
#   before execution.
#
# Outputs:
# - Data frames of estimated rejection rates and size-adjusted power
#   for each simulation setting
# - Multi-panel plots corresponding to Figures 10-14
# - Printed intermediate simulation summaries in the console
#
# Author: Aljosa Marjanovic
# Date: 10.03.2026
# Project: Two sample averaged Wilcoxon over interpoint distances
###############################################################################

source("Utilities.R")


################################################################################
# Fig. 10: Uniform distributions on spheres
# Left: concentric spheres (scale problem)
# Right: intersecting spheres (location problem)
################################################################################

set.seed(1)

MC <- 10
ns <- c(25, 50, 75, 100)
P1_unif <- 0.52

################################################################################
# wrappers for competitor tests
################################################################################

get_p_dd <- function(X, Y) {
  mWilcoxonTest(X, Y)$p.value
}

get_p_energy <- function(X, Y) {
  eqdist.etest(rbind(X, Y), sizes = c(nrow(X), nrow(Y)), R = 1000)$p.value
}

get_p_spatial <- function(X, Y) {
  sr.loc.test(X, Y, score = "rank", cond.n = 10)$p.value
}

get_p_knn <- function(X, Y) {
  knn_test(X, Y)$p_value
}

get_p_mmd1 <- function(X, Y) {
  mmd_test(X, Y, kernel = "rbfdot", sigma = 1)$p.value
}

get_p_mmd_opt <- function(X, Y, sigma_grid = NULL) {
  Z <- rbind(X, Y)
  D <- as.matrix(dist(Z))
  dvec <- D[upper.tri(D)]
  
  if (is.null(sigma_grid)) {
    dmed <- median(dvec)
    sigma_grid <- dmed * c(0.25, 0.5, 1, 2, 4)
  }
  
  pvals <- sapply(sigma_grid, function(s) {
    mmd_test(X, Y, kernel = "rbfdot", sigma = s)$p.value
  })
  
  min(pvals)
}

knn_test <- function(X, Y, k = 3, n_perm = 1000) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  data <- rbind(X, Y)
  labels <- c(rep(1, n1), rep(2, n2))
  
  knn_result <- get.knn(data, k = k)
  
  neighbor_labels <- matrix(labels[knn_result$nn.index], nrow = n1 + n2, ncol = k)
  same_class_counts <- rowSums(neighbor_labels == labels)
  test_statistic <- sum(same_class_counts)
  
  permuted_stats <- replicate(n_perm, {
    permuted_labels <- sample(labels)
    permuted_neighbor_labels <- matrix(
      permuted_labels[knn_result$nn.index],
      nrow = n1 + n2,
      ncol = k
    )
    sum(rowSums(permuted_neighbor_labels == permuted_labels))
  })
  
  p_value <- mean(permuted_stats >= test_statistic)
  
  list(statistic = test_statistic, p_value = p_value)
}

pval_marginal_ranks <- function(X, Y) {
  grp <- factor(c(rep("X", nrow(X)), rep("Y", nrow(Y))))
  Z   <- as.data.frame(rbind(X, Y))
  names(Z) <- paste0("V", seq_len(ncol(Z)))
  Z$grp <- grp
  
  form <- as.formula(
    paste(paste(names(Z)[-length(names(Z))], collapse = "|"), "~ grp")
  )
  
  fit <- npmv::nonpartest(
    form,
    Z,
    permtest   = FALSE,
    tests      = c(1, 0, 0, 0),
    plots      = FALSE,
    releffects = FALSE
  )
  
  as.numeric(fit[1, "P-value"])
}

################################################################################
# Left panel: concentric spheres
# X ~ U(S^2(0,1)), Y ~ U(S^2(0,1.1))
################################################################################

fig10_left <- data.frame(
  n = ns,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  indices <- generate_quadruples(n, m)
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (i in 1:MC) {
    X <- generate_uniform_sphere(n, 3, radius = 1)
    Y <- generate_uniform_sphere(m, 3, radius = 1.1)
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    Z <- rbind(X, Y)
    D <- as.matrix(dist(Z))
   
    
    beta_P1 <- beta_P1 + as.numeric(pnorm(abs(W), lower.tail = FALSE) < 0.025)
    beta_DD <- beta_DD + as.numeric(get_p_dd(X, Y) < 0.05)
    beta_Energy <- beta_Energy + as.numeric(get_p_energy(X, Y) < 0.05)
    beta_kNN <- beta_kNN + as.numeric(get_p_knn(X, Y) < 0.05)
    beta_MMD1 <- beta_MMD1 + as.numeric(get_p_mmd1(X, Y) < 0.05)
    beta_MMD <- beta_MMD + as.numeric(get_p_mmd_opt(X, Y) < 0.05)
    beta_Spatial <- beta_Spatial + as.numeric(get_p_spatial(X, Y) < 0.05)
    beta_MR <- beta_MR + as.numeric(pval_marginal_ranks(X, Y) < 0.05)
  }
  
  fig10_left$P1[k] <- beta_P1 / MC
  fig10_left$DD[k] <- beta_DD / MC
  fig10_left$Energy[k] <- beta_Energy / MC
  fig10_left$kNN[k] <- beta_kNN / MC
  fig10_left$MMD1[k] <- beta_MMD1 / MC
  fig10_left$MMD[k] <- beta_MMD / MC
  fig10_left$Spatial[k] <- beta_Spatial / MC
  fig10_left$MR[k] <- beta_MR / MC
  
  print("fig10_left")
  print(fig10_left[k, ])
  print("--------------------")
}

################################################################################
# Right panel: intersecting spheres
# X ~ U(S^2(c,1)), Y ~ U(S^2(-c,1)), c = (0.1, 0, 0)
################################################################################

fig10_right <- data.frame(
  n = ns,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  indices <- generate_quadruples(n, m)
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (i in 1:MC) {
    X <- generate_sphere_samples(n, 3, radius = 1, center = c(0.1, 0, 0))
    Y <- generate_sphere_samples(m, 3, radius = 1, center = c(-0.1, 0, 0))
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    Z <- rbind(X, Y)
    D <- as.matrix(dist(Z))
   
    
    beta_P1 <- beta_P1 + as.numeric(pnorm(abs(W), lower.tail = FALSE) < 0.025)
    beta_DD <- beta_DD + as.numeric(get_p_dd(X, Y) < 0.05)
    beta_Energy <- beta_Energy + as.numeric(get_p_energy(X, Y) < 0.05)
    beta_kNN <- beta_kNN + as.numeric(get_p_knn(X, Y) < 0.05)
    beta_MMD1 <- beta_MMD1 + as.numeric(get_p_mmd1(X, Y) < 0.05)
    beta_MMD <- beta_MMD + as.numeric(get_p_mmd_opt(X, Y) < 0.05)
    beta_Spatial <- beta_Spatial + as.numeric(get_p_spatial(X, Y) < 0.05)
    beta_MR <- beta_MR + as.numeric(pval_marginal_ranks(X, Y) < 0.05)
  }
  
  fig10_right$P1[k] <- beta_P1 / MC
  fig10_right$DD[k] <- beta_DD / MC
  fig10_right$Energy[k] <- beta_Energy / MC
  fig10_right$kNN[k] <- beta_kNN / MC
  fig10_right$MMD1[k] <- beta_MMD1 / MC
  fig10_right$MMD[k] <- beta_MMD / MC
  fig10_right$Spatial[k] <- beta_Spatial / MC
  fig10_right$MR[k] <- beta_MR / MC
  
  print("fig10_right")
  print(fig10_right[k, ])
  print("--------------------")
}

fig10_results <- list(
  concentric_spheres = fig10_left,
  intersecting_spheres = fig10_right
)

fig10_results

###############################################################################
# Figure 10
###############################################################################

fig10_left_plot  <- fig10_results$concentric_spheres
fig10_right_plot <- fig10_results$intersecting_spheres

fig10_legend_order <- c("DD", "kNN", "P1", "MMD1", "Energy", "MMD", "Spatial", "MR")

make_fig10_df <- function(df) {
  df |>
    dplyr::rename(
      P1 = P1,
      Spatial = Spatial,
      DD = DD,
      MMD1 = MMD1,
      MMD = MMD,
      kNN = kNN,
      Energy = Energy,
      MR = MR
    ) |>
    tidyr::pivot_longer(
      cols = -n,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig10_legend_order))
}

plot_left_df  <- make_fig10_df(fig10_left_plot)
plot_right_df <- make_fig10_df(fig10_right_plot)

fig10_cols <- c(
  DD = "blue",
  kNN = "#C77CFF",
  P1 = "red",
  MMD1 = "black",
  Energy = "#F4B6C2",
  MMD = "#00CD00",
  Spatial = "black",
  MR = "#8B4513"
)

fig10_types <- c(
  DD = "dashed",
  kNN = "solid",
  P1 = "solid",
  MMD1 = "solid",
  Energy = "solid",
  MMD = "dotdash",
  Spatial = "dashed",
  MR = "twodash"
)

make_fig10_panel <- function(df, panel_title) {
  ggplot(df, aes(x = n, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig10_cols, breaks = fig10_legend_order) +
    scale_linetype_manual(values = fig10_types, breaks = fig10_legend_order) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = c(25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Sample size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig10_panel(
  plot_left_df,
  "Power: Uniform on concentric spheres, small effect"
)

p2 <- make_fig10_panel(
  plot_right_df,
  "Power: Uniform on a sphere, shifted spheres (intersection)"
)

p1 | p2


################################################################################
# Uniform ellipses problem
# Left: power for circle vs ellipse
# Right: level for X, Y ~ U(S^2(0,1))
################################################################################

################################################################################
# Left panel: X ~ U(S^1(0,1)), Y ~ U(E(1.1, 0.6))
################################################################################

fig11_left <- data.frame(
  n = ns,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  indices <- generate_quadruples(n, m)
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (i in 1:MC) {
    X <- generate_uniform_sphere(n, d = 2, radius = 1)
    Y <- generate_ellipse_samples(m, a = 1.1, e = 0.6)
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    Z <- rbind(X, Y)
    D <- as.matrix(dist(Z))
    
    
    beta_P1 <- beta_P1 + as.numeric(pnorm(abs(W), lower.tail = FALSE) < 0.025)
    beta_DD <- beta_DD + as.numeric(get_p_dd(X, Y) < 0.05)
    beta_Energy <- beta_Energy + as.numeric(get_p_energy(X, Y) < 0.05)
    beta_kNN <- beta_kNN + as.numeric(get_p_knn(X, Y) < 0.05)
    beta_MMD1 <- beta_MMD1 + as.numeric(get_p_mmd1(X, Y) < 0.05)
    beta_MMD <- beta_MMD + as.numeric(get_p_mmd_opt(X, Y) < 0.05)
    beta_Spatial <- beta_Spatial + as.numeric(get_p_spatial(X, Y) < 0.05)
    beta_MR <- beta_MR + as.numeric(pval_marginal_ranks(X, Y) < 0.05)
  }
  
  fig11_left$P1[k] <- beta_P1 / MC
  fig11_left$DD[k] <- beta_DD / MC
  fig11_left$Energy[k] <- beta_Energy / MC
  fig11_left$kNN[k] <- beta_kNN / MC
  fig11_left$MMD1[k] <- beta_MMD1 / MC
  fig11_left$MMD[k] <- beta_MMD / MC
  fig11_left$Spatial[k] <- beta_Spatial / MC
  fig11_left$MR[k] <- beta_MR / MC
  
  print("fig11_left")
  print(fig11_left[k, ])
  print("--------------------")
}

################################################################################
# Right panel: X, Y ~ U(S^2(0,1))
################################################################################

fig11_right <- data.frame(
  n = ns,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  indices <- generate_quadruples(n, m)
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (i in 1:MC) {
    X <- generate_uniform_sphere(n, d = 3, radius = 1)
    Y <- generate_uniform_sphere(m, d = 3, radius = 1)
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    Z <- rbind(X, Y)
    D <- as.matrix(dist(Z))
    
    
    beta_P1 <- beta_P1 + as.numeric(pnorm(abs(W), lower.tail = FALSE) < 0.025)
    beta_DD <- beta_DD + as.numeric(get_p_dd(X, Y) < 0.05)
    beta_Energy <- beta_Energy + as.numeric(get_p_energy(X, Y) < 0.05)
    beta_kNN <- beta_kNN + as.numeric(get_p_knn(X, Y) < 0.05)
    beta_MMD1 <- beta_MMD1 + as.numeric(get_p_mmd1(X, Y) < 0.05)
    beta_MMD <- beta_MMD + as.numeric(get_p_mmd_opt(X, Y) < 0.05)
    beta_Spatial <- beta_Spatial + as.numeric(get_p_spatial(X, Y) < 0.05)
    beta_MR <- beta_MR + as.numeric(pval_marginal_ranks(X, Y) < 0.05)
  }
  
  fig11_right$P1[k] <- beta_P1 / MC
  fig11_right$DD[k] <- beta_DD / MC
  fig11_right$Energy[k] <- beta_Energy / MC
  fig11_right$kNN[k] <- beta_kNN / MC
  fig11_right$MMD1[k] <- beta_MMD1 / MC
  fig11_right$MMD[k] <- beta_MMD / MC
  fig11_right$Spatial[k] <- beta_Spatial / MC
  fig11_right$MR[k] <- beta_MR / MC
  
  print("fig11_right")
  print(fig11_right[k, ])
  print("--------------------")
}

fig11_results <- list(
  ellipse_problem = fig11_left,
  sphere_level = fig11_right
)

fig11_results


###############################################################################
# Figure 11
###############################################################################

fig11_left_plot  <- fig11_results$ellipse_problem
fig11_right_plot <- fig11_results$sphere_level

fig11_legend_order <- c("DD", "kNN", "P1", "MMD1", "Energy", "MMD", "Spatial", "MR")

make_fig11_df <- function(df) {
  df |>
    dplyr::rename(
      P1 = P1,
      Spatial = Spatial,
      DD = DD,
      MMD1 = MMD1,
      MMD = MMD,
      kNN = kNN,
      Energy = Energy,
      MR = MR
    ) |>
    tidyr::pivot_longer(
      cols = -n,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig11_legend_order))
}

plot_left_df  <- make_fig11_df(fig11_left_plot)
plot_right_df <- make_fig11_df(fig11_right_plot)

fig11_cols <- c(
  DD = "blue",
  kNN = "#C77CFF",
  P1 = "red",
  MMD1 = "black",
  Energy = "#F4B6C2",
  MMD = "#00CD00",
  Spatial = "black",
  MR = "#8B4513"
)

fig11_types <- c(
  DD = "dashed",
  kNN = "solid",
  P1 = "solid",
  MMD1 = "solid",
  Energy = "solid",
  MMD = "dotdash",
  Spatial = "dashed",
  MR = "twodash"
)

make_fig11_panel <- function(df, panel_title) {
  ggplot(df, aes(x = n, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig11_cols, breaks = fig11_legend_order) +
    scale_linetype_manual(values = fig11_types, breaks = fig11_legend_order) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = c(25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Sample size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig11_panel(
  plot_left_df,
  "Power: Uniform ellipses problem"
)

p2 <- make_fig11_panel(
  plot_right_df,
  "Level: Uniform on the sphere"
)

gridExtra::grid.arrange(p1, p2, nrow = 1)

###############################################################################
# Size-adjusted power, normal location and normal scale
###############################################################################

alpha <- 0.05
B <- 1250

Delta_grid <- seq(0, 1.4, by = 0.2)

n <- 25
m <- 25
d <- 3
P1_norm <- 0.535

fig12_left <- data.frame(
  Delta = Delta_grid,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)


indices <- generate_quadruples(n, m)

###############################################################################
# Null calibration
###############################################################################

p0_P1 <- numeric(B)
p0_DD <- numeric(B)
p0_Energy <- numeric(B)
p0_kNN <- numeric(B)
p0_MMD1 <- numeric(B)
p0_MMD <- numeric(B)
p0_Spatial <- numeric(B)
p0_MR <- numeric(B)

for (b in 1:B) {
  
  X <- matrix(rnorm(n*d), n, d)
  Y <- matrix(rnorm(m*d), m, d)
  
  dist_mat_XX <- as.matrix(dist(X))
  dist_mat_XY <- as.matrix(dist(rbind(X,Y)))
  
  T <- outer(dist_mat_XX, dist_mat_XY,
             FUN = function(x,y) as.numeric(x <= y))
  
  W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
    ((m+n-6)*(n-3)*(n-4)*length(T[indices])*(2*P1_norm-1)/8)^(1/2)
  
  p0_P1[b]      <- 2*pnorm(-abs(W))
  p0_DD[b]      <- get_p_dd(X,Y)
  p0_Energy[b]  <- get_p_energy(X,Y)
  p0_kNN[b]     <- get_p_knn(X,Y)
  p0_MMD1[b]    <- get_p_mmd1(X,Y)
  p0_MMD[b]     <- get_p_mmd_opt(X,Y)
  p0_Spatial[b] <- get_p_spatial(X,Y)
  p0_MR[b]      <- pval_marginal_ranks(X,Y)
  
}

###############################################################################
# Size-adjusted thresholds
###############################################################################

thr_P1      <- quantile(p0_P1, alpha, na.rm=TRUE)
thr_DD      <- quantile(p0_DD, alpha, na.rm=TRUE)
thr_Energy  <- quantile(p0_Energy, alpha, na.rm=TRUE)
thr_kNN     <- quantile(p0_kNN, alpha, na.rm=TRUE)
thr_MMD1    <- quantile(p0_MMD1, alpha, na.rm=TRUE)
thr_MMD     <- quantile(p0_MMD, alpha, na.rm=TRUE)
thr_Spatial <- quantile(p0_Spatial, alpha, na.rm=TRUE)
thr_MR      <- quantile(p0_MR, alpha, na.rm=TRUE)

###############################################################################
# Power under location alternatives
###############################################################################

for (k in seq_along(Delta_grid)) {
  
  Delta <- Delta_grid[k]
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n*d), n, d)
    Y <- matrix(rnorm(m*d), m, d)
    
    Y <- sweep(Y, 2, rep(Delta,d), "+")
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X,Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x,y) as.numeric(x <= y))
    
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m+n-6)*(n-3)*(n-4)*length(T[indices])*(2*P1_norm-1)/8)^(1/2)
    
    p_P1      <- 2*pnorm(-abs(W))
    p_DD      <- get_p_dd(X,Y)
    p_Energy  <- get_p_energy(X,Y)
    p_kNN     <- get_p_knn(X,Y)
    p_MMD1    <- get_p_mmd1(X,Y)
    p_MMD     <- get_p_mmd_opt(X,Y)
    p_Spatial <- get_p_spatial(X,Y)
    p_MR      <- pval_marginal_ranks(X,Y)
    
    beta_P1      <- beta_P1      + as.numeric(p_P1 < thr_P1)
    beta_DD      <- beta_DD      + as.numeric(p_DD < thr_DD)
    beta_Energy  <- beta_Energy  + as.numeric(p_Energy < thr_Energy)
    beta_kNN     <- beta_kNN     + as.numeric(p_kNN < thr_kNN)
    beta_MMD1    <- beta_MMD1    + as.numeric(p_MMD1 < thr_MMD1)
    beta_MMD     <- beta_MMD     + as.numeric(p_MMD < thr_MMD)
    beta_Spatial <- beta_Spatial + as.numeric(p_Spatial < thr_Spatial)
    beta_MR      <- beta_MR      + as.numeric(p_MR < thr_MR)
    
  }
  
  fig12_left$P1[k]      <- beta_P1/B
  fig12_left$DD[k]      <- beta_DD/B
  fig12_left$Energy[k]  <- beta_Energy/B
  fig12_left$kNN[k]     <- beta_kNN/B
  fig12_left$MMD1[k]    <- beta_MMD1/B
  fig12_left$MMD[k]     <- beta_MMD/B
  fig12_left$Spatial[k] <- beta_Spatial/B
  fig12_left$MR[k]      <- beta_MR/B
  
  print(fig12_left[k,])
}

fig12_results <- list(location_problem = fig12_left)


###############################################################################
# Size-adjusted power under scale alternatives
###############################################################################

s_grid <- seq(1, 0.1, by = -0.1)

fig12_right <- data.frame(
  s = s_grid,
  P1 = NA_real_,
  DD = NA_real_,
  Energy = NA_real_,
  kNN = NA_real_,
  MMD1 = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_,
  MR = NA_real_
)

for (k in seq_along(s_grid)) {
  
  s <- s_grid[k]
  
  beta_P1 <- 0
  beta_DD <- 0
  beta_Energy <- 0
  beta_kNN <- 0
  beta_MMD1 <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  beta_MR <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n*d), n, d)
    Y <- matrix(rnorm(m*d), m, d)
    
    Y <- Y * sqrt(s)
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X,Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x,y) as.numeric(x <= y))
    
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m+n-6)*(n-3)*(n-4)*length(T[indices])*(2*P1_norm-1)/8)^(1/2)
    
    p_P1      <- 2*pnorm(-abs(W))
    p_DD      <- get_p_dd(X,Y)
    p_Energy  <- get_p_energy(X,Y)
    p_kNN     <- get_p_knn(X,Y)
    p_MMD1    <- get_p_mmd1(X,Y)
    p_MMD     <- get_p_mmd_opt(X,Y)
    p_Spatial <- get_p_spatial(X,Y)
    p_MR      <- pval_marginal_ranks(X,Y)
    
    beta_P1      <- beta_P1      + as.numeric(p_P1 < thr_P1)
    beta_DD      <- beta_DD      + as.numeric(p_DD < thr_DD)
    beta_Energy  <- beta_Energy  + as.numeric(p_Energy < thr_Energy)
    beta_kNN     <- beta_kNN     + as.numeric(p_kNN < thr_kNN)
    beta_MMD1    <- beta_MMD1    + as.numeric(p_MMD1 < thr_MMD1)
    beta_MMD     <- beta_MMD     + as.numeric(p_MMD < thr_MMD)
    beta_Spatial <- beta_Spatial + as.numeric(p_Spatial < thr_Spatial)
    beta_MR      <- beta_MR      + as.numeric(p_MR < thr_MR)
    
  }
  
  fig12_right$P1[k]      <- beta_P1/B
  fig12_right$DD[k]      <- beta_DD/B
  fig12_right$Energy[k]  <- beta_Energy/B
  fig12_right$kNN[k]     <- beta_kNN/B
  fig12_right$MMD1[k]    <- beta_MMD1/B
  fig12_right$MMD[k]     <- beta_MMD/B
  fig12_right$Spatial[k] <- beta_Spatial/B
  fig12_right$MR[k]      <- beta_MR/B
  
  print(fig12_right[k,])
}

fig12_results$scale_problem <- fig12_right

fig12_results

###############################################################################
# Figure 12
###############################################################################

fig12_left_plot  <- fig12_results$location_problem
fig12_right_plot <- fig12_results$scale_problem

fig12_legend_order <- c("DD", "kNN", "P1", "MMD1", "Energy", "MMD", "Spatial", "MR")

make_fig12_df_left <- function(df) {
  df |>
    dplyr::rename(
      P1 = P1,
      Spatial = Spatial,
      DD = DD,
      MMD1 = MMD1,
      MMD = MMD,
      kNN = kNN,
      Energy = Energy,
      MR = MR
    ) |>
    tidyr::pivot_longer(
      cols = -Delta,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig12_legend_order))
}

make_fig12_df_right <- function(df) {
  df |>
    dplyr::rename(
      P1 = P1,
      Spatial = Spatial,
      DD = DD,
      MMD1 = MMD1,
      MMD = MMD,
      kNN = kNN,
      Energy = Energy,
      MR = MR
    ) |>
    tidyr::pivot_longer(
      cols = -s,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig12_legend_order))
}

plot_left_df  <- make_fig12_df_left(fig12_left_plot)
plot_right_df <- make_fig12_df_right(fig12_right_plot)

fig12_cols <- c(
  DD = "blue",
  kNN = "#C77CFF",
  P1 = "red",
  MMD1 = "black",
  Energy = "#F4B6C2",
  MMD = "#00CD00",
  Spatial = "black",
  MR = "#8B4513"
)

fig12_types <- c(
  DD = "dashed",
  kNN = "solid",
  P1 = "solid",
  MMD1 = "solid",
  Energy = "solid",
  MMD = "dotdash",
  Spatial = "dashed",
  MR = "twodash"
)

make_fig12_panel_left <- function(df, panel_title) {
  ggplot(df, aes(x = Delta, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig12_cols, breaks = fig12_legend_order) +
    scale_linetype_manual(values = fig12_types, breaks = fig12_legend_order) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = seq(0, 1.4, by = 0.2)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Effect size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

make_fig12_panel_right <- function(df, panel_title) {
  ggplot(df, aes(x = s, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig12_cols, breaks = fig12_legend_order) +
    scale_linetype_manual(values = fig12_types, breaks = fig12_legend_order) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
    ) +
    scale_x_reverse(breaks = seq(1, 0.1, by = -0.1)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Scale ratio",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig12_panel_left(
  plot_left_df,
  "Size-adjusted power: normal location problem"
)

p2 <- make_fig12_panel_right(
  plot_right_df,
  "Size-adjusted power: normal scale problem"
)

gridExtra::grid.arrange(p1, p2, nrow = 1)

###############################################################################
# Size-adjusted power for different l_p metrics
###############################################################################

Delta_grid_13 <- seq(0, 3, by = 0.25)

n <- 25
m <- 25
d <- 100

P1_l1   <- 0.538
P1_l2   <- 0.537
P1_l10  <- 0.53
P1_linf <- 0.52

fig13_left <- data.frame(
  Delta = Delta_grid_13,
  l1 = NA_real_,
  l2 = NA_real_,
  l10 = NA_real_,
  linf = NA_real_
)

indices <- generate_quadruples(n, m)

###############################################################################
# Null calibration
###############################################################################

p0_l1   <- numeric(B)
p0_l2   <- numeric(B)
p0_l10  <- numeric(B)
p0_linf <- numeric(B)

for (b in 1:B) {
  
  X <- matrix(rnorm(n * d), n, d)
  Y <- matrix(rnorm(m * d), m, d)
  
  dist_mat_XX_l1 <- as.matrix(dist(X, method = "manhattan"))
  dist_mat_XY_l1 <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
  
  T_l1 <- outer(dist_mat_XX_l1, dist_mat_XY_l1,
                FUN = function(x, y) as.numeric(x <= y))
  
  W_l1 <- (sum(T_l1[indices]) - 0.5 * length(T_l1[indices])) /
    ((m + n - 6) * (n - 3) * (n - 4) * length(T_l1[indices]) * (2 * P1_l1 - 1) / 8)^(1/2)
  
  
  dist_mat_XX_l2 <- as.matrix(dist(X, method = "euclidean"))
  dist_mat_XY_l2 <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
  
  T_l2 <- outer(dist_mat_XX_l2, dist_mat_XY_l2,
                FUN = function(x, y) as.numeric(x <= y))
  
  W_l2 <- (sum(T_l2[indices]) - 0.5 * length(T_l2[indices])) /
    ((m + n - 6) * (n - 3) * (n - 4) * length(T_l2[indices]) * (2 * P1_l2 - 1) / 8)^(1/2)
  
  
  dist_mat_XX_l10 <- as.matrix(dist(X, method = "minkowski", p = 10))
  dist_mat_XY_l10 <- as.matrix(dist(rbind(X, Y), method = "minkowski", p = 10))
  
  T_l10 <- outer(dist_mat_XX_l10, dist_mat_XY_l10,
                 FUN = function(x, y) as.numeric(x <= y))
  
  W_l10 <- (sum(T_l10[indices]) - 0.5 * length(T_l10[indices])) /
    ((m + n - 6) * (n - 3) * (n - 4) * length(T_l10[indices]) * (2 * P1_l10 - 1) / 8)^(1/2)
  
  
  dist_mat_XX_linf <- as.matrix(dist(X, method = "maximum"))
  dist_mat_XY_linf <- as.matrix(dist(rbind(X, Y), method = "maximum"))
  
  T_linf <- outer(dist_mat_XX_linf, dist_mat_XY_linf,
                  FUN = function(x, y) as.numeric(x <= y))
  
  W_linf <- (sum(T_linf[indices]) - 0.5 * length(T_linf[indices])) /
    ((m + n - 6) * (n - 3) * (n - 4) * length(T_linf[indices]) * (2 * P1_linf - 1) / 8)^(1/2)
  
  
  p0_l1[b]   <- 2 * pnorm(-abs(W_l1))
  p0_l2[b]   <- 2 * pnorm(-abs(W_l2))
  p0_l10[b]  <- 2 * pnorm(-abs(W_l10))
  p0_linf[b] <- 2 * pnorm(-abs(W_linf))
}

thr_l1   <- quantile(p0_l1, alpha, na.rm = TRUE)
thr_l2   <- quantile(p0_l2, alpha, na.rm = TRUE)
thr_l10  <- quantile(p0_l10, alpha, na.rm = TRUE)
thr_linf <- quantile(p0_linf, alpha, na.rm = TRUE)

###############################################################################
# Size-adjusted power under location alternatives
###############################################################################

for (k in seq_along(Delta_grid_13)) {
  
  Delta <- Delta_grid_13[k]
  
  beta_l1   <- 0
  beta_l2   <- 0
  beta_l10  <- 0
  beta_linf <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n * d), n, d)
    Y <- matrix(rnorm(m * d), m, d)
    
    Y[, 1] <- Y[, 1] + Delta
    
    dist_mat_XX_l1 <- as.matrix(dist(X, method = "manhattan"))
    dist_mat_XY_l1 <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
    
    T_l1 <- outer(dist_mat_XX_l1, dist_mat_XY_l1,
                  FUN = function(x, y) as.numeric(x <= y))
    
    W_l1 <- (sum(T_l1[indices]) - 0.5 * length(T_l1[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l1[indices]) * (2 * P1_l1 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_l2 <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY_l2 <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    T_l2 <- outer(dist_mat_XX_l2, dist_mat_XY_l2,
                  FUN = function(x, y) as.numeric(x <= y))
    
    W_l2 <- (sum(T_l2[indices]) - 0.5 * length(T_l2[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l2[indices]) * (2 * P1_l2 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_l10 <- as.matrix(dist(X, method = "minkowski", p = 10))
    dist_mat_XY_l10 <- as.matrix(dist(rbind(X, Y), method = "minkowski", p = 10))
    
    T_l10 <- outer(dist_mat_XX_l10, dist_mat_XY_l10,
                   FUN = function(x, y) as.numeric(x <= y))
    
    W_l10 <- (sum(T_l10[indices]) - 0.5 * length(T_l10[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l10[indices]) * (2 * P1_l10 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_linf <- as.matrix(dist(X, method = "maximum"))
    dist_mat_XY_linf <- as.matrix(dist(rbind(X, Y), method = "maximum"))
    
    T_linf <- outer(dist_mat_XX_linf, dist_mat_XY_linf,
                    FUN = function(x, y) as.numeric(x <= y))
    
    W_linf <- (sum(T_linf[indices]) - 0.5 * length(T_linf[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_linf[indices]) * (2 * P1_linf - 1) / 8)^(1/2)
    
    
    p_l1   <- 2 * pnorm(-abs(W_l1))
    p_l2   <- 2 * pnorm(-abs(W_l2))
    p_l10  <- 2 * pnorm(-abs(W_l10))
    p_linf <- 2 * pnorm(-abs(W_linf))
    
    beta_l1   <- beta_l1   + as.numeric(p_l1   < thr_l1)
    beta_l2   <- beta_l2   + as.numeric(p_l2   < thr_l2)
    beta_l10  <- beta_l10  + as.numeric(p_l10  < thr_l10)
    beta_linf <- beta_linf + as.numeric(p_linf < thr_linf)
  }
  
  fig13_left$l1[k]   <- beta_l1 / B
  fig13_left$l2[k]   <- beta_l2 / B
  fig13_left$l10[k]  <- beta_l10 / B
  fig13_left$linf[k] <- beta_linf / B
  
  print(fig13_left[k, ])
}

###############################################################################
# Figure 13 (right):
# Power degradation with growing dimension
###############################################################################

dims_13 <- c(10, 25, 50, 75, 100, 125, 150)

n <- 50
m <- 50
Delta <- 0.5
P1_norm <- 0.54

fig13_right <- data.frame(
  d = dims_13,
  DD = NA_real_,
  kNN = NA_real_,
  P1 = NA_real_,
  Energy = NA_real_,
  MMD = NA_real_,
  Spatial = NA_real_
)

for (k in seq_along(dims_13)) {
  
  d <- dims_13[k]
  indices <- generate_quadruples(n, m)
  
  beta_DD <- 0
  beta_kNN <- 0
  beta_P1 <- 0
  beta_Energy <- 0
  beta_MMD <- 0
  beta_Spatial <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n * d), n, d)
    Y <- matrix(rnorm(m * d), m, d)
    
    Y[, 1] <- Y[, 1] + Delta
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_norm - 1) / 8)^(1/2)
    
    p_P1 <- 2 * pnorm(-abs(W))
    
    p_DD <- tryCatch(get_p_dd(X, Y), error = function(e) NA_real_)
    p_kNN <- tryCatch(get_p_knn(X, Y), error = function(e) NA_real_)
    p_Energy <- tryCatch(get_p_energy(X, Y), error = function(e) NA_real_)
    p_MMD <- tryCatch(get_p_mmd1(X, Y), error = function(e) NA_real_)
    p_Spatial <- tryCatch(get_p_spatial(X, Y), error = function(e) NA_real_)
    
    if (is.null(p_DD) || length(p_DD) == 0 || !is.finite(p_DD)) p_DD <- 1
    if (is.null(p_kNN) || length(p_kNN) == 0 || !is.finite(p_kNN)) p_kNN <- 1
    if (is.null(p_Energy) || length(p_Energy) == 0 || !is.finite(p_Energy)) p_Energy <- 1
    if (is.null(p_MMD) || length(p_MMD) == 0 || !is.finite(p_MMD)) p_MMD <- 1
    if (is.null(p_Spatial) || length(p_Spatial) == 0 || !is.finite(p_Spatial)) p_Spatial <- 1
    
    beta_DD <- beta_DD + as.numeric(p_DD < alpha)
    beta_kNN <- beta_kNN + as.numeric(p_kNN < alpha)
    beta_P1 <- beta_P1 + as.numeric(p_P1 < alpha)
    beta_Energy <- beta_Energy + as.numeric(p_Energy < alpha)
    beta_MMD <- beta_MMD + as.numeric(p_MMD < alpha)
    beta_Spatial <- beta_Spatial + as.numeric(p_Spatial < alpha)
  }
  
  fig13_right$DD[k] <- beta_DD / B
  fig13_right$kNN[k] <- beta_kNN / B
  fig13_right$P1[k] <- beta_P1 / B
  fig13_right$Energy[k] <- beta_Energy / B
  fig13_right$MMD[k] <- beta_MMD / B
  fig13_right$Spatial[k] <- beta_Spatial / B
  
  print(fig13_right[k, ])
}

fig13_results$dimensionality <- fig13_right

fig13_results <- list(
  lp_size_adjusted = fig13_left,
  dimensionality = fig13_right
)

###############################################################################
# Figure 13
###############################################################################

fig13_left_plot  <- fig13_results$lp_size_adjusted
fig13_right_plot <- fig13_results$dimensionality

fig13_left_legend_order  <- c("l1", "l2", "l10", "max")
fig13_right_legend_order <- c("DD", "kNN", "P1", "Energy", "MMD", "Spatial")

make_fig13_df_left <- function(df) {
  df |>
    tidyr::pivot_longer(
      cols = -Delta,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig13_left_legend_order))
}

make_fig13_df_right <- function(df) {
  df |>
    tidyr::pivot_longer(
      cols = -d,
      names_to = "Method",
      values_to = "Power"
    ) |>
    dplyr::mutate(Method = factor(Method, levels = fig13_right_legend_order))
}

plot_left_df  <- make_fig13_df_left(fig13_left_plot)
plot_right_df <- make_fig13_df_right(fig13_right_plot)

fig13_left_cols <- c(
  l1  = "red",
  l2  = "blue",
  l10 = "black",
  max = "#00CD00"
)

fig13_left_types <- c(
  l1  = "dashed",
  l2  = "solid",
  l10 = "dotdash",
  max = "twodash"
)

fig13_right_cols <- c(
  DD      = "blue",
  kNN     = "#C77CFF",
  P1      = "red",
  Energy  = "#F4B6C2",
  MMD     = "#00CD00",
  Spatial = "black"
)

fig13_right_types <- c(
  DD      = "dashed",
  kNN     = "solid",
  P1      = "solid",
  Energy  = "solid",
  MMD     = "dotdash",
  Spatial = "dashed"
)

make_fig13_panel_left <- function(df, panel_title) {
  ggplot(df, aes(x = Delta, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig13_left_cols, breaks = fig13_left_legend_order) +
    scale_linetype_manual(values = fig13_left_types, breaks = fig13_left_legend_order) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 1, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Effect size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

make_fig13_panel_right <- function(df, panel_title) {
  ggplot(df, aes(x = d, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig13_right_cols, breaks = fig13_right_legend_order) +
    scale_linetype_manual(values = fig13_right_types, breaks = fig13_right_legend_order) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = c(10, 25, 50, 75, 100, 125, 150)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = panel_title,
      x = "Dimension",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig13_panel_left(
  plot_left_df,
  "Size-adjusted power: different l_p metrics"
)

p2 <- make_fig13_panel_right(
  plot_right_df,
  "Power degradation with growing dimension"
)

gridExtra::grid.arrange(p1, p2, nrow = 1)

###############################################################################
###############################################################################
# Power curves for l_p induced metrics p = 1, 2, 10, inf, and Canberra
# dim = 100, sparse normal location problem, varying sample size
###############################################################################
###############################################################################
ns_14 <- c(25, 50, 75, 100)

d <- 100
Delta <- 0.5


fig14_left <- data.frame(
  n = ns_14,
  l1 = NA_real_,
  l2 = NA_real_,
  l10 = NA_real_,
  max = NA_real_,
  Canberra = NA_real_
)

for (k in seq_along(ns_14)) {
  
  n <- ns_14[k]
  m <- n
  indices <- generate_quadruples(n, m)
  
  beta_l1 <- 0
  beta_l2 <- 0
  beta_l10 <- 0
  beta_max <- 0
  beta_canberra <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n * d), n, d)
    Y <- matrix(rnorm(m * d), m, d)
    
    Y[, 1] <- Y[, 1] + Delta
    
    
    dist_mat_XX_l1 <- as.matrix(dist(X, method = "manhattan"))
    dist_mat_XY_l1 <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
    
    T_l1 <- outer(dist_mat_XX_l1, dist_mat_XY_l1,
                  FUN = function(x, y) as.numeric(x <= y))
    
    W_l1 <- (sum(T_l1[indices]) - 0.5 * length(T_l1[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l1[indices]) * (2 * 0.53 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_l2 <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY_l2 <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    T_l2 <- outer(dist_mat_XX_l2, dist_mat_XY_l2,
                  FUN = function(x, y) as.numeric(x <= y))
    
    W_l2 <- (sum(T_l2[indices]) - 0.5 * length(T_l2[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l2[indices]) * (2 * 0.53 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_l10 <- as.matrix(dist(X, method = "minkowski", p = 10))
    dist_mat_XY_l10 <- as.matrix(dist(rbind(X, Y), method = "minkowski", p = 10))
    
    T_l10 <- outer(dist_mat_XX_l10, dist_mat_XY_l10,
                   FUN = function(x, y) as.numeric(x <= y))
    
    W_l10 <- (sum(T_l10[indices]) - 0.5 * length(T_l10[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_l10[indices]) * (2 * 0.53 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_max <- as.matrix(dist(X, method = "maximum"))
    dist_mat_XY_max <- as.matrix(dist(rbind(X, Y), method = "maximum"))
    
    T_max <- outer(dist_mat_XX_max, dist_mat_XY_max,
                   FUN = function(x, y) as.numeric(x <= y))
    
    W_max <- (sum(T_max[indices]) - 0.5 * length(T_max[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_max[indices]) * (2 * 0.52 - 1) / 8)^(1/2)
    
    
    dist_mat_XX_can <- as.matrix(dist(X, method = "canberra"))
    dist_mat_XY_can <- as.matrix(dist(rbind(X, Y), method = "canberra"))
    
    T_can <- outer(dist_mat_XX_can, dist_mat_XY_can,
                   FUN = function(x, y) as.numeric(x <= y))
    
    W_can <- (sum(T_can[indices]) - 0.5 * length(T_can[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_can[indices]) * (2 * 0.51 - 1) / 8)^(1/2)
    
    
    p_l1 <- 2 * pnorm(-abs(W_l1))
    p_l2 <- 2 * pnorm(-abs(W_l2))
    p_l10 <- 2 * pnorm(-abs(W_l10))
    p_max <- 2 * pnorm(-abs(W_max))
    p_can <- 2 * pnorm(-abs(W_can))
    
    beta_l1 <- beta_l1 + as.numeric(p_l1 < alpha)
    beta_l2 <- beta_l2 + as.numeric(p_l2 < alpha)
    beta_l10 <- beta_l10 + as.numeric(p_l10 < alpha)
    beta_max <- beta_max + as.numeric(p_max < alpha)
    beta_canberra <- beta_canberra + as.numeric(p_can < alpha)
  }
  
  fig14_left$l1[k] <- beta_l1 / B
  fig14_left$l2[k] <- beta_l2 / B
  fig14_left$l10[k] <- beta_l10 / B
  fig14_left$max[k] <- beta_max / B
  fig14_left$Canberra[k] <- beta_canberra / B
  
  print(fig14_left[k, ])
}

###############################################################################
# Figure 14 (right):
# Estimated power for sample size 50 for different choices of p
# p = 1, 1.2, ..., 5, Inf
###############################################################################

n <- 50
m <- 50
d <- 100
Delta <- 0.5
P1_lp <- 0.53

p_grid <- c(seq(1, 5, by = 0.2), Inf)

fig14_right <- data.frame(
  p_plot = c(seq(1, 5, by = 0.2), 5.8),
  p_label = c(as.character(seq(1, 5, by = 0.2)), "Inf"),
  Power = NA_real_
)

indices <- generate_quadruples(n, m)

for (k in seq_along(p_grid)) {
  
  p_now <- p_grid[k]
  beta_p <- 0
  
  for (b in 1:B) {
    
    X <- matrix(rnorm(n * d), n, d)
    Y <- matrix(rnorm(m * d), m, d)
    
    Y[, 1] <- Y[, 1] + Delta
    
    if (is.infinite(p_now)) {
      dist_mat_XX <- as.matrix(dist(X, method = "maximum"))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "maximum"))
    } else if (p_now == 1) {
      dist_mat_XX <- as.matrix(dist(X, method = "manhattan"))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
    } else if (p_now == 2) {
      dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    } else {
      dist_mat_XX <- as.matrix(dist(X, method = "minkowski", p = p_now))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "minkowski", p = p_now))
    }
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    W <- (sum(T[indices]) - 0.5 * length(T[indices])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T[indices]) * (2 * P1_lp - 1) / 8)^(1/2)
    
    p_val <- 2 * pnorm(-abs(W))
    beta_p <- beta_p + as.numeric(p_val < alpha)
  }
  
  fig14_right$Power[k] <- beta_p / B
  print(fig14_right[k, ])
}

fig14_results <- list(
  metric_curves = fig14_left,
  metric_parameter = fig14_right
)


###############################################################################
# Figure 14
###############################################################################

fig14_left_plot  <- fig14_results$metric_curves
fig14_right_plot <- fig14_results$metric_parameter

fig14_left_legend_order <- c("l1", "l2", "l10", "max", "Canberra")

make_fig14_panel_left <- function(df, panel_title) {
  ggplot(df, aes(x = n, y = Power, color = Method, linetype = Method)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fig14_left_cols, breaks = fig14_left_legend_order) +
    scale_linetype_manual(values = fig14_left_types, breaks = fig14_left_legend_order) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, title = NULL),
      linetype = guide_legend(nrow = 1, byrow = TRUE, title = NULL)
    ) +
    scale_x_continuous(breaks = c(25, 50, 75, 100)) +
    scale_y_continuous(
      limits = c(0, 0.32),
      breaks = c(0.0, 0.1, 0.2, 0.3),
      labels = c("0.0", "0.1", "0.2", "0.3"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = panel_title,
      x = "Sample size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

make_fig14_panel_right <- function(df, panel_title) {
  df_finite <- df |> dplyr::filter(p_label != "Inf")
  df_inf    <- df |> dplyr::filter(p_label == "Inf")
  
  ggplot() +
    geom_line(
      data = df_finite,
      aes(x = p_plot, y = Power, group = 1),
      linewidth = 0.8,
      color = "blue"
    ) +
    geom_point(
      data = df_finite,
      aes(x = p_plot, y = Power),
      size = 2,
      color = "blue"
    ) +
    geom_point(
      data = df_inf,
      aes(x = p_plot, y = Power),
      size = 2,
      color = "blue"
    ) +
    scale_x_continuous(
      breaks = c(seq(1, 5, by = 0.4), df_inf$p_plot),
      labels = c(as.character(seq(1, 5, by = 0.4)), "Inf")
    ) +
    scale_y_continuous(
      limits = c(0.10, 0.205),
      breaks = c(0.100, 0.125, 0.150, 0.175, 0.200),
      labels = c("0.100", "0.125", "0.150", "0.175", "0.200"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = panel_title,
      x = "p",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig14_panel_left(
  plot_left_df,
  "Power curves for selected metrics"
)

p2 <- make_fig14_panel_right(
  fig14_right_plot,
  "Estimated power for different choices of p"
)

gridExtra::grid.arrange(p1, p2, nrow = 1)
