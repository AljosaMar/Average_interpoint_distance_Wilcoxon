################################################################################
# Script: Simulations.R
#
# Purpose:
# Reproduce the simulation experiments for the article
# "A TWO-SAMPLE TEST BASED ON AVERAGED WILCOXON RANK SUMS OVER INTERPOINT DISTANCES" 
# (https://arxiv.org/pdf/2408.10570), including empirical level and power studies 
# for the averaged interpoint Wilcoxon statistic under several distributional settings, 
# corresponding to Figures 3-8.
#
# Simulations included:
#   1) empirical level under the null hypothesis for normal, Cauchy, uniform,
#      and dependent normal settings,
#   2) power under location alternatives,
#   3) power under scale alternatives,
#   4) power comparisons in dependent and non-Gaussian settings,
#   5) comparison with benchmark procedures where applicable
#      (Hotelling's T2, Cramér test, and independent-distance Wilcoxon).
#
# Analysis overview:
# - For each setting, samples X and Y are generated under the specified model.
# - Interpoint distance matrices are computed using Euclidean or Manhattan
#   distance, depending on the scenario.
# - The averaged interpoint Wilcoxon statistic is evaluated for:
#     P1 = 1,
#     P1 = 2/3,
#     oracle P1 values specific to the simulation setting.
# - Rejection frequencies are estimated over Monte Carlo repetitions.
# - Results are collected into data frames and plotted as level/power curves.
#
# Required packages:
# - ggplot2
# - dplyr
# - tidyr
# - patchwork
# - ICSNP
# - cramer
# - MASS
#
# Required source file:
# - Utilities.R
#   This script assumes that helper functions such as
#   `generate_quadruples()` and `ind_wilcox_reject()` are available in the
#   R session or sourced from Utilities.R before execution.
#
# Important notes:
# - Random-number generation is controlled via `set.seed()`.
# - The script is organized by figure blocks corresponding to the supplemental
#   simulation study.
#
# Because the sample sizes are fixed within each simulation block, we do not
# call the general AW_interpoint_test() function from Utilities.R inside every
# Monte Carlo iteration. Instead, the index sets and normalizing constants are
# implemented directly in the simulation code to avoid repeated reconstruction
# and reduce computational overhead. For the same reason, the corresponding P1
# values are hard-coded for each setting.
#
# Outputs:
# - Data frames of estimated rejection rates under each simulation setting
# - Multi-panel plots for empirical level and power comparisons
# - Printed intermediate simulation summaries in the console
#
# Author: Aljosa Marjanovic
# Date: 10.03.2026
# Project: Two sample averaged Wilcoxon over interpoint distances
################################################################################

source("Utilities.R")


set.seed(1)

################################################################################
# Empirical level 
################################################################################

MC <- 100
ns <- c(10, 25, 50, 75, 100)

# Oracle values taken from Tables 5-9
P1_oracle_norm_l2   <- 0.538320
P1_oracle_norm_l1   <- 0.534561
P1_oracle_cauchy_l2 <- 0.571805
P1_oracle_cauchy_l1 <- 0.573040
P1_oracle_unif_l2   <- 0.517824
P1_oracle_dep_l2    <- 0.536000

Sigma_dep <- outer(1:10, 1:10, function(i, j) 1 / (1 + abs(i - j)))

null_results <- list()

# Monte Carlo design:
# - samples are balanced throughout (n = m)
# - each scenario is evaluated over the sample sizes in `ns`
# - rejection rates are estimated as the proportion of MC runs with rejection
# - unless stated otherwise (by dimensionality studies), the data dimension
#   is fixed within a block

################################################################################
# Normal, l2, dimension 10
################################################################################

# Test labels used throughout:
# - P11: AW test with variance bound P1 = 1
# - P12: AW test with variance bound P1 = 2/3
# - P13: AW test with scenario-specific oracle value for P1

df_normal_l2 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)


# We inline the AW statistic here instead of calling AW_interpoint_test()
# because n and m are fixed within each simulation block. This avoids
# repeatedly reconstructing the index set (calling generate_quadrupels(n,m)) 
# and normalizing terms inside every Monte Carlo iteration, significantly saving
# computation time.
 
for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- matrix(rnorm(n * 10, mean = 0, sd = 1), ncol = 10)
    Y <- matrix(rnorm(m * 10, mean = 0, sd = 1), ncol = 10)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_norm_l2 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  df_normal_l2$P11[k] <- beta_P11 / MC
  df_normal_l2$P12[k] <- beta_P12 / MC
  df_normal_l2$P13[k] <- beta_P13 / MC
  
  print("normal_l2")
  print(df_normal_l2[k, ])
  print("--------------------")
}

null_results$normal_l2 <- df_normal_l2

################################################################################
# Normal, l1, dimension 10
################################################################################

df_normal_l1 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- matrix(rnorm(n * 10, mean = 0, sd = 1), ncol = 10)
    Y <- matrix(rnorm(m * 10, mean = 0, sd = 1), ncol = 10)
    
    dist_mat_XX <- as.matrix(dist(X, method = "manhattan"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
   
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_norm_l1 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  df_normal_l1$P11[k] <- beta_P11 / MC
  df_normal_l1$P12[k] <- beta_P12 / MC
  df_normal_l1$P13[k] <- beta_P13 / MC
  
  print("normal_l1")
  print(df_normal_l1[k, ])
  print("--------------------")
}

null_results$normal_l1 <- df_normal_l1

################################################################################
# Cauchy, l2, dimension 10
################################################################################

df_cauchy_l2 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- matrix(rcauchy(n * 10, location = 0, scale = 1), ncol = 10)
    Y <- matrix(rcauchy(m * 10, location = 0, scale = 1), ncol = 10)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
   
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_cauchy_l2 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  df_cauchy_l2$P11[k] <- beta_P11 / MC
  df_cauchy_l2$P12[k] <- beta_P12 / MC
  df_cauchy_l2$P13[k] <- beta_P13 / MC
  
  print("cauchy_l2")
  print(df_cauchy_l2[k, ])
  print("--------------------")
}

null_results$cauchy_l2 <- df_cauchy_l2

################################################################################
# Cauchy, l1, dimension 10
################################################################################

df_cauchy_l1 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- matrix(rcauchy(n * 10, location = 0, scale = 1), ncol = 10)
    Y <- matrix(rcauchy(m * 10, location = 0, scale = 1), ncol = 10)
    
    dist_mat_XX <- as.matrix(dist(X, method = "manhattan"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "manhattan"))

    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_cauchy_l1 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  df_cauchy_l1$P11[k] <- beta_P11 / MC
  df_cauchy_l1$P12[k] <- beta_P12 / MC
  df_cauchy_l1$P13[k] <- beta_P13 / MC
  
  print("cauchy_l1")
  print(df_cauchy_l1[k, ])
  print("--------------------")
}

null_results$cauchy_l1 <- df_cauchy_l1

################################################################################
# Uniform on [0,1]^2, l2, with Hotelling and Cramer
################################################################################

df_uniform_l2 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  beta_hot <- 0
  beta_cra <- 0
  
  for (i in 1:MC) {
    X <- matrix(runif(n * 2, min = 0, max = 1), ncol = 2)
    Y <- matrix(runif(m * 2, min = 0, max = 1), ncol = 2)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
  
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_unif_l2 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
    
    beta_hot <- beta_hot + as.numeric(ICSNP::HotellingsT2(X, Y)$p.value < 0.05)
    beta_cra <- beta_cra + as.numeric(cramer::cramer.test(X, Y, replicates = 499)$p.value < 0.05)
  }
  
  df_uniform_l2$P11[k] <- beta_P11 / MC
  df_uniform_l2$P12[k] <- beta_P12 / MC
  df_uniform_l2$P13[k] <- beta_P13 / MC
  df_uniform_l2$hotelling[k] <- beta_hot / MC
  df_uniform_l2$cramer[k] <- beta_cra / MC
  
  print("uniform_l2")
  print(df_uniform_l2[k, ])
  print("--------------------")
}

null_results$uniform_l2 <- df_uniform_l2

################################################################################
# Dependent normal, l2, dimension 10
################################################################################

df_normal_dep_l2 <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- MASS::mvrnorm(n, mu = rep(0, 10), Sigma = Sigma_dep)
    Y <- MASS::mvrnorm(m, mu = rep(0, 10), Sigma = Sigma_dep)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_dep_l2 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  df_normal_dep_l2$P11[k] <- beta_P11 / MC
  df_normal_dep_l2$P12[k] <- beta_P12 / MC
  df_normal_dep_l2$P13[k] <- beta_P13 / MC
  
  print("normal_dep_l2")
  print(df_normal_dep_l2[k, ])
  print("--------------------")
}

null_results$normal_dep_l2 <- df_normal_dep_l2

null_results

################################################################################
# Null Graphs
################################################################################
fig3_df <- bind_rows(
  null_results$normal_l2      |> mutate(panel = "Level: Normal distribution in l2"),
  null_results$normal_l1      |> mutate(panel = "Level: Normal distribution in l1"),
  null_results$cauchy_l2      |> mutate(panel = "Level: Cauchy(0,1) distribution in l2"),
  null_results$cauchy_l1      |> mutate(panel = "Level: Cauchy(0,1) distribution in l1"),
  null_results$uniform_l2     |> mutate(panel = "Level: U(0,1) in l2"),
  null_results$normal_dep_l2  |> mutate(panel = "Level: Normal under dependency in l2")
) |>
  pivot_longer(
    cols = -c(n, panel),
    names_to = "test",
    values_to = "level"
  ) |>
  mutate(
    test = factor(
      test,
      levels = c("P13", "P12", "P11", "hotelling", "cramer"),
      labels = c("P1=0.54", "P1=2/3", "P1=1", "Hotellings T2", "Cramer")
    )
  )

cols_fig3 <- c(
  "P1=0.54" = "red",
  "P1=2/3" = "red",
  "P1=1" = "red",
  "Hotellings T2" = "black",
  "Cramer" = "green"
)

types_fig3 <- c(
  "P1=0.54" = "solid",
  "P1=2/3" = "dashed",
  "P1=1" = "dotted",
  "Hotellings T2" = "solid",
  "Cramer" = "solid"
)

make_fig3_panel <- function(panel_name) {
  
  df_sub <- filter(fig3_df, panel == panel_name)
  
  if (panel_name == "Level: U(0,1) in l2") {
    df_sub <- filter(df_sub, test %in% c("P1=0.54", "Hotellings T2", "Cramer"))
    legend_breaks <- c("P1=0.54", "Hotellings T2", "Cramer")
  } else {
    df_sub <- filter(df_sub, test %in% c("P1=0.54", "P1=2/3", "P1=1"))
    legend_breaks <- c("P1=0.54", "P1=2/3", "P1=1")
  }
  
  ggplot(df_sub, aes(x = n, y = level, color = test, linetype = test)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = cols_fig3, breaks = legend_breaks) +
    scale_linetype_manual(values = types_fig3, breaks = legend_breaks) +
    scale_x_continuous(breaks = c(25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 0.52), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    labs(
      title = panel_name,
      x = "Sample size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_fig3_panel("Level: Normal distribution in l2")
p2 <- make_fig3_panel("Level: Normal distribution in l1")
p3 <- make_fig3_panel("Level: Cauchy(0,1) distribution in l2")
p4 <- make_fig3_panel("Level: Cauchy(0,1) distribution in l1")
p5 <- make_fig3_panel("Level: U(0,1) in l2")
p6 <- make_fig3_panel("Level: Normal under dependency in l2")

((p1 | p2) / (p3 | p4) / (p5 | p6))

################################################################################
# Normal location and scale problems
################################################################################

dim <- 10
MC  <- 100
ns  <- c(25, 50, 75, 100)
mus <- c(0.3, 0.6, 0.9)

P1_oracle_l2_D10 <- 0.538   # Table 5, N(0_D,I_D), l2, d=10

results <- list()




for (mu in mus) {
  
  df_mu <- data.frame(
    n = ns,
    P11 = NA_real_,     # P_1      (P1=1)
    P12 = NA_real_,     # P_2/3    (P1=2/3)
    P13 = NA_real_,     # P1 oracle (P1=0.538320)
    wilcoxon_ind = NA_real_
  )
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- n
    
    indices <- generate_quadruples(n, m)
    
    beta_P11 <- 0
    beta_P12 <- 0
    beta_P13 <- 0
    beta_wil <- 0
    
    for (i in 1:MC) {
      X <- matrix(rnorm(dim * n, mean = 0,  sd = 1), ncol = dim)
      Y <- matrix(rnorm(dim * m, mean = mu, sd = 1), ncol = dim)
      
      # AW statistic via XX vs XY block 
      dist_mat_XX <- as.matrix(dist(X))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
   
      
      T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) { as.numeric(x <= y) })
      
      S <- sum(T[indices])
      L <- length(T[indices])
      
      # P_1 (P1=1)
      W1 <- (S - 0.5*L) / (( (n + m - 6) * (n - 3) * (n - 4) * L * (2*1 - 1) / 8 )^(1/2))
      beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
      
      # P_2/3
      W2 <- (S - 0.5*L) / (( (n + m - 6) * (n - 3) * (n - 4) * L * (2*(2/3) - 1) / 8 )^(1/2))
      beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
      
      # P_oracle (Table 5)
      Wo <- (S - 0.5*L) / (( (n + m - 6) * (n - 3) * (n - 4) * L * (2*P1_oracle_l2_D10 - 1) / 8 )^(1/2))
      beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
      
      # independent Wilcoxon competitor on independent distances
      beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
    }
    
    df_mu$P11[k] <- beta_P11 / MC
    df_mu$P12[k] <- beta_P12 / MC
    df_mu$P13[k] <- beta_P13 / MC
    df_mu$wilcoxon_ind[k] <- beta_wil / MC
    
    print(paste("mu =", mu, " n =", n, " done at ", Sys.time()))
    print(df_mu[k, ])
    print("--------------------")
  }
  
  results[[paste0("mu_", mu)]] <- df_mu
}

# results now contains three data.frames:
results


# SCALE (variance) alternatives
scales <- c(0.8, 0.6, 0.4)

results_scale <- list()

for (s in scales) {
  
  df_s <- data.frame(
    n = ns,
    P11 = NA_real_,
    P12 = NA_real_,
    P13 = NA_real_,
    wilcoxon_ind = NA_real_
  )
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- n
    
    indices <- generate_quadruples(n, m)
    
    beta_P11 <- 0
    beta_P12 <- 0
    beta_P13 <- 0
    beta_wil <- 0
    
    for (i in 1:MC) {
      X <- matrix(rnorm(dim * n, mean = 0, sd = 1), ncol = dim)
      Y <- matrix(rnorm(dim * m, mean = 0, sd = sqrt(s)), ncol = dim)
      
      dist_mat_XX <- as.matrix(dist(X))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y)))
    
      
      T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) { as.numeric(x <= y) })
      
      S <- sum(T[indices])
      L <- length(T[indices])
      
      W1 <- (S - 0.5*L) / (((n + m - 6) * (n - 3) * (n - 4) * L * (2*1 - 1) / 8)^(1/2))
      beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
      
      W2 <- (S - 0.5*L) / (((n + m - 6) * (n - 3) * (n - 4) * L * (2*(2/3) - 1) / 8)^(1/2))
      beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
      
      Wo <- (S - 0.5*L) / (((n + m - 6) * (n - 3) * (n - 4) * L * (2*P1_oracle_l2_D10 - 1) / 8)^(1/2))
      beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
      
      beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
    }
    
    df_s$P11[k] <- beta_P11 / MC
    df_s$P12[k] <- beta_P12 / MC
    df_s$P13[k] <- beta_P13 / MC
    df_s$wilcoxon_ind[k] <- beta_wil / MC
  }
  
  results_scale[[paste0("s_", s)]] <- df_s
}


results_scale
################################################################################
# Construct plot_df for Figures 4-5 block
# (works with the existing results and results_scale objects in Simulations.R)
################################################################################
plot_df <- bind_rows(
  results$mu_0.3 |> mutate(effect = "small effect",  problem = "location problem"),
  results$mu_0.6 |> mutate(effect = "medium effect", problem = "location problem"),
  results$mu_0.9 |> mutate(effect = "large effect",  problem = "location problem"),
  results_scale$s_0.8 |> mutate(effect = "small effect",  problem = "scale problem"),
  results_scale$s_0.6 |> mutate(effect = "medium effect", problem = "scale problem"),
  results_scale$s_0.4 |> mutate(effect = "large effect",  problem = "scale problem")
) |>
  pivot_longer(
    cols = c(P13, P12, P11, wilcoxon_ind),
    names_to = "test",
    values_to = "power"
  ) |>
  mutate(
    test = recode(test,
                  P13 = "P1=0.54",
                  P12 = "P1=2/3",
                  P11 = "P1=1",
                  wilcoxon_ind = "Wilcoxon"),
    test = factor(test,
                  levels = c("P1=0.54","P1=2/3","P1=1","Wilcoxon"))
  )

################################################################################
# Figures 4-5
################################################################################

cols  <- c("P1=0.54" = "red", "P1=2/3" = "red", "P1=1" = "red", "Wilcoxon" = "blue")
types <- c("P1=0.54" = "solid", "P1=2/3" = "dashed", "P1=1" = "dotted", "Wilcoxon" = "dotdash")

make_panel <- function(eff, prob) {
  ggplot(filter(plot_df, effect == eff, problem == prob),
         aes(x = n, y = power, color = test, linetype = test)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = cols) +
    scale_linetype_manual(values = types) +
    scale_x_continuous(breaks = c(25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = paste("Power:", prob, eff),
      x = "Sample size",
      y = "Rejection rate"
    ) +
    theme_gray(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_line(linewidth = 0.25)
    )
}

p1 <- make_panel("small effect",  "location problem")
p2 <- make_panel("small effect",  "scale problem")
p3 <- make_panel("medium effect", "location problem")
p4 <- make_panel("medium effect", "scale problem")
p5 <- make_panel("large effect",  "location problem")
p6 <- make_panel("large effect",  "scale problem")

((p1 | p2) / (p3 | p4) / (p5 | p6)) 


################################################################################
# Power: Dependent dimensions and uniform distribution (Figure 6)
################################################################################

ns <- c(25, 50, 75, 100)

P1_oracle_dep_l2  <- 0.536
P1_oracle_unif_l2 <- 0.52

Sigma_dep <- outer(1:10, 1:10, function(i, j) 1 / (1 + abs(i - j)))

###############################################################################
# Left panel: dependent normal location problem in R^10
# X ~ N(0, Sigma), Y ~ N(DeltaM, Sigma), DeltaM = (0.6, ..., 0.6)
###############################################################################

fig6_dep <- data.frame(
  n = ns,
  P11 = NA_real_,
  P12 = NA_real_,
  P13 = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P11 <- 0
  beta_P12 <- 0
  beta_P13 <- 0
  
  for (i in 1:MC) {
    X <- MASS::mvrnorm(n, mu = rep(0, 10),   Sigma = Sigma_dep)
    Y <- MASS::mvrnorm(m, mu = rep(0.6, 10), Sigma = Sigma_dep)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))

    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W1 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 1 - 1) / 8)^(1/2))
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * (2/3) - 1) / 8)^(1/2))
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_dep_l2 - 1) / 8)^(1/2))
    
    beta_P11 <- beta_P11 + as.numeric(pnorm(abs(W1), lower.tail = FALSE) < 0.025)
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
  }
  
  fig6_dep$P11[k] <- beta_P11 / MC
  fig6_dep$P12[k] <- beta_P12 / MC
  fig6_dep$P13[k] <- beta_P13 / MC
  
  print("fig6_dep")
  print(fig6_dep[k, ])
  print("--------------------")
}

fig6_dep


###############################################################################
# Right panel: uniform location problem in R^2
# X ~ U([0,1]^2), Y ~ U([0.2,1.2]^2)
# Compare only AW-oracle, Hotelling, and Cramer
###############################################################################

fig6_unif <- data.frame(
  n = ns,
  P13 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_
)

for (k in 1:length(ns)) {
  n <- ns[k]
  m <- n
  
  indices <- generate_quadruples(n, m)
  
  beta_P13 <- 0
  beta_hot <- 0
  beta_cra <- 0
  
  for (i in 1:MC) {
    X <- matrix(runif(n * 2, min = 0.0, max = 1.0), ncol = 2)
    Y <- matrix(runif(m * 2, min = 0.2, max = 1.2), ncol = 2)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))

    
    T <- outer(dist_mat_XX, dist_mat_XY, FUN = function(x, y) as.numeric(x <= y))
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    Wo <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * P1_oracle_unif_l2 - 1) / 8)^(1/2))
    
    beta_P13 <- beta_P13 + as.numeric(pnorm(abs(Wo), lower.tail = FALSE) < 0.025)
    beta_hot <- beta_hot + as.numeric(ICSNP::HotellingsT2(X, Y)$p.value < 0.05)
    beta_cra <- beta_cra + as.numeric(cramer::cramer.test(X, Y, replicates = 499)$p.value < 0.05)
  }
  
  fig6_unif$P13[k] <- beta_P13 / MC
  fig6_unif$hotelling[k] <- beta_hot / MC
  fig6_unif$cramer[k] <- beta_cra / MC
  
  print("fig6_unif")
  print(fig6_unif[k, ])
  print("--------------------")
}

fig6_unif


fig6_results <- list(
  dependent_location = fig6_dep,
  uniform_location = fig6_unif
)

fig6_results

###############################################################################
# Plot Figure 6
###############################################################################

fig6_dep_plot <- fig6_dep |>
  pivot_longer(
    cols = c(P13, P12, P11),
    names_to = "test",
    values_to = "power"
  ) |>
  mutate(
    test = factor(
      test,
      levels = c("P13", "P12", "P11"),
      labels = c("P1=0.535", "P1=2/3", "P1=1")
    )
  )

cols_dep <- c(
  "P1=0.535" = "red",
  "P1=2/3"   = "red",
  "P1=1"     = "red"
)

types_dep <- c(
  "P1=0.535" = "solid",
  "P1=2/3"   = "dashed",
  "P1=1"     = "dotted"
)

p_dep <- ggplot(fig6_dep_plot,
                aes(x = n, y = power, color = test, linetype = test)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_dep) +
  scale_linetype_manual(values = types_dep) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(
    title = "Power: Dependent dimensions",
    x = "Sample size",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

###############################################################################
# Right panel: uniform location problem
###############################################################################

fig6_unif_plot <- fig6_unif |>
  pivot_longer(
    cols = c(P13, hotelling, cramer),
    names_to = "test",
    values_to = "power"
  ) |>
  mutate(
    test = factor(
      test,
      levels = c("P13", "hotelling", "cramer"),
      labels = c("P1=0.52", "hotelling", "cramer")
    )
  )

cols_unif <- c(
  "P1=0.52"   = "red",
  "hotelling" = "black",
  "cramer"    = "green"
)

types_unif <- c(
  "P1=0.52"   = "dashed",
  "hotelling" = "solid",
  "cramer"    = "dotted"
)

p_unif <- ggplot(fig6_unif_plot,
                 aes(x = n, y = power, color = test, linetype = test)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_unif) +
  scale_linetype_manual(values = types_unif) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(
    title = "Power: Uniform distribution",
    x = "Sample size",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

###############################################################################
# Combined figure: 1 row, 2 columns
###############################################################################

(p_dep | p_unif)

################################################################################
# Dimensionality simulations, comparison with Hotelling and Cramer
################################################################################

MC <- 300

dims_top <- c(5, 10, 50, 200)
ns_top   <- c(25, 50, 75, 100)

dims_bottom <- c(10, 25, 50, 75, 100, 125, 150)
n_fixed <- 50
m_fixed <- 50

################################################################################
# Top-left: normal location problem, one shifted component
# X ~ N(0_D, I_D)
# Y ~ N((0.6, 0, ..., 0), I_D)
# x-axis = sample size
################################################################################

fig8_top_location <- list()

for (dim in dims_top) {
  
  df_dim <- data.frame(
    n = ns_top,
    P12 = NA_real_
  )
  
  for (k in 1:length(ns_top)) {
    n <- ns_top[k]
    m <- n
    
    indices <- generate_quadruples(n, m)
    
    beta_P12 <- 0
    
    for (i in 1:MC) {
      X <- matrix(rnorm(dim * n, mean = 0, sd = 1), ncol = dim)
      Y <- matrix(rnorm(dim * m, mean = 0, sd = 1), ncol = dim)
      Y[, 1] <- Y[, 1] + 0.6
      
      dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
   
      
      T <- outer(dist_mat_XX, dist_mat_XY,
                 FUN = function(x, y) { as.numeric(x <= y) })
      
      S <- sum(T[indices])
      L <- length(T[indices])
      
      W2 <- (S - 0.5 * L) /
        (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.54 - 1) / 8)^(1/2))
      
      beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    }
    
    df_dim$P12[k] <- beta_P12 / MC
    
    print(paste("fig8 top-left: dim =", dim, "n =", n))
    print(df_dim[k, ])
    print("--------------------")
  }
  
  fig8_top_location[[paste0("dim_", dim)]] <- df_dim
}

fig8_top_location


################################################################################
# Top-right: normal scale problem, one scaled component
# X ~ N(0_D, I_D)
# Y ~ N(0_D, Sigma^D_0.5), Sigma^D_0.5 = diag(0.5, 1, ..., 1)
# x-axis = sample size
################################################################################

fig8_top_scale <- list()

for (dim in dims_top) {
  
  Sigma_scale <- diag(c(0.5, rep(1, dim - 1)))
  
  df_dim <- data.frame(
    n = ns_top,
    P12 = NA_real_
  )
  
  for (k in 1:length(ns_top)) {
    n <- ns_top[k]
    m <- n
    
    indices <- generate_quadruples(n, m)
    
    beta_P12 <- 0
    
    for (i in 1:MC) {
      X <- matrix(rnorm(dim * n, mean = 0, sd = 1), ncol = dim)
      Y <- MASS::mvrnorm(m, mu = rep(0, dim), Sigma = Sigma_scale)
      
      dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
      dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
      
      T <- outer(dist_mat_XX, dist_mat_XY,
                 FUN = function(x, y) { as.numeric(x <= y) })
      
      S <- sum(T[indices])
      L <- length(T[indices])
      
      W2 <- (S - 0.5 * L) /
        (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.54 - 1) / 8)^(1/2))
      
      beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    }
    
    df_dim$P12[k] <- beta_P12 / MC
    
    print(paste("fig8 top-right: dim =", dim, "n =", n))
    print(df_dim[k, ])
    print("--------------------")
  }
  
  fig8_top_scale[[paste0("dim_", dim)]] <- df_dim
}

fig8_top_scale


################################################################################
# Bottom-left: normal location problem, one shifted component
# X ~ N(0_D, I_D)
# Y ~ N((0.6, 0, ..., 0), I_D)
# x-axis = dimension, n = m = 50 fixed
# compare: P12, Hotelling, Cramer, independent Wilcoxon
################################################################################

fig8_bottom_normal <- data.frame(
  dimension = dims_bottom,
  P12 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_,
  wilcoxon_ind = NA_real_
)

for (k in 1:length(dims_bottom)) {
  dim <- dims_bottom[k]
  n <- n_fixed
  m <- m_fixed
  
  indices <- generate_quadruples(n, m)
  
  beta_P12 <- 0
  beta_hot <- 0
  beta_cra <- 0
  beta_wil <- 0
  
  for (i in 1:MC) {
    X <- matrix(rnorm(dim * n, mean = 0, sd = 1), ncol = dim)
    Y <- matrix(rnorm(dim * m, mean = 0, sd = 1), ncol = dim)
    Y[, 1] <- Y[, 1] + 0.9
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
  
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) { as.numeric(x <= y) })
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.534 - 1) / 8)^(1/2))
    
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    
    hot_p <- tryCatch(
      ICSNP::HotellingsT2(X, Y)$p.value,
      error = function(e) NA_real_
    )
    beta_hot <- beta_hot + ifelse(length(hot_p) == 0 || is.na(hot_p) || is.null(hot_p), 0,
                                  as.numeric(hot_p < 0.05))
    
    cra_p <- tryCatch(
      cramer::cramer.test(X, Y, replicates = 499)$p.value,
      error = function(e) NA_real_
    )
    beta_cra <- beta_cra + ifelse(length(cra_p) == 0 || is.na(cra_p) || is.null(cra_p), 0,
                                  as.numeric(cra_p < 0.05))
    
    beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
  }
  
  fig8_bottom_normal$P12[k] <- beta_P12 / MC
  fig8_bottom_normal$hotelling[k] <- beta_hot / MC
  fig8_bottom_normal$cramer[k] <- beta_cra / MC
  fig8_bottom_normal$wilcoxon_ind[k] <- beta_wil / MC
  
  print("fig8 bottom-left")
  print(fig8_bottom_normal[k, ])
  print("--------------------")
}

fig8_bottom_normal


################################################################################
# Bottom-right: Cauchy location problem, one shifted component
# X ~ Cauchy(0,1)^D
# Y ~ Cauchy(1,1) x Cauchy(0,1)^(D-1)
# x-axis = dimension, n = m = 50 fixed
# compare: P12, Hotelling, Cramer, independent Wilcoxon
################################################################################

fig8_bottom_cauchy <- data.frame(
  dimension = dims_bottom,
  P12 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_,
  wilcoxon_ind = NA_real_
)

for (k in 1:length(dims_bottom)) {
  dim <- dims_bottom[k]
  n <- n_fixed
  m <- m_fixed
  
  indices <- generate_quadruples(n, m)
  
  beta_P12 <- 0
  beta_hot <- 0
  beta_cra <- 0
  beta_wil <- 0
  
  for (i in 1:MC) {
    X <- matrix(rcauchy(dim * n, location = 0, scale = 1), ncol = dim)
    Y <- matrix(rcauchy(dim * m, location = 0, scale = 1), ncol = dim)
    Y[, 1] <- rcauchy(m, location = 1, scale = 1)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
  
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) { as.numeric(x <= y) })
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.57 - 1) / 8)^(1/2))
    
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    
    hot_p <- tryCatch(
      ICSNP::HotellingsT2(X, Y)$p.value,
      error = function(e) NA_real_
    )
    beta_hot <- beta_hot + ifelse(length(hot_p) == 0 || is.na(hot_p) || is.null(hot_p), 0,
                                  as.numeric(hot_p < 0.05))
    
    cra_p <- tryCatch(
      cramer::cramer.test(X, Y, replicates = 499)$p.value,
      error = function(e) NA_real_
    )
    beta_cra <- beta_cra + ifelse(length(cra_p) == 0 || is.na(cra_p) || is.null(cra_p), 0,
                                  as.numeric(cra_p < 0.05))
    
    beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
  }
  
  fig8_bottom_cauchy$P12[k] <- beta_P12 / MC
  fig8_bottom_cauchy$hotelling[k] <- beta_hot / MC
  fig8_bottom_cauchy$cramer[k] <- beta_cra / MC
  fig8_bottom_cauchy$wilcoxon_ind[k] <- beta_wil / MC
  
  print("fig8 bottom-right")
  print(fig8_bottom_cauchy[k, ])
  print("--------------------")
}

fig8_bottom_cauchy

################################################################################
# Bottom-left: normal location problem, first 10 shifted components
# X ~ N(0_D, I_D)
# Y: first 10 components shifted by +0.6, remaining components unchanged
# x-axis = dimension, n = m = 50 fixed
# compare: P12, Hotelling, Cramer, independent Wilcoxon
################################################################################

fig8_bottom_normal <- data.frame(
  dimension = dims_bottom,
  P12 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_,
  wilcoxon_ind = NA_real_
)

for (k in 1:length(dims_bottom)) {
  dim <- dims_bottom[k]
  n <- n_fixed
  m <- m_fixed
  
  indices <- generate_quadruples(n, m)
  
  beta_P12 <- 0
  beta_hot <- 0
  beta_cra <- 0
  beta_wil <- 0
  
  for (i in 1:MC) {
    X <- matrix(rnorm(dim * n, mean = 0, sd = 1), ncol = dim)
    Y <- matrix(rnorm(dim * m, mean = 0, sd = 1), ncol = dim)
    Y[, 1:10] <- Y[, 1:10] + 0.5
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) { as.numeric(x <= y) })
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.54 - 1) / 8)^(1/2))
    
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    
    hot_p <- tryCatch(
      ICSNP::HotellingsT2(X, Y)$p.value,
      error = function(e) NA_real_
    )
    beta_hot <- beta_hot + ifelse(length(hot_p) == 0 || is.na(hot_p) || is.null(hot_p), 0,
                                  as.numeric(hot_p < 0.05))
    
    cra_p <- tryCatch(
      cramer::cramer.test(X, Y, replicates = 499)$p.value,
      error = function(e) NA_real_
    )
    beta_cra <- beta_cra + ifelse(length(cra_p) == 0 || is.na(cra_p) || is.null(cra_p), 0,
                                  as.numeric(cra_p < 0.05))
    
    beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
  }
  
  fig8_bottom_normal$P12[k] <- beta_P12 / MC
  fig8_bottom_normal$hotelling[k] <- beta_hot / MC
  fig8_bottom_normal$cramer[k] <- beta_cra / MC
  fig8_bottom_normal$wilcoxon_ind[k] <- beta_wil / MC
  
  print("fig8 bottom-left")
  print(fig8_bottom_normal[k, ])
  print("--------------------")
}

fig8_bottom_normal


################################################################################
# Bottom-right: Cauchy location problem, first 10 shifted components
# X ~ Cauchy(0,1)^D
# Y: first 10 components are Cauchy(1,1), remaining components Cauchy(0,1)
# x-axis = dimension, n = m = 50 fixed
# compare: P12, Hotelling, Cramer, independent Wilcoxon
################################################################################

fig8_bottom_cauchy <- data.frame(
  dimension = dims_bottom,
  P12 = NA_real_,
  hotelling = NA_real_,
  cramer = NA_real_,
  wilcoxon_ind = NA_real_
)

for (k in 1:length(dims_bottom)) {
  dim <- dims_bottom[k]
  n <- n_fixed
  m <- m_fixed
  
  indices <- generate_quadruples(n, m)
  
  beta_P12 <- 0
  beta_hot <- 0
  beta_cra <- 0
  beta_wil <- 0
  
  for (i in 1:MC) {
    X <- matrix(rcauchy(dim * n, location = 0, scale = 1), ncol = dim)
    Y <- matrix(rcauchy(dim * m, location = 0, scale = 1), ncol = dim)
    Y[, 1:10] <- matrix(rcauchy(m * 10, location = 1, scale = 1), ncol = 10)
    
    dist_mat_XX <- as.matrix(dist(X, method = "euclidean"))
    dist_mat_XY <- as.matrix(dist(rbind(X, Y), method = "euclidean"))
    
    T <- outer(dist_mat_XX, dist_mat_XY,
               FUN = function(x, y) { as.numeric(x <= y) })
    
    S <- sum(T[indices])
    L <- length(T[indices])
    
    W2 <- (S - 0.5 * L) /
      (((n + m - 6) * (n - 3) * (n - 4) * L * (2 * 0.57 - 1) / 8)^(1/2))
    
    beta_P12 <- beta_P12 + as.numeric(pnorm(abs(W2), lower.tail = FALSE) < 0.025)
    
    hot_p <- tryCatch(
      ICSNP::HotellingsT2(X, Y)$p.value,
      error = function(e) NA_real_
    )
    beta_hot <- beta_hot + ifelse(length(hot_p) == 0 || is.na(hot_p) || is.null(hot_p), 0,
                                  as.numeric(hot_p < 0.05))
    
    cra_p <- tryCatch(
      cramer::cramer.test(X, Y, replicates = 499)$p.value,
      error = function(e) NA_real_
    )
    beta_cra <- beta_cra + ifelse(length(cra_p) == 0 || is.na(cra_p) || is.null(cra_p), 0,
                                  as.numeric(cra_p < 0.05))
    
    beta_wil <- beta_wil + ind_wilcox_reject(X, Y)
  }
  
  fig8_bottom_cauchy$P12[k] <- beta_P12 / MC
  fig8_bottom_cauchy$hotelling[k] <- beta_hot / MC
  fig8_bottom_cauchy$cramer[k] <- beta_cra / MC
  fig8_bottom_cauchy$wilcoxon_ind[k] <- beta_wil / MC
  
  print("fig8 bottom-right")
  print(fig8_bottom_cauchy[k, ])
  print("--------------------")
}

fig8_bottom_cauchy

################################################################################

fig8_results <- list(
  top_location = fig8_top_location,
  top_scale = fig8_top_scale,
  bottom_normal = fig8_bottom_normal,
  bottom_cauchy = fig8_bottom_cauchy
)

fig8_results

################################################################################
# Figure 8
################################################################################
# Top-left: location, x-axis = sample size
################################################################################

fig8_top_loc_plot <- bind_rows(lapply(names(fig8_top_location), function(nm) {
  fig8_top_location[[nm]] |>
    mutate(dimension_group = sub("dim_", "dim", nm))
}))

fig8_top_loc_plot <- fig8_top_loc_plot |>
  mutate(
    dimension_group = factor(
      dimension_group,
      levels = c("dim5", "dim10", "dim50", "dim200")
    )
  )

cols_top <- c(
  "dim5"   = "red",
  "dim10"  = "green",
  "dim50"  = "black",
  "dim200" = "blue"
)

types_top <- c(
  "dim5"   = "solid",
  "dim10"  = "dashed",
  "dim50"  = "dotted",
  "dim200" = "dotdash"
)

p1 <- ggplot(
  fig8_top_loc_plot,
  aes(x = n, y = P12, color = dimension_group, linetype = dimension_group)
) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_top) +
  scale_linetype_manual(values = types_top) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    labels = function(x) sprintf("%.2f", x)
  ) +
  labs(
    title = "Dimensionality: location problem",
    x = "Sample size",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

################################################################################
# Top-right: scale, x-axis = sample size
################################################################################

fig8_top_scale_plot <- bind_rows(lapply(names(fig8_top_scale), function(nm) {
  fig8_top_scale[[nm]] |>
    mutate(dimension_group = sub("dim_", "dim", nm))
}))

fig8_top_scale_plot <- fig8_top_scale_plot |>
  mutate(
    dimension_group = factor(
      dimension_group,
      levels = c("dim5", "dim10", "dim50", "dim200")
    )
  )

p2 <- ggplot(
  fig8_top_scale_plot,
  aes(x = n, y = P12, color = dimension_group, linetype = dimension_group)
) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_top) +
  scale_linetype_manual(values = types_top) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    labels = function(x) sprintf("%.2f", x)
  ) +
  labs(
    title = "Dimensionality: scale problem",
    x = "Sample size",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

################################################################################
# Bottom-left: normal location, x-axis = dimension
################################################################################

fig8_bottom_normal_plot <- fig8_bottom_normal |>
  pivot_longer(
    cols = c(P12, hotelling, cramer, wilcoxon_ind),
    names_to = "test",
    values_to = "power"
  ) |>
  mutate(
    test = recode(
      test,
      P12 = "P1 = 2/3",
      hotelling = "hotelling",
      cramer = "cramer",
      wilcoxon_ind = "wil"
    ),
    test = factor(
      test,
      levels = c("cramer", "hotelling", "wil", "P1 = 2/3")
    )
  )

cols_bottom <- c(
  "cramer"    = "red",
  "hotelling" = "green",
  "wil"       = "black",
  "P1 = 2/3"  = "blue"
)

types_bottom <- c(
  "cramer"    = "solid",
  "hotelling" = "dashed",
  "wil"       = "dotted",
  "P1 = 2/3"  = "dotdash"
)

p3 <- ggplot(
  fig8_bottom_normal_plot,
  aes(x = dimension, y = power, color = test, linetype = test)
) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_bottom) +
  scale_linetype_manual(values = types_bottom) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 100, 125, 150)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(
    title = "Power degradation: normal location",
    x = "Dimension",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

################################################################################
# Bottom-right: Cauchy location, x-axis = dimension
################################################################################

fig8_bottom_cauchy_plot <- fig8_bottom_cauchy |>
  pivot_longer(
    cols = c(P12, hotelling, cramer, wilcoxon_ind),
    names_to = "test",
    values_to = "power"
  ) |>
  mutate(
    test = recode(
      test,
      P12 = "P1 = 2/3",
      hotelling = "hotelling",
      cramer = "cramer",
      wilcoxon_ind = "wil"
    ),
    test = factor(
      test,
      levels = c("cramer", "hotelling", "wil", "P1 = 2/3")
    )
  )

p4 <- ggplot(
  fig8_bottom_cauchy_plot,
  aes(x = dimension, y = power, color = test, linetype = test)
) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_bottom) +
  scale_linetype_manual(values = types_bottom) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 100, 125, 150)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(
    title = "Power degradation: Cauchy location",
    x = "Dimension",
    y = "Rejection rate"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.25)
  )

((p1 | p2) / (p3 | p4))

