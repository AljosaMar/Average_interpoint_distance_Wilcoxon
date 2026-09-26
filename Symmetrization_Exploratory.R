

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

source("Utilities_v2.R")

################################################################################
# Fig. 15: Symmetrizations
# AW = original averaged Wilcoxon test
# Sym1 = statistic (22), permutation test
# Sym2 = statistic (23), permutation test
################################################################################

set.seed(1)

MC <- 1000
B_perm <- 250
ns <- c(25, 50, 75, 100)
alpha <- 0.05
P1_unif <- 0.52

################################################################################
# Left panel: concentric spheres
# X ~ U(S^2(0,1)), Y ~ U(S^2(0,1.1))
################################################################################

fig15_left <- data.frame(
  n = ns,
  AW = NA_real_,
  Sym1 = NA_real_,
  Sym2 = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  
  # AW indices: same style as in Further_Comparisons.R
  indices_aw_xy <- generate_quadruples(n, m)
  indices_aw_yx <- generate_quadruples(m, n)
  
  # For Sym1:
  # ordered pairs (i1,i2) from X excluding i5
  pairs_X_excl <- vector("list", n)
  for (i5 in 1:n) {
    idx <- setdiff(1:n, i5)
    tmp <- as.matrix(expand.grid(i1 = idx, i2 = idx))
    pairs_X_excl[[i5]] <- tmp[tmp[, 1] != tmp[, 2], , drop = FALSE]
  }
  
  # ordered pairs (i3,i4) from Y excluding i6
  pairs_Y_excl <- vector("list", m)
  for (i6 in 1:m) {
    idx <- setdiff(1:m, i6)
    tmp <- as.matrix(expand.grid(i3 = idx, i4 = idx))
    pairs_Y_excl[[i6]] <- tmp[tmp[, 1] != tmp[, 2], , drop = FALSE]
  }
  
  denom_sym1 <- n * (n - 1) * (n - 2) * m * (m - 1) * (m - 2)
  
  beta_AW <- 0
  beta_Sym1 <- 0
  beta_Sym2 <- 0
  
  for (i in 1:MC) {
    X <- generate_uniform_sphere(n, 3, radius = 1)
    Y <- generate_uniform_sphere(m, 3, radius = 1.1)
    
    ###########################################################################
    # AW (original test, asymptotic p-value)
    ###########################################################################
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_all_xy <- as.matrix(dist(rbind(X, Y)))
    
    T_aw_xy <- outer(
      dist_mat_XX,
      dist_mat_all_xy,
      FUN = function(x, y) as.numeric(x <= y)
    )
    
    W_aw <- (sum(T_aw_xy[indices_aw_xy]) - 0.5 * length(T_aw_xy[indices_aw_xy])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_aw_xy[indices_aw_xy]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    p_aw <- 2 * pnorm(-abs(W_aw))
    
    ###########################################################################
    # Sym1 observed statistic
    # T_s1 = average of I{ d(Xi1,Xi2) + d(Yi3,Yi4) < 2 d(Xi5,Yi6) }
    ###########################################################################
    
    dist_mat_YY <- as.matrix(dist(Y))
    dist_mat_XY <- dist_mat_all_xy[1:n, (n + 1):(n + m), drop = FALSE]
    
    dxx_excl <- vector("list", n)
    for (i5 in 1:n) {
      dxx_excl[[i5]] <- dist_mat_XX[pairs_X_excl[[i5]]]
    }
    
    dyy_excl_sorted <- vector("list", m)
    for (i6 in 1:m) {
      dyy_excl_sorted[[i6]] <- sort(dist_mat_YY[pairs_Y_excl[[i6]]], method = "quick")
    }
    
    T_sym1_obs <- 0
    for (i5 in 1:n) {
      for (i6 in 1:m) {
        thr <- 2 * dist_mat_XY[i5, i6]
        T_sym1_obs <- T_sym1_obs +
          sum(findInterval(thr - dxx_excl[[i5]], dyy_excl_sorted[[i6]], rightmost.closed = TRUE))
      }
    }
    T_sym1_obs <- T_sym1_obs / denom_sym1
    
    ###########################################################################
    # Sym2 observed statistic
    # T_s2 = 1/2 { T(X,Y) + T(Y,X) }
    ###########################################################################
    
    dist_mat_all_yx <- as.matrix(dist(rbind(Y, X)))
    
    T_aw_yx <- outer(
      dist_mat_YY,
      dist_mat_all_yx,
      FUN = function(x, y) as.numeric(x <= y)
    )
    
    T_xy_obs <- sum(T_aw_xy[indices_aw_xy]) / length(T_aw_xy[indices_aw_xy])
    T_yx_obs <- sum(T_aw_yx[indices_aw_yx]) / length(T_aw_yx[indices_aw_yx])
    T_sym2_obs <- 0.5 * (T_xy_obs + T_yx_obs)
    
    ###########################################################################
    # Permutation p-values for Sym1 and Sym2
    ###########################################################################
    
    Z <- rbind(X, Y)
    T_sym1_perm <- numeric(B_perm)
    T_sym2_perm <- numeric(B_perm)
    
    for (b in 1:B_perm) {
      id <- sample.int(n + m, n + m, replace = FALSE)
      Xb <- Z[id[1:n], , drop = FALSE]
      Yb <- Z[id[(n + 1):(n + m)], , drop = FALSE]
      
      dist_mat_XX_b <- as.matrix(dist(Xb))
      dist_mat_YY_b <- as.matrix(dist(Yb))
      
      dist_mat_all_xy_b <- as.matrix(dist(rbind(Xb, Yb)))
      dist_mat_XY_b <- dist_mat_all_xy_b[1:n, (n + 1):(n + m), drop = FALSE]
      
      # Sym1 permutation statistic
      dxx_excl_b <- vector("list", n)
      for (i5 in 1:n) {
        dxx_excl_b[[i5]] <- dist_mat_XX_b[pairs_X_excl[[i5]]]
      }
      
      dyy_excl_sorted_b <- vector("list", m)
      for (i6 in 1:m) {
        dyy_excl_sorted_b[[i6]] <- sort(dist_mat_YY_b[pairs_Y_excl[[i6]]], method = "quick")
      }
      
      T_sym1_b <- 0
      for (i5 in 1:n) {
        for (i6 in 1:m) {
          thr_b <- 2 * dist_mat_XY_b[i5, i6]
          T_sym1_b <- T_sym1_b +
            sum(findInterval(thr_b - dxx_excl_b[[i5]], dyy_excl_sorted_b[[i6]], rightmost.closed = TRUE))
        }
      }
      T_sym1_perm[b] <- T_sym1_b / denom_sym1
      
      # Sym2 permutation statistic
      T_aw_xy_b <- outer(
        dist_mat_XX_b,
        dist_mat_all_xy_b,
        FUN = function(x, y) as.numeric(x <= y)
      )
      
      dist_mat_all_yx_b <- as.matrix(dist(rbind(Yb, Xb)))
      T_aw_yx_b <- outer(
        dist_mat_YY_b,
        dist_mat_all_yx_b,
        FUN = function(x, y) as.numeric(x <= y)
      )
      
      T_xy_b <- sum(T_aw_xy_b[indices_aw_xy]) / length(T_aw_xy_b[indices_aw_xy])
      T_yx_b <- sum(T_aw_yx_b[indices_aw_yx]) / length(T_aw_yx_b[indices_aw_yx])
      T_sym2_perm[b] <- 0.5 * (T_xy_b + T_yx_b)
    }
    
    p_sym1 <- (1 + sum(abs(T_sym1_perm - mean(T_sym1_perm)) >= abs(T_sym1_obs - mean(T_sym1_perm)))) / (B_perm + 1)
    p_sym2 <- (1 + sum(abs(T_sym2_perm - mean(T_sym2_perm)) >= abs(T_sym2_obs - mean(T_sym2_perm)))) / (B_perm + 1)
    
    beta_AW <- beta_AW + as.numeric(p_aw < alpha)
    beta_Sym1 <- beta_Sym1 + as.numeric(p_sym1 < alpha)
    beta_Sym2 <- beta_Sym2 + as.numeric(p_sym2 < alpha)
  }
  
  fig15_left$AW[k] <- beta_AW / MC
  fig15_left$Sym1[k] <- beta_Sym1 / MC
  fig15_left$Sym2[k] <- beta_Sym2 / MC
  
  print("fig15_left")
  print(fig15_left[k, ])
  print("--------------------")
}

################################################################################
# Right panel: intersecting spheres
# X ~ U(S^2(c,1)), Y ~ U(S^2(-c,1)), c = (0.1, 0, 0)
################################################################################

fig15_right <- data.frame(
  n = ns,
  AW = NA_real_,
  Sym1 = NA_real_,
  Sym2 = NA_real_
)

for (k in seq_along(ns)) {
  n <- ns[k]
  m <- n
  
  indices_aw_xy <- generate_quadruples(n, m)
  indices_aw_yx <- generate_quadruples(m, n)
  
  pairs_X_excl <- vector("list", n)
  for (i5 in 1:n) {
    idx <- setdiff(1:n, i5)
    tmp <- as.matrix(expand.grid(i1 = idx, i2 = idx))
    pairs_X_excl[[i5]] <- tmp[tmp[, 1] != tmp[, 2], , drop = FALSE]
  }
  
  pairs_Y_excl <- vector("list", m)
  for (i6 in 1:m) {
    idx <- setdiff(1:m, i6)
    tmp <- as.matrix(expand.grid(i3 = idx, i4 = idx))
    pairs_Y_excl[[i6]] <- tmp[tmp[, 1] != tmp[, 2], , drop = FALSE]
  }
  
  denom_sym1 <- n * (n - 1) * (n - 2) * m * (m - 1) * (m - 2)
  
  beta_AW <- 0
  beta_Sym1 <- 0
  beta_Sym2 <- 0
  
  for (i in 1:MC) {
    X <- generate_sphere_samples(n, 3, radius = 1, center = c(0.1, 0, 0))
    Y <- generate_sphere_samples(m, 3, radius = 1, center = c(-0.1, 0, 0))
    
    ###########################################################################
    # AW (original test, asymptotic p-value)
    ###########################################################################
    
    dist_mat_XX <- as.matrix(dist(X))
    dist_mat_all_xy <- as.matrix(dist(rbind(X, Y)))
    
    T_aw_xy <- outer(
      dist_mat_XX,
      dist_mat_all_xy,
      FUN = function(x, y) as.numeric(x <= y)
    )
    
    W_aw <- (sum(T_aw_xy[indices_aw_xy]) - 0.5 * length(T_aw_xy[indices_aw_xy])) /
      ((m + n - 6) * (n - 3) * (n - 4) * length(T_aw_xy[indices_aw_xy]) * (2 * P1_unif - 1) / 8)^(1/2)
    
    p_aw <- 2 * pnorm(-abs(W_aw))
    
    ###########################################################################
    # Sym1 observed statistic
    ###########################################################################
    
    dist_mat_YY <- as.matrix(dist(Y))
    dist_mat_XY <- dist_mat_all_xy[1:n, (n + 1):(n + m), drop = FALSE]
    
    dxx_excl <- vector("list", n)
    for (i5 in 1:n) {
      dxx_excl[[i5]] <- dist_mat_XX[pairs_X_excl[[i5]]]
    }
    
    dyy_excl_sorted <- vector("list", m)
    for (i6 in 1:m) {
      dyy_excl_sorted[[i6]] <- sort(dist_mat_YY[pairs_Y_excl[[i6]]], method = "quick")
    }
    
    T_sym1_obs <- 0
    for (i5 in 1:n) {
      for (i6 in 1:m) {
        thr <- 2 * dist_mat_XY[i5, i6]
        T_sym1_obs <- T_sym1_obs +
          sum(findInterval(thr - dxx_excl[[i5]], dyy_excl_sorted[[i6]], rightmost.closed = TRUE))
      }
    }
    T_sym1_obs <- T_sym1_obs / denom_sym1
    
    ###########################################################################
    # Sym2 observed statistic
    ###########################################################################
    
    dist_mat_all_yx <- as.matrix(dist(rbind(Y, X)))
    
    T_aw_yx <- outer(
      dist_mat_YY,
      dist_mat_all_yx,
      FUN = function(x, y) as.numeric(x <= y)
    )
    
    T_xy_obs <- sum(T_aw_xy[indices_aw_xy]) / length(T_aw_xy[indices_aw_xy])
    T_yx_obs <- sum(T_aw_yx[indices_aw_yx]) / length(T_aw_yx[indices_aw_yx])
    T_sym2_obs <- 0.5 * (T_xy_obs + T_yx_obs)
    
    ###########################################################################
    # Permutation p-values for Sym1 and Sym2
    ###########################################################################
    
    Z <- rbind(X, Y)
    T_sym1_perm <- numeric(B_perm)
    T_sym2_perm <- numeric(B_perm)
    
    for (b in 1:B_perm) {
      id <- sample.int(n + m, n + m, replace = FALSE)
      Xb <- Z[id[1:n], , drop = FALSE]
      Yb <- Z[id[(n + 1):(n + m)], , drop = FALSE]
      
      dist_mat_XX_b <- as.matrix(dist(Xb))
      dist_mat_YY_b <- as.matrix(dist(Yb))
      
      dist_mat_all_xy_b <- as.matrix(dist(rbind(Xb, Yb)))
      dist_mat_XY_b <- dist_mat_all_xy_b[1:n, (n + 1):(n + m), drop = FALSE]
      
      # Sym1 permutation statistic
      dxx_excl_b <- vector("list", n)
      for (i5 in 1:n) {
        dxx_excl_b[[i5]] <- dist_mat_XX_b[pairs_X_excl[[i5]]]
      }
      
      dyy_excl_sorted_b <- vector("list", m)
      for (i6 in 1:m) {
        dyy_excl_sorted_b[[i6]] <- sort(dist_mat_YY_b[pairs_Y_excl[[i6]]], method = "quick")
      }
      
      T_sym1_b <- 0
      for (i5 in 1:n) {
        for (i6 in 1:m) {
          thr_b <- 2 * dist_mat_XY_b[i5, i6]
          T_sym1_b <- T_sym1_b +
            sum(findInterval(thr_b - dxx_excl_b[[i5]], dyy_excl_sorted_b[[i6]], rightmost.closed = TRUE))
        }
      }
      T_sym1_perm[b] <- T_sym1_b / denom_sym1
      
      # Sym2 permutation statistic
      T_aw_xy_b <- outer(
        dist_mat_XX_b,
        dist_mat_all_xy_b,
        FUN = function(x, y) as.numeric(x <= y)
      )
      
      dist_mat_all_yx_b <- as.matrix(dist(rbind(Yb, Xb)))
      T_aw_yx_b <- outer(
        dist_mat_YY_b,
        dist_mat_all_yx_b,
        FUN = function(x, y) as.numeric(x <= y)
      )
      
      T_xy_b <- sum(T_aw_xy_b[indices_aw_xy]) / length(T_aw_xy_b[indices_aw_xy])
      T_yx_b <- sum(T_aw_yx_b[indices_aw_yx]) / length(T_aw_yx_b[indices_aw_yx])
      T_sym2_perm[b] <- 0.5 * (T_xy_b + T_yx_b)
    }
    
    p_sym1 <- (1 + sum(abs(T_sym1_perm - mean(T_sym1_perm)) >= abs(T_sym1_obs - mean(T_sym1_perm)))) / (B_perm + 1)
    p_sym2 <- (1 + sum(abs(T_sym2_perm - mean(T_sym2_perm)) >= abs(T_sym2_obs - mean(T_sym2_perm)))) / (B_perm + 1)
    
    beta_AW <- beta_AW + as.numeric(p_aw < alpha)
    beta_Sym1 <- beta_Sym1 + as.numeric(p_sym1 < alpha)
    beta_Sym2 <- beta_Sym2 + as.numeric(p_sym2 < alpha)
  }
  
  fig15_right$AW[k] <- beta_AW / MC
  fig15_right$Sym1[k] <- beta_Sym1 / MC
  fig15_right$Sym2[k] <- beta_Sym2 / MC
  
  print("fig15_right")
  print(fig15_right[k, ])
  print("--------------------")
}

fig15_results <- list(
  concentric_spheres = fig15_left,
  intersecting_spheres = fig15_right
)

fig15_results

#################################################################################
################################################################################

fig15_left_long <- fig15_left %>%
  rename(`P1=0.52` = AW) %>%
  pivot_longer(cols = -n, names_to = "Test", values_to = "Power")

fig15_right_long <- fig15_right %>%
  rename(`P1=0.52` = AW) %>%
  pivot_longer(cols = -n, names_to = "Test", values_to = "Power")

p15_left <- ggplot(fig15_left_long, aes(x = n, y = Power, color = Test, linetype = Test)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("P1=0.52" = "red", "Sym1" = "green", "Sym2" = "blue")) +
  scale_linetype_manual(values = c("P1=0.52" = "solid", "Sym1" = "solid", "Sym2" = "solid")) +
  labs(
    title = "Power: uniform on concentric spheres",
    x = "Sample size",
    y = "Power",
    color = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(breaks = fig15_left$n) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

p15_right <- ggplot(fig15_right_long, aes(x = n, y = Power, color = Test, linetype = Test)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("P1=0.52" = "red", "Sym1" = "green", "Sym2" = "blue")) +
  scale_linetype_manual(values = c("P1=0.52" = "solid", "Sym1" = "solid", "Sym2" = "solid")) +
  labs(
    title = "Power: uniform on a sphere, shifted spheres",
    x = "Sample size",
    y = "Power",
    color = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(breaks = fig15_right$n) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

p15_left + p15_right