###############################################################################
# File: Utilities.R
#
# Purpose:
# Helper functions for interpoint-distance two-sample testing and simulation.
# This file provides utility functions used in the real-data example and
# in related simulation experiments.
#
# Contents:
# - generate_quadruples():
#     Construct index quadruples for the averaged interpoint Wilcoxon statistic.
# - AW_interpoint_test():
#     Compute the averaged interpoint Wilcoxon test statistic and p-value for
#     comparing within-group and between-group distance distributions.
# - ind_wilcox_reject():
#     Perform a Wilcoxon rank-sum comparison on independently sampled XX and XY
#     distances and return a binary rejection indicator.
# - generate_sphere_samples(), generate_uniform_sphere():
#     Generate points uniformly on a sphere in Euclidean space.
# - generate_ellipse_samples():
#     Generate points on a planar ellipse.
#
# Notes:
# - Input data are coerced to matrices where needed.
# - AW_interpoint_test() is used by Real_Data_Example.R.
# - This file is intended to be sourced before running the main analysis script:
#     source("Utilities.R")
#
# Outputs:
# - Test functions return either a list with statistic/p-value or a rejection
#   indicator.
# - Sample-generation functions return matrices of simulated coordinates.
#
# Author: Aljosa Marjanovic
# Date: 10.03.2026
# Project:  Two sample averaged Wilcoxon over interpoint distances
###############################################################################

library(tidyverse)
library(gridExtra)
library(readxl)
library(edgeR)
library(patchwork)
library(ICSNP)
library(cramer)

library(energy)      # eqdist.etest
library(SpatialNP)   # sr.loc.test
library(DepthProc)   # mWilcoxonTest
library(kernlab)     # rbfdot
library(Ecume)       # mmd_test
library(FNN)         # get.knn
library(npmv)

###############################################################################
# Index construction for interpoint distance comparisons
###############################################################################

generate_quadruples <- function(n, m) {
  combinations <- expand.grid(
    a = 1:(n - 1),
    b = 2:n,
    c = 1:n,
    d = (n + 1):(n + m)
  )
  
  filtered_combinations <- subset(combinations, a < b & c != a & c != b)
  
  cbind(filtered_combinations$a,
        filtered_combinations$b,
        filtered_combinations$c,
        filtered_combinations$d)
}

###############################################################################
# Averaged interpoint Wilcoxon test
###############################################################################

AW_interpoint_test <- function(X, Y, method = "euclidean", P1 = 2/3) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n <- ncol(X)
  m <- ncol(Y)
  
  Z <- cbind(X, Y)
  D <- as.matrix(dist(t(Z), method = method))
  
  indices <- generate_quadruples(n, m)
  
  T <- outer(D, D, FUN = function(x, y) as.numeric(x <= y))
  
  W <- (sum(T[indices]) - 0.5 * nrow(indices)) /
    sqrt((m + n - 6) * (n - 3) * (n - 4) * nrow(indices) * (2 * P1 - 1) / 8)
  
  p <- pnorm(abs(W), lower.tail = FALSE)
  
  list(statistic = W, p.value = p)
}

###############################################################################
# Wilcoxon test on independently sampled distances
# Compares an independent sample of XX distances vs an independent sample 
# of XY distances
###############################################################################

ind_wilcox_reject <- function(X, Y, method = "euclidean", alpha = 0.05) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n <- nrow(X)
  m <- nrow(Y)
  
  r <- min(floor(n / 3), m)
  if (r < 1) return(0)
  
  # choose 3r distinct X-observations:
  # 2r for r disjoint XX pairs, and r more for XY distances
  ix <- sample(seq_len(n), 3 * r, replace = FALSE)
  
  x1   <- ix[seq(1, 2 * r, by = 2)]
  x2   <- ix[seq(2, 2 * r, by = 2)]
  ix_b <- ix[(2 * r + 1):(3 * r)]
  
  # choose r distinct Y-observations
  iy_b <- sample(seq_len(m), r, replace = FALSE)
  
  # independent XX and XY distances
  dxx <- numeric(r)
  dxy <- numeric(r)
  
  for (i in seq_len(r)) {
    dxx[i] <- as.numeric(dist(rbind(X[x1[i], , drop = FALSE], X[x2[i], , drop = FALSE]),method = method))
    
    dxy[i] <- as.numeric(dist(rbind(X[ix_b[i], , drop = FALSE], Y[iy_b[i], , drop = FALSE]),method = method))
  }
  
  pval <- wilcox.test(dxx, dxy, alternative = "two.sided", exact = FALSE)$p.value
  as.numeric(pval < alpha)
}

###############################################################################
# Uniform sampling on a sphere
###############################################################################

generate_sphere_samples <- function(n, dim, radius = 1, center = NULL) {
  samples <- matrix(rnorm(n * dim), nrow = n, ncol = dim)
  samples <- samples / sqrt(rowSums(samples^2))
  samples <- radius * samples
  
  if (is.null(center)) {
    center <- rep(0, dim)
  }
  
  samples <- sweep(samples, 2, center, "+")
  return(samples)
}

###############################################################################
# Convenience wrapper for centered sphere sampling
###############################################################################

generate_uniform_sphere <- function(n, d, radius = 1) {
  generate_sphere_samples(n = n, dim = d, radius = radius, center = rep(0, d))
}

###############################################################################
# Sampling points on a planar ellipse
###############################################################################

generate_ellipse_samples <- function(n, a, e) {
  theta <- runif(n, 0, 2 * pi)
  b <- a * sqrt(1 - e^2)
  samples <- cbind(a * cos(theta), b * sin(theta))
  return(samples)
}


