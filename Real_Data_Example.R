###############################################################################
# Script: Real_Data_Example.R
#
# Purpose:
# Reproduce the example from the article "A TWO-SAMPLE TEST BASED ON AVERAGED 
# WILCOXON RANK SUMS OVER INTERPOINT DISTANCES" (https://arxiv.org/pdf/2408.10570) 
# using the GSE71661 miRNA count dataset, including:
#
#   1) differential expression analysis with edgeR,
#   2) stepwise feature selection based on the top-ranked miRNAs,
#   3) averaged interpoint Wilcoxon tests for dXX vs dXY and dYY vs dXY,
#   4) summary table of test statistics and p-values,
#   5) histogram plots of the corresponding distance distributions.
#
# Data:
# - Input file: GSE71661_processed_data.xlsx
# - Source: GEO accession GSE71661
# - Biological setting: pooled lung cancer vs pooled healthy/benign samples
#
# Analysis overview:
# - Raw miRNA counts are imported from Excel.
# - Samples are pooled into two classes.
# - Differential expression is performed with edgeR.
# - The top miRNAs are added stepwise by dimension.
# - Euclidean-distance averaged Wilcoxon tests are computed for:
#     H0,1: dXX = dXY
#     H0,2: dYY = dXY
#
# Required packages:
# - readxl
# - edgeR
# - ggplot2
# - gridExtra
#
# Important notes:
# - This script assumes that the custom function `AW_interpoint_test()`
#   is available in the R session or sourced from a separate file.
# - Internet access is required only if the dataset is downloaded directly
#   from GEO.
#
# Outputs:
# - Top tags of differential expression analysis, log-fold changes, p-values
# - Printed summary table with W statistics and p-values
# - Histogram panel for dXX, dXY, and dYY
#
# Author: Aljosa Marjanovic
# Date: 10.03.2026
# Project: Two sample averaged Wilcoxon over interpoint distances
###############################################################################

source("Utilities.R")


###########################################################################
# Data download + import + formatting
###########################################################################

# GEO processed data file used for the GSE71661  example
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71661&format=file&file=GSE71661%5Fprocessed%5Fdata%2Exlsx"
destfile <- "GSE71661_processed_data.xlsx"

if (!file.exists(destfile)) {
  message("Downloading input file: ", destfile)
  download.file(
    url = url,
    destfile = destfile,
    mode = "wb",
    method = "libcurl"
  )
} else {
  message("Using existing local file: ", destfile)
}

raw <- read_excel(destfile, sheet = 1)

miRNA_names <- raw[, 1]
Data <- raw[, -1]
Data <- as.data.frame(Data)
rownames(Data) <- as.vector(t(miRNA_names))

###########################################################################
# Pooled group labels
###########################################################################

Groups <- c(rep(1, times = 9), rep(2, times = 10),
            rep(1, times = 9), rep(2, times = 10),
            rep(1, times = 9), rep(2, times = 10))

Data <- rbind(Groups, Data)
rownames(Data) <- c("Groups", as.vector(t(miRNA_names)))

###########################################################################
# Differential Expression Analysis
###########################################################################

sub <- Data[-1, ]
rownames(sub) <- as.vector(t(miRNA_names))

sub <- as.matrix(sub)
mode(sub) <- "numeric"

data <- DGEList(counts = sub, group = factor(Groups))

filter <- rowSums(cpm(data) > 1) >= 1
d <- data[filter, , keep.lib.sizes = FALSE]

d <- calcNormFactors(d)

design.mat <- model.matrix(~ Groups)

d <- estimateDisp(d, design.mat)
test <- exactTest(d)

tT <- edgeR::topTags(test, n = 15, adjust.method = "BH")
tT

select <- rownames(tT$table)  # most differentially expressed miRNAs

###########################################################################
# Averaged interpoint Wilcoxon
###########################################################################


Table2 <- data.frame(
  Dimensions = paste0("1..", 2:10),
  W1 = NA_real_,
  p1 = NA_real_,
  W2 = NA_real_,
  p2 = NA_real_
)

###########################################################################
# XX vs XY
###########################################################################

Data1 <- cpm(as.matrix(Data[-1, ]))
rownames(Data1) <- as.vector(t(miRNA_names))

# AW_interpoint_test() is defined in Utilities.R and computes the
# averaged interpoint Wilcoxon test statistic for within-group vs between-group
# distance comparisons.

for (i in 2:10) {
  
  sel <- select[1:i]
  Data_sub <- Data1[sel, , drop = FALSE]
  
  X <- Data_sub[, 1:27, drop = FALSE]
  Y <- Data_sub[, 28:57, drop = FALSE]
  
  res_xy <- AW_interpoint_test(X, Y, method = "euclidean", P1 = 0.54)
  
  Table2$W1[i - 1] <- res_xy$statistic
  Table2$p1[i - 1] <- res_xy$p.value
}

###########################################################################
# YY vs XY
###########################################################################

Data_reordered <- Data
groups <- as.numeric(as.matrix(Data_reordered["Groups", ]))
column_order <- order(groups, decreasing = TRUE)
Data_reordered <- Data_reordered[, column_order]
rownames(Data_reordered) <- c("Groups", as.vector(t(miRNA_names)))

G <- factor(as.numeric(as.matrix(Data_reordered["Groups", ])))

Data2 <- cpm(as.matrix(Data_reordered[-1, ]))
rownames(Data2) <- as.vector(t(miRNA_names))

for (i in 2:10) {
  
  sel <- select[1:i]
  Data_sub <- Data2[sel, , drop = FALSE]
  
  X <- Data_sub[, G == 1, drop = FALSE]
  Y <- Data_sub[, G == 2, drop = FALSE]
  
  res_yx <- AW_interpoint_test(Y, X, method = "euclidean", P1 = 0.54)
  
  Table2$W2[i - 1] <- res_yx$statistic
  Table2$p2[i - 1] <- res_yx$p.value
}

###########################################################################
# Output
###########################################################################

Table2_print <- Table2

Table2_print$W1 <- sprintf("%.6f", Table2_print$W1)
Table2_print$p1 <- formatC(Table2_print$p1, format = "f", digits = 5)
Table2_print$W2 <- sprintf("%.6f", Table2_print$W2)
Table2_print$p2 <- formatC(Table2_print$p2, format = "f", digits = 5)

print(Table2_print, row.names = FALSE)


#########################################################################
# Histogram plots
##########################################################################

sel_hist <- select[1:6]

Data_hist <- cpm(as.matrix(Data[-1, ]))
rownames(Data_hist) <- as.vector(t(miRNA_names))
Data_hist <- Data_hist[sel_hist, , drop = FALSE]

# According to the pooled labels in the script:
# group 1 = 27 samples, group 2 = 30 samples
# and in the paper X = healthy, Y = sick
X_hist <- Data_hist[, Groups == 1, drop = FALSE]  # healthy, 27
Y_hist <- Data_hist[, Groups == 2, drop = FALSE]  # sick, 30

# same construction as in AW_interpoint_test(): dist(t(cbind(X, Y)))
Z_hist <- cbind(X_hist, Y_hist)
dist_mat <- as.matrix(dist(t(Z_hist), method = "euclidean"))

dxx <- unique(as.vector(dist_mat[1:27,1:27]))
dyy <- unique(as.vector(dist_mat[28:57,28:57]))
dxy <- unique(as.vector(dist_mat[1:27,28:57]))

p1 = ggplot(data.frame(distances = dxx), aes(x = distances, y = ..density..)) +
  geom_histogram(bins = 100, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "dXX",
       x = "Distance",
       y = "Density") +
  xlim(0, 8000) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)
  )


p2 = ggplot(data.frame(distances = dxy), aes(x = distances, y = ..density..)) +
  geom_histogram(bins = 100, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "dXY",
       x = "Distance",
       y = "Density") +
  xlim(0, 8000) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)
  )


p3 = ggplot(data.frame(distances = dyy), aes(x = distances, y = ..density..)) +
  geom_histogram(bins = 100, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "dYY",
       x = "Distance",
       y = "Density") +
  xlim(0, 8000) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 18)
  )



g = grid.arrange(p1, p2, p3, ncol=3)
