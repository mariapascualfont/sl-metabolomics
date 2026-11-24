##############################################################################
# SL_Project.R
# Metabolomics data: pre-processing, normality assessment and PCA
#
# Project structure:
#   Project/
#     |-- data/     original xlsx file
#     |-- r_files/  this script (set it as working directory when running it)
#     |-- results/  saved objects (.rds) and tables (.csv)
#     |-- rmd/      R Markdown report
#     |-- images/   figures (vector format: SVG)
##############################################################################


### 0. Installing packages and project setup ###

library(readxl)
library(ggplot2)
library(mdatools)

## Output folders (relative to r_files/)
img_dir <- "../images"
res_dir <- "../results"
dir.create(img_dir, showWarnings = FALSE)
dir.create(res_dir, showWarnings = FALSE)

## Helper: draws a plot straight into an SVG file (vector format, so it keeps
## the best possible resolution at any zoom level)
save_plot <- function(file, plot_expr, width = 8, height = 6) {
  svg(file, width = width, height = height)
  on.exit(dev.off())
  plot_expr
  invisible(file)
}


### 1. Data Pre-Processing ###

## 1.1 Data Frame Creation ##

## Read table and delete columns from previous statistical analysis
df <- read_excel("../data/41598_2018_29592_MOESM2_ESM.xlsx")
df_clean <- df[, 15:ncol(df), drop = FALSE]
df_clean <- as.data.frame(df_clean)

row.names(df_clean) <- df_clean[[1]]
df_clean <- df_clean[, -1, drop = FALSE]

## Retrieve relevant metadata information (sex and label)
sample_labels <- as.character(unlist(df_clean["label", ]))
sex_chr <- as.character(unlist(df_clean["comment", ]))


## Delete rows from unnecessary metadata information
metabolomics <- df_clean[!rownames(df_clean) %in% c("mx class","mx sample", "comment",
                                                    "species","organ","treatment","label"), ]

## Transpose, rename rows and convert measures to numeric
metabolomics <- t(metabolomics)
metabolomics <- as.data.frame(metabolomics)
rownames(metabolomics) <- sample_labels
metabolomics[] <- lapply(metabolomics, function(x) {
  as.numeric(as.character(x))
})

## Add outcome variable (sex)
metabolomics$sex <- factor(sex_chr)
sex <- metabolomics$sex

## Save the clean data frame
saveRDS(metabolomics, file.path(res_dir, "metabolomics_clean_df.rds"))


## 1.2 Basic Exploration ##

# Only the first columns are shown, the table is very wide
head(metabolomics[, c(1:5, ncol(metabolomics))])
str(metabolomics, list.len = 5)
table(metabolomics$sex)

# Matrix with the metabolite measures only (outcome removed)
X <- metabolomics[, colnames(metabolomics) != "sex"]


## 1.3 Missing Values, Zeros and Log Transformation ##

# Checking the data for missing values and zero/negatives (for the log)
sum(is.na(X))  
sum(X <= 0, na.rm = TRUE) # There are two 0 values         
dim(X)

# Note: 61 Female / 60 Male samples (121 total). The paper's text says 60 + 60,
# but its own Fig. 1A also shows 121 samples (e.g. X132_121), so all are kept.

# Check the values of compound 109708 and creatinine
idx <- which(X <= 0, arr.ind = TRUE)
zero_values <- data.frame(sample = rownames(X)[idx[, 1]],
                          metabolite = colnames(X)[idx[, 2]],
                          value = X[idx])
zero_values
# These two values are the same ones found in the paper (page. 4), 
# which they replaced it by 1.

# So following the paper we replace those values by 1
X_fix <- X
X_fix[X_fix == 0] <- 1
Xlog <- log10(X_fix)

## Save the corrected and log-transformed matrices
saveRDS(X_fix, file.path(res_dir, "X_fix_zeros_replaced.rds"))
saveRDS(Xlog,  file.path(res_dir, "Xlog_log10_transformed.rds"))
write.csv(zero_values, file.path(res_dir, "zero_values_replaced.csv"), row.names = FALSE)


## 1.4 Descriptive Summary ##

# Summary
summ <- data.frame(
  mean   = sapply(X, mean),
  sd     = sapply(X, sd),
  median = sapply(X, median),
  skew   = sapply(X, function(x) mean(scale(x)^3))
)

head(summ) # they are very skewed, another reason to do the log transformation
write.csv(summ, file.path(res_dir, "descriptive_summary_raw.csv"))


## 1.5 Gaussianity Assessment ##

# Gaussian distribution with Shapiro-Wilk test
shapiro_p <- function(M) sapply(M, function(x) shapiro.test(x)$p.value)
p_raw <- shapiro_p(X_fix)
p_log <- shapiro_p(Xlog)

c(raw_gaussian = mean(p_raw > 0.05),     # paper states --> ~23% Gaussian (77% non-normal)
  log_gaussian = mean(p_log > 0.05))

# Here we obtain raw(12% gaussian) and log(22.2% gaussian)
# 77.8% of metabolites fail the test
# Close to the 77% they report

# Skewness before and after the transformation
skew_raw <- sapply(X_fix, function(x) mean(scale(x)^3))
skew_log <- sapply(Xlog,  function(x) mean(scale(x)^3))
c(median_abs_skew_raw = median(abs(skew_raw)),
  median_abs_skew_log = median(abs(skew_log))) # We reduced the skew median from ~2 to ~0.69

# Store the per-metabolite normality results
normality <- data.frame(shapiro_p_raw = p_raw, shapiro_p_log = p_log,
                        skew_raw = skew_raw, skew_log = skew_log)
write.csv(normality, file.path(res_dir, "normality_shapiro_skewness.csv"))


## 1.6 Q-Q Plots ##

# Some Q-Q plot to see comparison (raw vs log10)
save_plot(file.path(img_dir, "qq_plots_raw_vs_log10_xylitol_glucose_xylose.svg"), {
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  for (met in c("xylitol", "glucose", "xylose")) {
    qqnorm(X_fix[[met]], main = paste(met, "raw"));   qqline(X_fix[[met]])
    qqnorm(Xlog[[met]],  main = paste(met, "log10")); qqline(Xlog[[met]])
  }
}, width = 8, height = 9)


### 2. Principal Component Analysis (PCA) ###

## 2.1 Extreme Values per Sample ##

# Samples with the most extreme values (|z| > 3 on log data)
Z <- scale(Xlog)
extreme_counts <- sort(rowSums(abs(Z) > 3), decreasing = TRUE)
extreme_counts[1:10]
write.csv(data.frame(sample = names(extreme_counts), n_extreme = extreme_counts),
          file.path(res_dir, "extreme_values_per_sample.csv"), row.names = FALSE)


## 2.2 PCA with prcomp ##

pca_base <- prcomp(Xlog, center = TRUE, scale. = TRUE)
var_exp  <- pca_base$sdev^2 / sum(pca_base$sdev^2)
round(100 * var_exp[1:5], 2)
# PC1 = 19.7%, PC2 = 5.9% (paper: 15.4% and 5.5%)
saveRDS(pca_base, file.path(res_dir, "pca_prcomp.rds"))


## 2.3 PCA Model (mdatools) ##

# PCA model
m <- pca(as.matrix(Xlog), ncomp = 8, center = TRUE, scale = TRUE)
saveRDS(m, file.path(res_dir, "pca_mdatools_model.rds"))


## 2.4 Variance Explained (Scree Plot) ##

# Scree plot: variance per component and cumulative variance
save_plot(file.path(img_dir, "pca_scree_variance_and_cumulative_variance.svg"), {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plotVariance(m, type = "h", show.labels = TRUE, labels = "values",
               main = "Variance Explained")
  plotCumVariance(m, type = "b", show.labels = TRUE,
                  main = "Cumulative Variance")
}, width = 12, height = 5)

cum <- cumsum(var_exp)
which(cum >= 0.50)[1]     # PCs for 50% of the variance
which(cum >= 0.80)[1]     # PCs for 80% of the variance
round(100 * cum[c(2, 5, 10)], 1)
# 10 PCs are needed for 50% of the variance and 40 PCs for 80%
# So it is not low-dimensional: after PC1 the variance decays slowly. 
# A 2D score plot therefore shows only about a quarter of the total variance, and a lack of
# sex separation there does not rule out separation in other directions.
# So the supervised analysis (PLS-DA) could be used.

# Variance table of all components
var_table <- data.frame(PC = seq_along(var_exp),
                        variance_pct   = 100 * var_exp,
                        cumulative_pct = 100 * cum)
write.csv(var_table, file.path(res_dir, "pca_variance_table.csv"), row.names = FALSE)


## 2.5 Score Plots ##

# Score plot
save_plot(file.path(img_dir, "pca_scores_pc1_vs_pc2.svg"), {
  plotScores(m, comp = c(1, 2), show.labels = FALSE,
             main = "Scores PC1 vs PC2")
})

# Score plot coloured by sex
save_plot(file.path(img_dir, "pca_scores_pc1_vs_pc2_coloured_by_sex.svg"), {
  plotScores(m, comp = c(1, 2), cgroup = sex, show.labels = FALSE,
             main = "Scores coloured by sex")
})

# With labels
save_plot(file.path(img_dir, "pca_scores_pc1_vs_pc2_with_sample_labels.svg"), {
  plotScores(m, comp = c(1, 2), cgroup = sex, show.labels = TRUE,
             main = "Scores with sample labels")
}, width = 10, height = 8)


## 2.6 Residuals ##

save_plot(file.path(img_dir, "pca_residuals_5_components.svg"), {
  plotResiduals(m, ncomp = 5, show.labels = FALSE,
                main = "Residuals (5 PCs)")
})


# Basic statistics
# Gaussianity?
#   Does it improve with transformation
# Outliers? FALTA AIXÔ CREC
# PCA
#   Score-plot
#   Scree-plot ?
#   How many dimensions retain the most variance
