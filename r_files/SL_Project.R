### 0. Installing packages ###
library(readxl)
library(ggplot2)
library(mdatools)

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

## Basic statistics 

head(metabolomics)
str(metabolomics)
table(metabolomics$sex)


X <- metabolomics[, colnames(metabolomics) != "sex"]

## 1.2 ___ ##

# Checking the data for missing values and zero/negatives (for the log)
sum(is.na(X))  
sum(X <= 0, na.rm = TRUE) # There are two 0 values         
dim(X)

# Note: 61 Female / 60 Male samples (121 total). The paper's text says 60 + 60,
# but its own Fig. 1A also shows 121 samples (e.g. X132_121), so all are kept.

# Check the values of compound 109708 and creatinine
idx <- which(X <= 0, arr.ind = TRUE)
data.frame(sample = rownames(X)[idx[, 1]],
           metabolite = colnames(X)[idx[, 2]],
           value = X[idx])
# These two values are the same ones found in the paper (page. 4), 
# which they replaced it by 1.

# So following the paper we replace those values by 1
X_fix <- X
X_fix[X_fix == 0] <- 1
Xlog <- log10(X_fix)

# Summary
summ <- data.frame(
  mean   = sapply(X, mean),
  sd     = sapply(X, sd),
  median = sapply(X, median),
  skew   = sapply(X, function(x) mean(scale(x)^3))
)

head(summ) # they are very skewed, another reason to do the log transformation


# Gaussian distribution with Shapiro-Wilk test
shapiro_p <- function(M) sapply(M, function(x) shapiro.test(x)$p.value)
p_raw <- shapiro_p(X_fix)
p_log <- shapiro_p(Xlog)

c(raw_gaussian = mean(p_raw > 0.05),     # paper states --> ~23% Gaussian (77% non-normal)
  log_gaussian = mean(p_log > 0.05))

# Here we obtain raw(12% gaussian) and log(22.2% gaussian)
# 77.8% of metabolites fail the test
# Close to the 77% they report

skew_raw <- sapply(X_fix, function(x) mean(scale(x)^3))
skew_log <- sapply(Xlog,  function(x) mean(scale(x)^3))
c(median_abs_skew_raw = median(abs(skew_raw)),
  median_abs_skew_log = median(abs(skew_log))) # We reduced the skew median from ~2 to ~0.69

# Some Q-Q plot to see comparison
png("qq_plots.png", width = 1600, height = 1800, res = 200)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (m in c("xylitol", "glucose", "xylose")) {
  qqnorm(X_fix[[m]], main = paste(m, "raw"));   qqline(X_fix[[m]])
  qqnorm(Xlog[[m]],  main = paste(m, "log10")); qqline(Xlog[[m]])
}
dev.off()

# ---------------- PCA ----------------------

# Samples with the most extreme values (|z| > 3 on log data)
Z <- scale(Xlog)
sort(rowSums(abs(Z) > 3), decreasing = TRUE)[1:10]


pca_base <- prcomp(Xlog, center = TRUE, scale. = TRUE)
var_exp  <- pca_base$sdev^2 / sum(pca_base$sdev^2)
round(100 * var_exp[1:5], 2)
# PC1 = 19.7%, PC2 = 5.9% (paper: 15.4% and 5.5%)

# PCA model
m <- pca(as.matrix(Xlog), ncomp = 8, center = TRUE, scale = TRUE)
                   
# Scree plot: variance per component and cumulative variance
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plotVariance(m, type = "h", show.labels = TRUE, labels = "values",
             main = "Variance Explained")
plotCumVariance(m, type = "b", show.labels = TRUE,
                main = "Cumulative Variance")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cum <- cumsum(var_exp)
which(cum >= 0.50)[1]     # PCs for 50% of the variance
which(cum >= 0.80)[1]     # PCs for 80% of the variance
round(100 * cum[c(2, 5, 10)], 1)
# 10 PCs are needed for 50% of the variance and 40 PCs for 80%
# So it is not low-dimensional: after PC1 the variance decays slowly. 
# A 2D score plot therefore shows only about a quarter of the total variance, and a lack of
# sex separation there does not rule out separation in other directions.
# So the supervised analysis (PLS-DA) could be used.

# Score plot
plotScores(m, comp = c(1, 2), show.labels = FALSE,
           main = "Scores PC1 vs PC2")
# Score plot coloured by sex
plotScores(m, comp = c(1, 2), cgroup = sex, show.labels = FALSE,
           main = "Scores coloured by sex")

# With labels
plotScores(m, comp = c(1, 2), cgroup = sex, show.labels = TRUE,
           main = "Scores with sample labels")
                   
# Residuals
plotResiduals(m, ncomp = 5, show.labels = FALSE,
              main = "Residuals (5 PCs)")
                   

# Basic statistics
# Gaussianity?
#   Does it improve with transformation
# Outliers? FALTA AIXÔ CREC
# PCA
#   Score-plot
#   Scree-plot ?
#   How many dimensions retain the most variance
