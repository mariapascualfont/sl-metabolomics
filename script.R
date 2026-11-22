library(readxl)
library(ggplot2)


# Read table and delete columns from previous statistical analysis

df <- read_excel("C:/Users/mirei/Documents/Bioinf_2/3r_any/1r_quatri/Statistics/Project/41598_2018_29592_MOESM2_ESM.xlsx")
df_clean <- df[, 15:ncol(df), drop = FALSE]
df_clean <- as.data.frame(df_clean)
row.names(df_clean) <- df_clean[[1]]
df_clean <- df_clean[, -1, drop = FALSE]

# Retrieve relevant metadata information (sex and label)
sample_labels <- as.character(unlist(df_clean["label", ]))
sex_chr <- as.character(unlist(df_clean["comment", ]))


# Delete rows from unnecessary metadata information
metabolomics <- df_clean[!rownames(df_clean) %in% c("mx class","mx sample", "comment",
                                                    "species","organ","treatment","label"), ]

# Transpose, rename rows and convert measures to numeric
metabolomics <- t(metabolomics)
metabolomics <- as.data.frame(metabolomics)
rownames(metabolomics) <- sample_labels
metabolomics[] <- lapply(metabolomics, function(x) {
  as.numeric(as.character(x))
})

# Add outcome variable (sex)
metabolomics$sex <- factor(sex_chr)
sex <- metabolomics$sex

# Basic statistics 

head(metabolomics)
str(metabolomics)
table(metabolomics$sex)


X   <- metabolomics[, colnames(metabolomics) != "sex"]


# Do some check in the data for missing values and zero/negatives (for the log)
sum(is.na(X))  
sum(X <= 0, na.rm = TRUE) # There are 2 0 values         
dim(X)
# Note: 61 Female / 60 Male samples (121 total). The paper's text says 60 + 60,
# but its own Fig. 1A also shows 121 samples (e.g. X132_121), so all are kept.

# Check the values of compound 109708 and creatinine,
idx <- which(X <= 0, arr.ind = TRUE)
data.frame(sample = rownames(X)[idx[, 1]],
           metabolite = colnames(X)[idx[, 2]],
           value = X[idx])
# This two values are the same ones found in the paper (page. 4), 
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







# Basic statistics
# Gaussianity?
#   Does it improve with transformation
# Outliers?
# PCA
#   Score-plot
#   Scree-plot ?
#   How many dimensions retain the most variance
