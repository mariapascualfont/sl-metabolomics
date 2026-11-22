library(readxl)
library(ggplot2)


# Read table and delete columns from previous statistical analysis

df <- read_excel("~/Desktop/Bioinformatics/S5/SL/project/41598_2018_29592_MOESM2_ESM.xlsx")
df_clean <- df[, 15:ncol(df), drop = FALSE]
df_clean <- as.data.frame(df_clean)
row.names(df_clean) <- df_clean[[1]]
df_clean <- df_clean[, -1, drop = FALSE]

# Retrieve relevant metadata information (sex and label)
sample_labels <- as.character(unlist(df_clean["label", ]))
sex <- as.character(unlist(df_clean["comment", ]))


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
metabolomics$sex <- factor(sex)

# Basic statistics 

head(metabolomics)
str(metabolomics)
table(metabolomics$sex)

X   <- metabolomics[, colnames(metabolomics) != "sex"]
sex <- metabolomics$sex

# Do some check in the data for missing values and zero/negatives (for the log)
sum(is.na(X))  
sum(X <= 0, na.rm = TRUE) # There are 2 negative values         
dim(X)

# Check the negative values of compound 109708 and creatinine,
idx <- which(X <= 0, arr.ind = TRUE)
data.frame(sample = rownames(X)[idx[, 1]],
           metabolite = colnames(X)[idx[, 2]],
           value = X[idx])
# This two values are the same ones found in the paper (page. 4), 
# which they replaced it by 1.



# Basic statistics
# Gaussianity?
#   Does it improve with transformation
# Outliers?
# PCA
#   Score-plot
#   Scree-plot ?
#   How many dimensions retain the most variance
