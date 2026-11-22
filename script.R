library(readxl)

df <- read_excel("~/Desktop/Bioinformatics/S5/SL/project/41598_2018_29592_MOESM2_ESM.xlsx")

#from column 15 onwards

df_clean <- df[, 15:ncol(df), drop = FALSE]

df_clean <- as.data.frame(df_clean)
row.names(df_clean) <- df_clean[[1]]
df_clean <- df_clean[, -1, drop = FALSE]