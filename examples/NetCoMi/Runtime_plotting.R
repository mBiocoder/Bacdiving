#Create random SILVA_input.txt file with 10000 rows
silva_ids <- c("AB681649",
               "AB121974",
               "EU847536",
               "D30768",
               "L35516",
               "AB681086",
               "AB052706",
               "AF144407",
               "AF363064",
               "AJ430586")

df <- data.frame(silva_ids) 
df <- df[rep(seq_len(nrow(df)), each = 1000), ]

#write as .txt file
write.table(df, "SILVA-input.txt", append = FALSE, sep = "\n",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


################################################################################

#Plot BacDiving runtime
library(ggplot2)

# Create dataframe with measured runtime data
runtime <- c(627.73, 405.99, 309.49, 385.75, 42.02, 116.66, 179.17, 236.76, 583.73, 146.39, 102.34, 4665.44)
ASV_number <- c(16797, 6863, 4047, 5534, 718, 1561, 1091, 22206, 14812, 1051, 944, 10000)
rows_after_NA_removal <- c(853, 454 , 381 , 571 , 50 , 256 , 171 , 327 , 836, 192 , 147 , 10000)
dataset_name <- c("Agp", "Fukui", "Hugerth", "Liu", "Lopresti", "Mars", "Nagel", "Pozuelo", "Zeber", "Zhu", "Zhuang", "SILVA-input")
data=data.frame(runtime, ASV_number, dataset_name, rows_after_NA_removal)

#Plot datasets
plot(rows_after_NA_removal, runtime, pch = 16, cex = 1.3, col = "blue", main = "Empirical runtime anaylsis of the function bacdive_call() \n using linear regression", xlab = "Input size after removing unknown species", ylab = "Runtime in sec")

#linear regression
lm(runtime ~ rows_after_NA_removal)
abline(lm(runtime ~ rows_after_NA_removal))
