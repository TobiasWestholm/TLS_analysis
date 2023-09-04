####################################################
########### TLS dataset analysis## #################
####################################################
library(visdat) # For visualizing data
library(ggplot2) # for nice plots
library(naniar) # For visualizing and working with missing data
library(GGally) # for combined plots
library(patchwork)
library(rstatix)
library(ggpubr)
library(EnhancedVolcano)

library(caret)
library(simputation) # For handling missing values
library(tidyverse) # Metapackage for dplyr and ggplot2
library("mice")



## Load dataset --------------------------------------------------------
# read.table is a function that allows to read multiple kinds of text files. It admits the following arguments:
df <- read.table(file = "TLS_COMBINED_ROI_measurements.csv", #Name of text file.
                 sep = "\t",                       #Separation character.
                 header = TRUE,                   #If column names are in the first row.
                 na.strings = "",               #Character to be marked as missing value.
                 stringsAsFactors = FALSE)         #żconvert string to factors

# data exploration
summary(df)
View(df)
head(df)

#-------------------------------------------------------------------------------------------------
#Pearson correlation pheatSets--------------------------------

corrSet <- df[,1:22]
# Replace dots with plus signs and underscores with spaces in the feature names
names(corrSet) <- gsub("\\.", "+", names(corrSet))
names(corrSet) <- gsub("_", " ", names(corrSet))

# combined spatial plot
ggpairs(corrSet[,c(3:8,15:18)], aes(alpha = 0.3))

# Mean regression
corrSet <- df[,1:22]
ggplot(data = corrSet, aes(x = CD45_mean_intensity, y = CD20_mean_intensity)) + 
  geom_point() +
  labs(title = "CD45-CD20 mean intensity correlation",
       y = "CD20 mean intensity",
       x = "CD45 mean intensity",) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  annotate("text", x = min(corrSet$CD45_mean_intensity), y = max(corrSet$CD20_mean_intensity),
           label = paste0("R-Square = ", round(summary(lm(CD20_mean_intensity ~ CD45_mean_intensity, data = corrSet))$r.squared, 3)), 
           size = 4, vjust = -4, hjust = 0) +
  annotate("text", x = min(corrSet$CD45_mean_intensity), y = max(corrSet$CD20_mean_intensity),
           label = paste0("p-value = ", round(summary(lm(CD20_mean_intensity ~ CD45_mean_intensity, data = corrSet))$coefficients[8], 5)), 
           size = 4, vjust = -2, hjust = 0)

#Median regression
ggplot(data = corrSet, aes(x = CD45_median_intensity, y = CD20_median_intensity)) + 
  geom_point() +
  labs(title = "CD45-CD20 median intensity correlation",
       y = "CD20 median intensity",
       x = "CD45 median intensity",) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  annotate("text", x = min(corrSet$CD45_median_intensity), y = max(corrSet$CD20_median_intensity),
           label = paste0("R-Square = ", round(summary(lm(CD20_median_intensity ~ CD45_median_intensity, data = corrSet))$r.squared, 3)), 
           size = 4, vjust = -3, hjust = 0) +
  annotate("text", x = min(corrSet$CD45_median_intensity), y = max(corrSet$CD20_median_intensity),
           label = paste0("p-value = ", round(summary(lm(CD20_median_intensity ~ CD45_median_intensity, data = corrSet))$coefficients[8], 5)), 
           size = 4, vjust = -1, hjust = 0)

#std regression
ggplot(data = corrSet, aes(x = CD45_intensity_standard_deviation, y = CD20_intensity_standard_deviation)) + 
  geom_point() +
  labs(title = "CD45-CD20 standard deviation correlation",
       y = "CD20 standard deviation intensity",
       x = "CD45 standard deviation intensity",) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  annotate("text", x = min(corrSet$CD45_intensity_standard_deviation), y = max(corrSet$CD20_intensity_standard_deviation),
           label = paste0("R-Square = ", round(summary(lm(CD20_intensity_standard_deviation ~ CD45_intensity_standard_deviation, data = corrSet))$r.squared, 3)), 
           size = 4, vjust = 0, hjust = 0) +
  annotate("text", x = min(corrSet$CD45_intensity_standard_deviation), y = max(corrSet$CD20_intensity_standard_deviation),
           label = paste0("p-value = ", round(summary(lm(CD20_intensity_standard_deviation ~ CD45_intensity_standard_deviation, data = corrSet))$coefficients[8], 5)), 
           size = 4, vjust = 2, hjust = 0)

#Cells per area regression
ggplot(data = corrSet, aes(x = CD45._cell_density, y = CD20._cell_density)) + 
  geom_point() +
  labs(title = expression("CD45+-CD20+ cell density correlation"),
       y = expression(paste("CD20+ cell density (", "µm" ^-2, ")")),
       x = expression(paste("CD45+ cell density (", "µm" ^-2, ")"))) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  annotate("text", x = min(corrSet$CD45._cell_density), y = max(corrSet$CD20._cell_density),
           label = paste0("R-Square = ", round(summary(lm(CD20._cell_density ~ CD45._cell_density, data = corrSet))$r.squared, 3)), 
           size = 4, vjust = -2, hjust = 0) +
  annotate("text", x = min(corrSet$CD45._cell_density), y = max(corrSet$CD20._cell_density),
           label = paste0("p-value = ", round(summary(lm(CD20._cell_density ~ CD45._cell_density, data = corrSet))$coefficients[8], 5)), 
           size = 4, vjust = 0, hjust = 0)

#Cell ratio regression
ggplot(data = corrSet, aes(x = CD45._ratio, y = CD20._ratio)) + 
  geom_point() +
  labs(title = "CD45+-CD20+ ratio correlation",
       y = "CD20+ ratio",
       x = "CD45+ ratio",) +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  annotate("text", x = min(corrSet$CD45._ratio), y = max(corrSet$CD20._ratio),
           label = paste0("R-Square = ", round(summary(lm(CD20._ratio ~ CD45._ratio, data = corrSet))$r.squared, 3)), 
           size = 4, vjust = 0, hjust = 0) +
  annotate("text", x = min(corrSet$CD45._ratio), y = max(corrSet$CD20._ratio),
           label = paste0("p-value = ", round(summary(lm(CD20._ratio ~ CD45._ratio, data = corrSet))$coefficients[8], 5)), 
           size = 4, vjust = 2, hjust = 0)







#-------------------------------------------------------------------------------------------------
#Pheatmap clustering--------------------------
library(pheatmap)

#fix the matrix for pheatmap (NU LOGARITMERAS ALL DATA ALLA NOLLVÄRDEN SÄTTS VÄLDIGT LÅGT)
pheatSet <- df[, -c(1, 2)] #alla variabler
#pheatSet <- df[, 3:22] #spatiala variabler
# Find rows with 0 values in "CD23_cell_density" and "CD23_ratio"
zero_indices_density <- which(pheatSet$CD23._cell_density == 0)
zero_indices_ratio <- which(pheatSet$CD23._ratio == 0)

# Replace 0 values with "CD23_mean_intensity" divided by max("CD23_mean_intensity")
pheatSet[zero_indices_density, "CD23._cell_density"] <- pheatSet[zero_indices_density, "CD23_mean_intensity"] / max(pheatSet$CD23_mean_intensity)*max(pheatSet$CD23._cell_density)
pheatSet[zero_indices_ratio, "CD23._ratio"] <- pheatSet[zero_indices_ratio, "CD23_mean_intensity"] / max(pheatSet$CD23_mean_intensity)*max(pheatSet$CD23._ratio)
pheatSet <- log(pheatSet)
#pheatSet <- df[, -c(1:23)] #barcodevariabler
rownames(pheatSet) <- paste0("ROI ", df$ROI)


library(data.table)
#transpose data frame
df_t <- transpose(pheatSet)

# Scale each row in the dataframe
df_t <- t(apply(df_t, 1, scale))

#redefine row and column names
rownames(df_t) <- gsub("\\.", "-", gsub("_", " ", colnames(pheatSet)))
#change the first 20 rows to + instead of -
rownames(df_t)[1:20] <- gsub("-", "+", rownames(df_t)[1:20])
colnames(df_t) <- gsub("\\.1", " (1)", gsub("\\.2", " (2)", rownames(pheatSet)))

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#plot a ROI-clustered heatmap
heatmap_result <- pheatmap(df_t, cluster_row = TRUE, cutree_cols = 4, 
                           clustering_callback = callback, angle_col = "90")
dendrogram_result <- heatmap_result$tree_col

# Plot the dendrogram to visualize column dissimilarity
plot(dendrogram_result)



#spatial analysis
pheatSet_spatial <- df_t[1:20,]
#plot a ROI-clustered heatmap
heatmap_result <- pheatmap(pheatSet_spatial, cluster_row = FALSE, cutree_cols = 2, clustering_callback = callback, angle_col = "90")
dendrogram_result <- heatmap_result$tree_col

# Plot the dendrogram to visualize column dissimilarity
plot(dendrogram_result)



#spatial analysis without CD23
#withoutCD23 <- pheatSet_spatial[-c(10:12, 17, 18),]
#plot a ROI-clustered heatmap
heatmap_result <- pheatmap(withoutCD23, cluster_row = FALSE, cutree_cols = 2, clustering_callback = callback, angle_col = "90")
dendrogram_result <- heatmap_result$tree_col

# Plot the dendrogram to visualize column dissimilarity
plot(dendrogram_result)



#probe analysis
pheatSet_barcode <- df_t[-c(1:20),]
#plot a ROI-clustered heatmap
heatmap_result <- pheatmap(pheatSet_barcode, cluster_row = FALSE, cutree_cols = 1, clustering_callback = callback, angle_col = "90")
dendrogram_result <- heatmap_result$tree_col

# Plot the dendrogram to visualize column dissimilarity
plot(dendrogram_result)


#for saving of arguments
# pheatmap(df_t, cutree_cols = 5, cutree_rows = 5, clustering_callback = callback, angle_col = "45", gaps_row = c(12, 20), 
#                            cutree_col = 1)

