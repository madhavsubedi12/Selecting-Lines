# Selecting-Lines

#Making PCA plot of 100 lines to be selected using kinship matrix and details of the accessions

### Required packages
install.packages("tidyverse")
install.packages("factoextra")
install.packages("ggplot2")

library(tidyverse)
library(factoextra)
library(ggplot2)


### Load your kinship matrix
kinship_data <- read.delim("KinshipLARGElist.relatedness2",header = TRUE, sep = "\t")
dim(kinship_data) #Since it was created from a vcf file, it is not in the matrixc format we need. 


### getting the individual lines being compared
individuals <- unique(c(kinship_data$INDV1, kinship_data$INDV2)) 
n <- length(individuals)
#creating a square matrix to fill
kinship_matrix <- matrix(0, nrow=n, ncol=n, dimnames=list(individuals, individuals))

### Filling kinship matrix as necessary
for (i in 1:nrow(kinship_data)) {
  indv1 <- kinship_data$INDV1[i]
  indv2 <- kinship_data$INDV2[i]
  phi <- kinship_data$RELATEDNESS_PHI[i]
  
*** Fill symmetrically
  kinship_matrix[indv1, indv2] <- phi
  kinship_matrix[indv2, indv1] <- phi
}


###Check the structure of the matrix
str(kinship_matrix)

### Ensure the matrix is numeric
kinship_matrix <- as.matrix(kinship_matrix)

# Perform PCA on the kinship matrix
pca <- prcomp(kinship_matrix, center = TRUE, scale. = TRUE)


### Get PCA results
pca_df <- as.data.frame(pca$x)

# Add additional information about the accessions if available
### CSV file with accession information

accession_info <- read.delim("Accessionscategorized.txt")

### Merging the two data

### For this sort kinship_matrix by row names alphabetically
pca_sorted <- pca_df[order(rownames(pca_df)), ]

### Ensure Accessions is character in accession_info
accession_info$Accessions <- as.character(accession_info$Accessions)

### Sort accession_info by Accessions column alphabetically
accession_info_sorted <- accession_info[order(accession_info$Accessions), ]

rownames(pca_sorted)==accession_info_sorted$Accessions #checking if the order is correct

### Finally merging the data
merged_data <- cbind(accession_info_sorted, pca_sorted)
nrow(merged_data)

# Plotting
### Define your color values
color_values <- c("blue", "black", "pink")  # Adjust based on your unique categories in 'Identified.by'

### Identify the top 10 most frequent categories in 'Origin'
top_origin <- names(sort(table(merged_data$Origin), decreasing = TRUE))[1:10]

### Create a new category 'Others' for less frequent categories
merged_data <- merged_data %>%
  mutate(Origin_grouped = ifelse(Origin %in% top_origin, as.character(Origin), "Others"))
### PCA plot
ggplot(merged_data, aes(x = PC1, y = PC2, color = Identified.by, shape = Origin_grouped)) +  
  geom_point(size = 3) +  
  labs(title = "PCA of 1981 Hop Accessions",  
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 2), "%)"),  
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 2), "%)")) +  
  theme_minimal() +  
  theme(  
    legend.position = "right",  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 14),  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 11),  
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey80"),  
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "grey90"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1) #creates border line) +  
  coord_cartesian() +  # Use coord_cartesian to ensure all data points are plotted  
  scale_color_manual(values = color_values) +  # Manual color scale  
  scale_shape_manual(values = 1:length(unique(merged_data$Origin_grouped))) +  # All unique shapes  
  guides( color = guide_legend(title.position = "top", title = "Identified.by"),  
    shape = guide_legend(title.position = "top", title = "Origin")  
  )  
 
