# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(conflicted)
library(ggplot2)
# read in dataset
data <- read_csv("GSE183947_fpkm.csv")

#get the metadata file
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE) # downloads and retrieve data from the GEO database
gse

# Extract the sample-level metadata (phenotypic data) from the first ExpressionSet in the GEO dataset
# - `gse[[1]]`: Access the first ExpressionSet object
# - `phenoData()`: Retrieves the phenotypic metadata slot (sample annotations)
# - `pData()`: Converts the phenotypic metadata into a data frame format for easier analysis
metadata <- pData(phenoData(gse[[1]]))
head(metadata)


metadata_mod <- metadata%>%
select(c(1,10,11,17))%>% # selects the desired column
rename(tissue = characteristics_ch1, metastasis = characteristics_ch1.1)%>% #renames the columns
mutate(tissue = gsub("tissue: ", "", tissue), metastasis = gsub("metastasis: ", "", metastasis)) #removes specific prefixes from the text values in each column
head(metadata_mod)


#reshape data
head(data)
# Convert gene expression data from wide format to long format for easier analysis
# - Rename the first column ("...1") to "gene" for clarity
data.long <- data%>%
  rename(gene = ...1 )%>%
  # Gather all other columns into two new columns: "samples" and "FPKM"
  # - "samples": Original column names (e.g., sample IDs)
  # - "FPKM": Expression values for each gene-sample pair
  gather(key = "samples", value = "FPKM", -gene)

# Merge additional sample metadata into the long-format gene expression data
# - Perform a left join to add columns from "metadata_mod" to "data.long"
# - Match rows where the "samples" column in "data.long" equals the "description" column in "metadata_mod"
data.long <- data.long%>%
  left_join(.,metadata_mod, by = c("samples" = "description"))


#explore data
# Filter the data to include only rows for BRCA1 and BRCA2 genes
data.long%>%
  filter(gene == "BRCA1" | gene == "BRCA2")%>%
# Group the data by gene and tissue to analyze expression levels within each group
  group_by(gene, tissue)%>%
# Calculate the mean and median FPKM values for each group
  summarise(mean_FPKM = mean(FPKM), median_FPKM = median(FPKM))%>%
  # Arrange the results in ascending order of mean FPKM
  arrange(mean_FPKM)
  
#visualize data
  #bar chart
  # Filter the dataset to include only rows where the gene is "BRCA1"
a <- data.long%>%# 
  filter(gene == "BRCA1")%>%
  # Pass the filtered data to ggplot for visualization
  ggplot(.,aes(x = samples, y = FPKM, fill = tissue)) + 
  # Create a bar plot (geom_col() is used to create bar charts)
    geom_col()
ggsave(a, filename = "barchart.pdf", width = 10, height = 8)
  #Density
b <- data.long%>%
  filter(gene == "BRCA1")%>%# Filter for rows where the gene is BRCA1
  ggplot(.,aes(x = FPKM, fill = tissue)) + # Map FPKM to x-axis and tissue to fill color
  geom_density(alpha = 0.5)# Create overlapping density curves with transparency (alpha)
ggsave(b, filename = "density.pdf", width = 10, height = 8)

  #boxplot
c <- data.long%>%
  filter(gene == "BRCA1")%>%# Filter for rows where the gene is BRCA1
  ggplot(.,aes(x = metastasis, y =FPKM)) +# Map metastasis status to x-axis and FPKM to y-axis
  geom_boxplot()# Create a box plot
ggsave(c, filename = "boxplot.pdf", width = 10, height = 8)

  #scatterplot
d <- data.long%>%
  filter(gene == "BRCA1" | gene == "BRCA2")%>%# Filter for BRCA1 and BRCA2 genes
  spread(key = gene, value = FPKM)%>%# Reshape data to have separate columns for BRCA1 and BRCA2 expression
  ggplot(.,aes(x = BRCA1, y = BRCA2, colour = tissue)) +# Map BRCA1 to x-axis, BRCA2 to y-axis, and tissue to color
  geom_point() +# Add points to represent samples
  geom_smooth(method = "lm", se = FALSE)# Add a linear regression line without confidence interval (se = FALSE)
ggsave(d, filename = "scatterplot.pdf", width = 10, height = 8)


  #heatmap
genes.of.interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")# List of genes to analyze
e <- data.long%>%
  filter(gene %in% genes.of.interest)%>%# Filter for rows where the gene is in the list of interest
  ggplot(.,aes(x = samples, y = gene, fill = FPKM)) +# Map samples to x-axis, genes to y-axis, and FPKM to tile color
  geom_tile() +# Create heatmap tiles
  scale_fill_gradient(low = "white", high = "blue")# Set color gradient from white (low expression) to blue (high expression)
ggsave(e, filename = "heatmap.pdf", width = 10, height = 8)  
