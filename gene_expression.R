# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(conflicted)

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
  
  
  
  
  