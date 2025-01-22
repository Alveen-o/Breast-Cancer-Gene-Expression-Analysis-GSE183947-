# **Breast Cancer Gene Expression Analysis (GSE183947)**

## **Overview**
This project involves the analysis of RNA-seq data from the **GSE183947** dataset, which focuses on breast cancer gene expression. The study includes differential expression analysis, data visualization, and comparisons of expression patterns for key cancer-associated genes across tissues and metastasis conditions.

## **Dataset**
- **Source**: [GEO Database - GSE183947](https://www.ncbi.nlm.nih.gov/geo/)
- **Data Description**:  
  - FPKM-normalized expression values for genes across breast cancer samples and normal tissues.
  - Metadata includes tissue type and metastasis status for each sample.

---

## **Project Workflow**
### 1. **Data Preprocessing**
- **Input Files**:
  - Gene expression data: `GSE183947_fpkm.csv`
  - Metadata retrieved using `GEOquery` in R.
- **Steps**:
  - Data reshaped from wide to long format for easier analysis using `tidyverse`.
  - Metadata cleaned and merged with gene expression data.

### 2. **Analysis**
- Filtered for genes of interest, including:
  - **BRCA1**, **BRCA2**, **TP53**, **ALK**, and **MYCN**.
- Calculated statistical summaries (mean, median) for FPKM values by gene and tissue.

### 3. **Visualizations**
Generated multiple plots using `ggplot2` to explore and visualize the data:

#### **Bar Plot**
- **Description**: Displays the expression of BRCA1 across samples grouped by tissue type.
- **File**: [Bar Chart](barchart.pdf)

#### **Density Plot**
- **Description**: Shows the distribution of BRCA1 expression levels across tissues.
- **File**: [Density Plot](density.pdf)

#### **Box Plot**
- **Description**: Compares BRCA1 expression levels across metastasis conditions.
- **File**: [Box Plot](boxplot.pdf)

#### **Scatter Plot**
- **Description**: Correlates the expression levels of BRCA1 and BRCA2 across samples, with points colored by tissue type.
- **File**: [Scatter Plot](scatterplot.pdf)

#### **Heatmap**
- **Description**: Visualizes the expression levels of selected genes (**BRCA1**, **BRCA2**, **TP53**, **ALK**, **MYCN**) across all samples.
- **File**: [Heatmap](heatmap.pdf)

---

## **Files**
### **Scripts**
- [R Script for Analysis](gene_expression.R): Contains all the R code used for data preprocessing, analysis, and visualization.

### **Data**
- `GSE183947_fpkm.csv`: Gene expression dataset.

### **Plots**
- [Bar Chart](barchart.pdf)
- [Density Plot](density.pdf)
- [Box Plot](boxplot.pdf)
- [Scatter Plot](scatterplot.pdf)
- [Heatmap](heatmap.pdf)

---

## **How to Reproduce**
1. **Dependencies**:
   - R version ≥ 4.0.0
   - Required R libraries:
     ```R
     install.packages(c("dplyr", "tidyverse", "GEOquery", "ggplot2", "conflicted"))
     ```

2. **Run the script**:
   ```R
   source("gene_expression.R")
   ```

---

## **Insights**
1. **Gene Expression Patterns**:
   - **BRCA1** and **BRCA2** showed distinct expression profiles in tumor and normal tissues.
   - **TP53** displayed consistent expression across most samples.
2. **Metastasis and Gene Expression**:
   - Significant differences in **BRCA1** expression between metastatic and non-metastatic samples.

---

## **License**
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## **Acknowledgments**
- GEO Database for providing the data.
- R community for the tools used in this analysis.

