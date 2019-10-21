# Shinydashboard-DESeq2
A Shinydashboard for visualizing RNAseq data analyzed with DESeq2
For demonstration purposes, the script downloads the countmatrix and phenodata from Mov10 data made available at the Harvard Chan Bioinformatic Core (https://github.com/hbctraining/DGE_workshop). This app uses some of the analytical tools employed by the workshop. Thus, for further details on explanation of the DESeq2 analysis, please refer to the HBC workshop link. The material and presentation of workshop was very insightful. The shinydashboard codes were synthesized from the many talented answers to questions on stackoverflow and the shinydashboard help website.

Any questions or feedback are welcome. Feel free to modify and repurpose the app. There are many other well developed analytical and visualization tools available.



libraries required for the app. So ensure that they are installed. Some are directly from CRAN and others are to be installed via BiocManager from Bioconductor.
## From CRAN
example:
install.packages("shiny")

## From Bioconductor
example:

if (!requireNamespace("BiocManager", quietly = TRUE))
    
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## Required Libraries
library(shiny)
library(shinydashboard)
library(shinyjs) # for dynamically displaying tab titles
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library("ggrepel") #Avoid overlapping labels
library(scales)
library(DT)
library(tidyr)

OTHER ATTACHED PACKAGES
jsonlite_1.6    Matrix_1.2-17      viridisLite_0.3.0    DT_0.9  
 shiny_1.4.0    shinyjs_1.0     shinydashboard_0.7.1    
ggfortify_0.4.7     data.table_1.12.4       tidyr_1.0.0  
ggrepel_0.8.1       dplyr_0.8.3     magrittr_1.5  
DESeq2_1.24.0         pheatmap_1.0.12     
SummarizedExperiment_1.14.1     DelayedArray_0.10.0     BiocParallel_1.18.1 
matrixStats_0.55.0      Biobase_2.44.0 
GenomicRanges_1.36.1    GenomeInfoDb_1.20.0     IRanges_2.18.3      BiocGenerics_0.30.0 
S4Vectors_0.22.1    RColorBrewer_1.1-2          ggplot2_3.2.1    scales_1.0.0          
