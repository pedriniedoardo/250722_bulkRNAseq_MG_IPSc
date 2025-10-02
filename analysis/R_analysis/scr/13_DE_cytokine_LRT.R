# AIM ---------------------------------------------------------------------
# run DE analysis using the LRT for the cytokine sample

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(limma)
library(ComplexHeatmap)
library(ashr)
library(AnnotationDbi)
library(AnnotationHub)
library(DEGreport)

# read in the data --------------------------------------------------------
# read in the biotype annotation
LUT_gene <- readRDS("../../out/object/annotation_GTF.rds") %>%
  dplyr::select(GENEID,SYMBOL) %>%
  distinct()

dds_full <- readRDS("../../out/object/dds_cytokine_LRT_filter.rds")
design_full <- readRDS("../../out/object/design_cytokine_LRT.rds")
vds_filter <- readRDS("../../out/object/vds_cytokine_LRT_filter.rds")

LUT_sample <- colData(dds_full) %>% 
  data.frame()

# LRT ---------------------------------------------------------------------

# build the reduce model
coldata <- colData(dds_full)
clone <- coldata$clone

# build the design
# design <- model.matrix(~ clone + batch + gender + treat)
design_reduced <- model.matrix(~ clone)
colnames(design_reduced)[1] <- c("intercept")

saveRDS(design_reduced,file = "../../out/object/design_reduced_cytokine_LRT.rds")

# Likelihood ratio test
dds_lrt <- DESeq(dds_full, test="LRT", reduced = design_reduced)

# save the LRT object
saveRDS(dds_lrt,"../../out/object/dds_cytokine_LRT_filter.rds")

# Extract results for LRT
res_LRT <- results(dds_lrt)

# resultsNames(dds)
resultsNames(dds_lrt)

# When filtering significant genes from the LRT we threshold only the padj column.
padj.cutoff <- 0.05

# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="GENEID") %>% 
  left_join(LUT_gene,by = "GENEID") %>%
  as_tibble() %>%
  arrange(padj)

# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  dplyr::filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)

# Identifying clusters of genes with shared expression profiles
# We now have this list of significant genes that we know are changing in some way across the three different sample groups.
# A good next step is to identify groups of genes that share a pattern of expression change across the sample groups (levels). To do this we will be using a clustering tool called degPatterns from the ‘DEGreport’ package. The degPatterns tool uses a hierarchical clustering approach based on pair-wise correlations between genes, then cuts the hierarchical tree to generate groups of genes with similar expression profiles. The tool cuts the tree in a way to optimize the diversity of the clusters, such that the variability inter-cluster > the variability intra-cluster.

# Before we begin clustering, we will first subset our vst transformed normalized counts to retain only the differentially expressed genes (padj < 0.05). In our case, it may take some time to run the clustering genes.
vsd_mat <- assay(vds_filter)
head(vsd_mat)

# Obtain vsd values for those significant genes
cluster_vsd <- vsd_mat[sigLRT_genes$GENEID, ]
head(cluster_vsd)

# The vsd transformed counts for the significant genes are input to degPatterns along with a few additional arguments:
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_vsd,
                        metadata = coldata %>% as.data.frame() %>% mutate(exposure = factor(exposure,levels = c("0","6","24"),labels = c("00","06","24"))),
                        time = "exposure",
                        # integer minimum number of genes in a group that will be return
                        minc = 2,
                        col = NULL)

# check the structure of the object
str(clusters,max.level = 1)

# Once the clustering is finished running, you will get your command prompt back in the console and you should see a figure appear in your plot window. The genes have been clustered into four different groups. For each group of genes, we have a boxplot illustrating expression change across the different sample groups. A line graph is overlayed to illustrate the trend in expression change.

# What type of data structure is the `clusters` output?
class(clusters)

# We can see what objects are stored in the list by using names(clusters). There is a dataframe stored inside. This is the main result so let’s take a look at it. The first column contains the genes, and the second column contains the cluster number to which they belong.

str(clusters$df)
dim(clusters$df)

# add the informations to the table of statistics
dim(sigLRT_genes)

sigLRT_genes_clusters <- sigLRT_genes %>%
  left_join(clusters$df, by = c("GENEID" = "genes")) %>%
  arrange(cluster, padj)

# save the table for the genes
sigLRT_genes_clusters %>%
  write_tsv("../../out/table/res_LRT_cytokine.tsv")

# split the table per group
sigLRT_genes_clusters %>%
  split(f = .$cluster)
