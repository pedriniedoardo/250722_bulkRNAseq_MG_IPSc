# AIM ---------------------------------------------------------------------
# comfirm the gender form the transcriptomic profile

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)
library(scales)
library(pals)

# plot genes with heatmap -------------------------------------------------
vds <- readRDS(file = "../../out/object/vds_unfilter.rds")

# import the annotation file for the genes
annotations <- readRDS("../../out/object/annotation_GTF.rds") %>%
  dplyr::select(GENEID,SYMBOL) %>%
  distinct()

df_GOI <- annotations %>%
  dplyr::filter(SYMBOL %in% c("XIST","DDX3Y","RPS4Y1","USP9Y"))
  
# gene_id <- rownames(assay(vds)) %in% df_GOI$GENEID
# mat <- assay(vds)[gene_id, ]

mat <- assay(vds)[df_GOI$GENEID, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat,useNames = TRUE)

# add the gene symbol to the matrix
rownames(mat2) <- df_GOI$SYMBOL

# build the annotation object  
meta <- colData(vds) %>% 
  data.frame()

LUT_sample_matrix <- data.frame(sample_matrix = colnames(mat2)) %>%
  left_join(meta,by = c("sample_matrix"="sample")) %>%
  mutate(test = paste0(clone,"_",treat))

# sample_ordered <- LUT_sample_matrix$clone

# define a LUT for the clones
color_id <- alphabet(length(unique(LUT_sample_matrix$clone)))
# check the colors
show_col(color_id)

# build the named vector
names(color_id) <- unique(LUT_sample_matrix$clone)

column_ha <- HeatmapAnnotation(clone = LUT_sample_matrix$clone,
                               gender = LUT_sample_matrix$gender,
                               treat = LUT_sample_matrix$treat,
                               col = list(clone = color_id,
                                          gender = c("M" = "blue", "F" = "pink"))) 

# change the name of the column in the matrix
colnames(mat2) <- LUT_sample_matrix$test

# row_ha <- rowAnnotation(class = rep(c("common B_D","offtarget B","offset D"),c(22,3,24)),
#                         col = list(class = c("common B_D" = "violet", "offtarget B" = "yellow","offset D"="brown")))

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "gender genes") 

pdf("../../out/plot/heatmap_GOI_gender.pdf",width = 8,height = 5) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()