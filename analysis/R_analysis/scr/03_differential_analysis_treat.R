# AIM ---------------------------------------------------------------------
# run DE analysis

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(limma)
library(ComplexHeatmap)
library(ashr)
library(AnnotationDbi)
library(AnnotationHub)

# read in the data --------------------------------------------------------
# read in the biotype annotation
LUT_gene <- readRDS("../../out/object/annotation_GTF.rds") %>%
  dplyr::select(GENEID,SYMBOL) %>%
  distinct()

ddsHTSeq_filter <- readRDS("../../out/object/dds_filter.rds")
design <- readRDS("../../out/object/design.rds")
vds_filter <- readRDS("../../out/object/vds_filter.rds")

LUT_sample <- colData(ddsHTSeq_filter) %>% 
  data.frame()

# differential expression analyisis ---------------------------------------
ddsHTSeq_filter <- DESeq(ddsHTSeq_filter)
# if needed is possible to check the distributions of the counts before and after the normalizatoin
# boxplot(log(counts(ddsHTSeq_structure,normalized = T)))
# boxplot(log(counts(ddsHTSeq_structure,normalized = F)))

# save the filtered object
saveRDS(ddsHTSeq_filter,"../../out/object/dds_filter_DESeq.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(CSFvsCTRL = treatCSF,
                          CytokinevsCTRL = treatcytokine,
                          levels = design)

res_cytokine <- results(ddsHTSeq_filter, contrast=contrast[,"CytokinevsCTRL"],alpha = 0.05)
res_CSF <- results(ddsHTSeq_filter, contrast=contrast[,"CSFvsCTRL"],alpha = 0.05)
summary(res_cytokine)
summary(res_CSF)

# add the gene symbols
list_df <- 
  list("cytokine_vs_ctrl" = res_cytokine,
       "CSF_vs_ctrl" = res_CSF) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "GENEID") %>%
      arrange(pvalue) %>%
      # add the symbol
      left_join(LUT_gene,by = c("GENEID"))
    df
  })

# save the tables
pmap(list(list_df,names(list_df)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# shrink ------------------------------------------------------------------  
res_cytokine_shr <- lfcShrink(ddsHTSeq_filter, res = res_cytokine, type = "ashr")
res_CSF_shr <- lfcShrink(ddsHTSeq_filter, res = res_CSF, type = "ashr")
summary(res_cytokine_shr)
summary(res_CSF_shr)

# add the gene symbols
list_df_shr <- 
  list("cytokine_vs_ctrl_shr" = res_cytokine_shr,
       "CSF_vs_ctrl_shr" = res_CSF_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "GENEID") %>%
      arrange(pvalue) %>%
      # add the symbol
      left_join(LUT_gene,by = c("GENEID"))
    df
  })


# save the tables
pmap(list(list_df_shr,names(list_df_shr)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
list_df$cytokine_vs_ctrl %>%
  data.frame()%>%
  dplyr::filter(baseMean>1)%>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  theme_bw()
ggsave("../../out/plot/histogram_pvalue_cytokine.pdf",width = 4,height = 3)

list_df$CSF_vs_ctrl %>%
  data.frame()%>%
  dplyr::filter(baseMean>1)%>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  theme_bw()
ggsave("../../out/plot/histogram_pvalue_CSF.pdf",width = 4,height = 3)

# PLOTTING RESULTS --------------------------------------------------------
# add the info of the genename
test_plot <- list_df$cytokine_vs_ctrl%>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(SYMBOL)),yes = 1,no = 0)) %>%
  dplyr::filter(!is.na(col))

test_plot %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot[test_plot$col==1,],
    aes(label = SYMBOL),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/plot/vulcano_plot_text_cytokine.pdf",width = 12,height = 12)
#
test_plot_shr <- list_df_shr$cytokine_vs_ctrl_shr%>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(SYMBOL)),yes = 1,no = 0)) %>%
  dplyr::filter(!is.na(col))

test_plot_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot_shr[test_plot_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot_shr[test_plot_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot_shr[test_plot_shr$col==1,],
    aes(label = SYMBOL),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/plot/vulcano_plot_text_cytokine_shr.pdf",width = 12,height = 12)

# plotMA(res, ylim = c(-5, 5))
list_df$cytokine_vs_ctrl %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) + theme_bw() + 
  theme(legend.position = "none")
ggsave("../../out/plot/MA_plot_cytokine.pdf",width = 4,height = 3)

# using the shrinked values
list_df_shr$cytokine_vs_ctrl_shr %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color))+geom_point(alpha=0.2) + 
  scale_x_log10() + 
  scale_color_manual(values = c("gray","red")) + 
  theme_bw() +
  theme(legend.position = "none")
ggsave("../../out/plot/MA_plot_shr_cytokine.pdf",width = 4,height = 3)

# the DEGs plot stringent
df_DEG_1 <- list_df$cytokine_vs_ctrl %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1)

mat_filter <- assay(vds_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# mat <- mat_filter[rownames(vds_filter) %in% df_DEG_1$GENEID, ]
mat <- mat_filter[df_DEG_1$GENEID, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)
# change the rownames with gene symbols
rownames(mat2) <- df_DEG_1$SYMBOL
#

sample_ordered <- data.frame(sample = colnames(mat2)) %>%
  left_join(ddsHTSeq_filter@colData %>%
              data.frame(),by="sample")
  

# update the column name of the matrix
colnames(mat2) <- sample_ordered$sample

column_ha <- HeatmapAnnotation(treat = sample_ordered$treat,
                               gender = sample_ordered$gender,
                               col = list(treat = c("untreated" = "gray", "CSF"="cyan","cytokine" = "black"),
                                          gender = c("M" = "blue", "F" = "pink"))) 

ht2 <- Heatmap(mat2, 
               name = "exp", 
               column_title = "treat",
               row_names_gp = gpar(fontsize = 5),
               top_annotation = column_ha, show_row_names = T
               # cluster_rows = F, 
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7))
) 
pdf("../../out/plot/heatmap_DEG_cytokine.pdf",width = 5,height = 6) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# the DEGs plot shr
df_DEG_2 <- list_df_shr$cytokine_vs_ctrl_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1)

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

# mat_shr <- mat_filter[rownames(vds_filter) %in% df_DEG_2$GENEID, ]
mat_shr <- mat_filter[df_DEG_2$GENEID, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr)
# change the rownames with gene symbols
rownames(mat2_shr) <- df_DEG_2$SYMBOL

#

sample_ordered_shr <- data.frame(sample = colnames(mat2_shr)) %>%
  left_join(ddsHTSeq_filter@colData %>%
              data.frame(),by="sample")

# update the column name of the matrix
colnames(mat2_shr) <- sample_ordered_shr$sample

column_ha_shr <- HeatmapAnnotation(treat = sample_ordered_shr$treat,
                                   gender = sample_ordered_shr$gender,
                                   col = list(treat = c("untreated" = "gray", "CSF"="cyan","cytokine" = "black"),
                                              gender = c("M" = "blue", "F" = "pink"))) 


ht2_shr <- Heatmap(mat2_shr, 
                   name = "exp", 
                   column_title = "treat shr",
                   row_names_gp = gpar(fontsize = 5),
                   top_annotation = column_ha_shr, show_row_names = T
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
) 
pdf("../../out/plot/heatmap_DEG_treat_shr.pdf",width = 5,height = 6) 
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/plot/ddsHTSeq_filter_dispersion.pdf",width = 5,height = 5) 
plotDispEsts(ddsHTSeq_filter)
dev.off()
