# AIM ---------------------------------------------------------------------
# run EnrichR with a panel of annotations

# libraries ---------------------------------------------------------------
library(enrichR)
library(tidyverse)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  dplyr::filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  dplyr::filter(str_detect(libraryName,pattern = "KEGG"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2022","Human_Gene_Atlas","Azimuth_Cell_Types_2021")

# GENE SELECTION ----------------------------------------------------------
# save the ranked object, also change the genes into genenames
full_df <- read_tsv("../../out/table/res_BMP9_vs_Mock_shr.txt")

list_res_tot <- list(BMP9_vs_mock = full_df %>%
                       dplyr::filter(padj<0.05,abs(log2FoldChange)>1) %>%
                       pull(symbol))

# query -------------------------------------------------------------------
list <- lapply(list_res_tot,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  #
  out_enrich
})

df_enrichr_annotation_enriched_tot <- list$BMP9_vs_mock %>%
  bind_rows(.id = "annotation")

df_enrichr_annotation_enriched_tot %>%
  write_tsv("../../out/table/enrichR_out_BMP9_vs_mock_shr.tsv")

library(scales)
df_enrichr_annotation_enriched_tot %>%
  group_by(annotation) %>%
  arrange(P.value) %>%
  dplyr::slice(1:20) %>%
  mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
  mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
  ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,ncol = 1,scales = "free")+theme_bw() +
  scale_color_gradientn(colors = c("red","blue"),
                        values = rescale(c(0,1)),
                        limits = c(0,0.2))+theme(strip.background = element_blank())
# scale_color_gradient(low = "red",high = "blue")
ggsave("../../out/plot/enrichR_out_BMP9_vs_mock_shr_plot.pdf",width = 7,height = 15)
