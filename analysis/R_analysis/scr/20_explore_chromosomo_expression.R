# AIM ---------------------------------------------------------------------
# the aim of the script is explore if the gene expression from one specific clone is unbalanced chromosomo wide, compared to the others clones

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")
library("AnnotationHub")
library("ensembldb")
library("tximport")
library("GenomicFeatures")
library("txdbmaker")
# library("BSgenome")
# library("BSgenome.Hsapiens.UCSC.hg19")
# library("karyoploteR")

# read in the data --------------------------------------------------------
# read in the fully processed dds object
dds <- readRDS("../../out/object/dds_filter_DESeq.rds")

# GTF gene anntation
gtf_file <- "../../data/genes_ensembl.gtf"

# Drop non-standard chromosome names
# Create a TranscriptDb object
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# Extract gene ranges and clean it
gene_ranges <- genes(txdb)
keep_seqs <- standardChromosomes(gene_ranges)
gene_ranges_clean <- keepSeqlevels(gene_ranges, keep_seqs, pruning.mode="coarse")

# wrangling ---------------------------------------------------------------
# pull the metadata
LUT_samples <- colData(dds) %>%
  data.frame()

# Extract normalized expression data from the dds object
# Using Variance Stabilizing Transformation is generally best for plotting as it normalizes for size factors and stabilizes variance.

# Run VST on your dds object
# vsd <- vst(dds, blind=T)

# Extract the expression matrix
# expr_matrix <- assay(vsd)
expr_matrix <- counts(dds, normalized=T)

# aggregate the expression by clone
# mean_expr <- rowMeans(expr_matrix)
# define the clones
clone_ids <- LUT_samples$clone %>% unique()

list_vsd_avg <- lapply(LUT_samples$clone %>% unique(), function(x){
  # index the column
  index_col <- LUT_samples$clone %in% x
  # print(index_col)
  
  # pull the column and generate the average
  vsd_avg <- data.frame(expr_matrix[,index_col]) %>%
    rowMeans()
  
  return(vsd_avg)
}) %>%
  setNames(clone_ids)

# add the avg vsd values to the range object
df_vsd_avg_order <- pmap(list(list_vsd_avg,names(list_vsd_avg)), function(val,nm){
  # Add the mean expression data to the gene_ranges object
  # Ensure the gene IDs (row names of mean_expr) match the names of gene_ranges
  match_order <- match(names(gene_ranges_clean), names(val))
  test <- val[match_order]
  
  # add the gene names
  names(test) <- names(gene_ranges_clean)
  
  # add the column name
  test <- as.data.frame(test)
  colnames(test) <- nm
  
  return(test)
}) %>%
  bind_cols()

# add the values to the range
# gene_ranges_clean2 <- gene_ranges_clean
mcols(gene_ranges_clean) <- cbind(mcols(gene_ranges_clean),df_vsd_avg_order)

# Remove any genes that is not a complete cases
# gene_ranges_clean <- gene_ranges_clean[!is.na(gene_ranges_clean$mean_vst_expr)]
gene_index <- mcols(gene_ranges_clean) %>% complete.cases()
gene_ranges_clean2 <- gene_ranges_clean[gene_index]

# plot boxplot per gene ---------------------------------------------------
# 1. Convert GRanges to a data frame
df_expr <- as.data.frame(gene_ranges_clean2) %>%
  # Create a column for chromosome names %>%
  mutate(chr = as.character(seqnames)) %>%
  mutate(chr = factor(chr,levels = c(1:22,"X","Y","MT"))) %>%
  pivot_longer(names_to = "clone",values_to = "mean_expr",c("OSR_001","OSR_002","CTRL8","RR25","CTRL4","RR16","RR24")) %>%
  mutate(clone = fct_relevel(clone,c("OSR_001","OSR_002"))) %>%
  # dplyr::filter(seqnames == 1) %>%
  mutate(mean_expr_log1p = log1p(mean_expr))

# 2. Create the boxplot
ggplot(df_expr, aes(y=chr, x=mean_expr_log1p,fill=clone)) +
  geom_boxplot(outlier.shape = NA) + # Remove individual outliers for clarity
  # geom_jitter(alpha=0.1, color="blue") + # Add jittered points
  theme_bw() +
  labs(title = "Gene Expression Distribution by Chromosome",
       x = "Chromosome",
       y = "Mean VST Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = unname(pals::trubetskoy()))
  # scale_x_continuous(trans = "log1p")


# make line plots ---------------------------------------------------------


# Get chromosome lengths (seqlengths) and convert to a data frame

# -------------------------------------------------------------------------
# 1. Convert GRanges to a data frame
df_ranges <- as.data.frame(gene_ranges_clean2) %>%
  dplyr::rename(chr = seqnames)

# 2. Calculate the maximum 'end' position for each chromosome
calculated_lengths <- df_ranges %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    Calculated_Length = max(end),
    .groups = 'drop'
  )

# Convert the result back to a named vector
final_chr_lengths <- setNames(calculated_lengths$Calculated_Length, calculated_lengths$chr)
seqlengths(gene_ranges_clean2) <- final_chr_lengths

# -------------------------------------------------------------------------

# 2. Calculate the cumulative start position of each chromosome
chr_len <- data.frame(chr = names(seqlengths(gene_ranges_clean2)),
                      len = seqlengths(gene_ranges_clean2))


# Calculate cumulative length of all chromosomes before the current one
chr_len_cum <- chr_len %>%
  mutate(len = as.numeric(len)) %>%
  # arrange(chr) %>%
  # mutate(test = cumsum(len)) %>%
  mutate(chr_start = lag(cumsum(len), default = 0)) %>%
  dplyr::select(chr, chr_start)

# 3. Join and calculate the final genome position (pos_cum)
df_final <- df_expr %>%
  left_join(chr_len_cum, by = "chr") %>%
  # The continuous genome position is the chromosome start + gene start
  mutate(pos_cum = start + chr_start)

# Identify the midpoint for labeling each chromosome on the x-axis
axis_df <- df_final %>%
  group_by(chr) %>%
  summarize(center = mean(pos_cum))

test <- df_final
  # dplyr::filter(clone == "OSR_001",seqnames == 1)
  # dplyr::filter(clone == "OSR_001")
# 
# df_final %>%
#   dplyr::filter(clone == "OSR_002",seqnames == "Y") %>%
#   arrange(pos_cum)

ggplot(test, aes(x = pos_cum, y = mean_expr_log1p,col = clone,group=clone)) +
  # geom_line()
  geom_smooth(method = "loess", 
              span = 0.1, # Controls smoothing intensity (lower = less smooth)
              formula = y ~ x,
              alpha = 0.05
              # color = "black", 
              # fill = "lightblue",
              # linetype = "dashed"
              ) +
  
  # # 3. Add vertical lines to delineate chromosomes
  # geom_vline(xintercept = chr_len_cum$chr_start, 
  #            color = "grey50", 
  #            linetype = "dashed", 
  #            linewidth = 0.5) +
  # 
  # # 4. Set custom x-axis breaks and labels
  # scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
  facet_wrap(~ seqnames,scales = "free") +
  
  # 5. Apply theme and labels
  labs(
    title = "Smoothed Gene Expression Across the Whole Genome",
    x = "Chromosome",
    y = "Mean Normalized Expression"
  ) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1)
  ) +
  scale_color_manual(values = unname(pals::trubetskoy()))
ggsave("../../out/plot/20_reads_across_genome.pdf",width = 20,height = 20)

# simulate ----------------------------------------------------------------
# simulate for one clone a big portion of genes missing (set the expression to 0)
test2 <- df_final %>%
  mutate(mean_expr_log1p = case_when(clone == "OSR_002" & seqnames == 16 & pos_cum > 2420000000 & pos_cum < 2440000000 ~ 
                                       0,
                                     T ~ mean_expr_log1p))
# dplyr::filter(clone == "OSR_001",seqnames == 1)
# dplyr::filter(clone == "OSR_001")
# 
# df_final %>%
#   dplyr::filter(clone == "OSR_002",seqnames == "Y") %>%
#   arrange(pos_cum)

ggplot(test2, aes(x = pos_cum, y = mean_expr_log1p,col = clone,group=clone)) +
  # geom_line()
  geom_smooth(method = "loess", 
              span = 0.1, # Controls smoothing intensity (lower = less smooth)
              formula = y ~ x,
              alpha = 0.05
              # color = "black", 
              # fill = "lightblue",
              # linetype = "dashed"
  ) +
  
  # # 3. Add vertical lines to delineate chromosomes
  # geom_vline(xintercept = chr_len_cum$chr_start, 
  #            color = "grey50", 
  #            linetype = "dashed", 
  #            linewidth = 0.5) +
  # 
  # # 4. Set custom x-axis breaks and labels
  # scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
  facet_wrap(~ seqnames,scales = "free") +
  
  # 5. Apply theme and labels
  labs(
    title = "Smoothed Gene Expression Across the Whole Genome",
    x = "Chromosome",
    y = "Mean Normalized Expression"
  ) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1)
  ) +
  scale_color_manual(values = unname(pals::trubetskoy()))
ggsave("../../out/plot/20_reads_across_genome_simulation.pdf",width = 20,height = 20)
