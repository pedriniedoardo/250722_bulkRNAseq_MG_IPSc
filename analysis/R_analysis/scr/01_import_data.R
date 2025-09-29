# AIM ---------------------------------------------------------------------
# the aim of the script is to read in the raw table of counts and build the dds object

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")
library("AnnotationHub")
library("ensembldb")
library("tximport")
library("GenomicFeatures")
library("txdbmaker")

# read in the data --------------------------------------------------------
# read in the expressino data
test <- readRDS("../../out/object/raw_count_salmon_GTF.rds")

# remember that the test object derive frmo the txi rounding
# txi2 <- readRDS("../../out/object/txi_GTF.rds")
# data2 <- txi2$counts %>% 
#   round()

# read in the metadata
# build the annotation besed on the sample metadata
LUT_samples <- read_csv("../../data/LUT_sample.csv")

# extract only the count information from the samples of interest
mat_exp <- test %>%
  as.data.frame() %>%
  dplyr::select(LUT_samples$sample_id) %>% 
  as.matrix()

# match the order of the sample in the matrix with the sample in the sample sheet
# build the metadata for the deseq object
coldata <- data.frame(sample = colnames(mat_exp)) %>% 
  left_join(LUT_samples,by = c("sample"="sample_id")) %>% 
  mutate(rowname = sample) %>% 
  column_to_rownames("rowname")

# save the table
write.csv(coldata,file = "../../data/LUT_samples_final.csv",row.names = T)

# define the model --------------------------------------------------------
clone <- coldata$clone
# gender <- coldata$gender
# batch <- coldata$batch
treat <- factor(coldata$treat,levels = c("untreated","CSF","cytokine"))

# build the design
# design <- model.matrix(~ clone + batch + gender + treat)
design <- model.matrix(~ clone + treat)
colnames(design)[1] <- c("intercept")

saveRDS(design,file = "../../out/object/design.rds")

# build the object --------------------------------------------------------
# is keeping only the objext in the lut_sample

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = coldata,
                              design = design)
# dds2 <- DESeqDataSetFromTximport(txi2,
#                                  colData = coldata,
#                                  design = design)

saveRDS(dds,file = "../../out/object/dds.rds")

# confirm dds and dds2 are teh same
# identical(dds,dds2)
