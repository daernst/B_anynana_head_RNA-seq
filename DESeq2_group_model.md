## Analyze differential gene expression between male and female adults:

```
### Load required packages ###
library(DESeq2))
library(apeglm)
library(dplyr)


### Read in table containing sample and read counts file information (sampleName, fileName, family, sex, stage, and grouping variable) ###
sampleTable <- read.csv("D:/Documents/B_anynana_Larva_Adult_Head_Study/sampleTable.csv")


### Build DESeqDataSet using htseq-count output files ###
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "D:/Documents/B_anynana_Larva_Adult_Head_Study/Revised_study_for_BMC_Genomics/htseq-count_output_CORRECT",
                                       design = ~ group)


### Set female adult as the reference level ###
ddsHTSeq$group <- relevel(ddsHTSeq$group, 
                          ref = "Female_Adult")


### Filter out all genes with <10 reads total across all samples ###
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


### Run differential expression analysis ###
dds_adult <- DESeq(ddsHTSeq)


### Shrink LFC to obtain more accurate effect sizes ###
res.lfcShrink_adult <- lfcShrink(dds_adult, 
                                 coef="group_Male_Adult_vs_Female_Adult", 
                                 type="apeglm")


### Get summary of results ###
summary(res.lfcShrink_adult, 
        alpha = 0.05)


### Get total number of differentially expressed genes ###
sum(res.lfcShrink_adult$padj < 0.05,
    na.rm=TRUE)


### Add functional annotations to results ###
annot <- read.csv("Bicyclus_anynana_v1_2_genes_only_b2g_functional_annotation_B2G-Desc_and_Best_Hit_3.csv")
res.lfcShrink_adult$SeqName <- row.names(res.lfcShrink_adult)
res.lfcShrink_adult <- as.data.frame(res.lfcShrink_adult) %>% 
  right_join(annot, ., by="SeqName")


### Order all results by FDR ###
res.lfcShrink_adult.ordered <- res.lfcShrink_adult[order(res.lfcShrink_adult$padj),]


### Extract all differentially expressed genes ###
res.lfcShrink_adult_05 <- res.lfcShrink_adult.ordered %>% 
  filter(padj < 0.05)


### Order all differentially expressed genes by effect size (the absolute value of log2FoldChange) ###
res.lfcShrink_adult_05_lfc <- res.lfcShrink_adult_05[order(-abs(res.lfcShrink_adult_05$log2FoldChange)),]


### Extract genes upregulated in male adults ###
res.lfcShrink_adult_05_male_up <- res.lfcShrink_adult_05 %>% 
  filter(log2FoldChange > 0)


### Extract genes upregulated in female adults ###
res.lfcShrink_adult_female_up <- res.lfcShrink_adult_05 %>% 
  filter(log2FoldChange < 0)
```

## Analyze differential gene expression between male and female larvae:

```
### Set female larva as the reference level ###
ddsHTSeq$group <- relevel(ddsHTSeq$group, 
                          ref = "Female_Larva")


### Run differential expression analysis ###
dds_larva <- DESeq(ddsHTSeq)


### Shrink LFC to obtain more accurate effect sizes ###
res.lfcShrink_larva <- lfcShrink(dds_larva, 
                                 coef="group_Male_Larva_vs_Female_Larva", 
                                 type="apeglm")


### Get summary of results ###
summary(res.lfcShrink_larva, 
        alpha = 0.05)


### Get total number of differentially expressed genes ###
sum(res.lfcShrink_larva$padj < 0.05,
    na.rm=TRUE)


### Add functional annotations to results ###
res.lfcShrink_larva$SeqName <- row.names(res.lfcShrink_larva)
res.lfcShrink_larva <- as.data.frame(res.lfcShrink_larva) %>% 
  right_join(annot, ., by="SeqName")


### Order all results by FDR ###
res.lfcShrink_larva.ordered <- res.lfcShrink_larva[order(res.lfcShrink_larva$padj),]


### Extract all differentially expressed genes ###
res.lfcShrink_larva_05 <- res.lfcShrink_larva.ordered %>% 
  filter(padj < 0.05)


### Order all differentially expressed genes by effect size (the absolute value of log2FoldChange) ###
res.lfcShrink_larva_05_lfc <- res.lfcShrink_larva_05[order(-abs(res.lfcShrink_larva_05$log2FoldChange)),]


### Extract genes upregulated in male larvae ###
res.lfcShrink_larva_05_male_up <- res.lfcShrink_larva_05 %>% 
  filter(log2FoldChange > 0)


### Extract genes upregulated in female larvae ###
res.lfcShrink_larva_05_female_up <- res.lfcShrink_larva_05 %>% 
  filter(log2FoldChange < 0)
```
