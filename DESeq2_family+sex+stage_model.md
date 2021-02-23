## Analyze differential gene expression between adult and larval *B. anynana* while controlling for differences in expression due to lineage and sex:

```
### Load required packages ###
library(DESeq2)
library(apeglm)
library(readr)
library(dplyr)


### Read in table containing sample and read counts file information (sampleName, fileName, family, sex, stage, and grouping variable) ###
sampleTable <- read_csv("D:/Documents/B_anynana_Larva_Adult_Head_Study/sampleTable.csv")


### Build DESeqDataSet using htseq-count output files ###
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "D:/Documents/B_anynana_Larva_Adult_Head_Study/htseq-count_output",
                                       design = ~ family + sex + stage)


### Set Larva as the reference level ###
ddsHTSeq$stage <- relevel(ddsHTSeq$stage, 
                          ref = "Larva")


### Filter out all genes with <10 reads total across all samples ###
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


### Run differential expression analysis ###
dds <- DESeq(ddsHTSeq)


### Shrink LFC to obtain more accurate effect sizes ###
res.lfcShrink <- lfcShrink(dds, 
                           coef="stage_Adult_vs_Larva", 
                           type="apeglm") 


### Get summary of results ###
summary(res.lfcShrink, 
        alpha = 0.05)


### Get total number of differentially expressed genes ###
sum(res.lfcShrink$padj < 0.05, 
    na.rm=TRUE)


### Add functional annotations to results ###
annot <- read_csv("Bicyclus_anynana_v1_2_genes_only_b2g_functional_annotation_B2G-Desc_and_Best_Hit_3.csv")
res.lfcShrink$SeqName <- row.names(res.lfcShrink)
res.lfcShrink <- as.data.frame(res.lfcShrink) %>% 
  right_join(annot,
             ., 
             by="SeqName")


### Order all results by FDR ###
res.lfcShrink.ordered <- res.lfcShrink[order(res.lfcShrink$padj),]


### Extract all differentially expressed genes ###
res.lfcShrink_05 <- res.lfcShrink.ordered %>% 
  filter(padj < 0.05)


### Order all differentially expressed genes by effect size (the absolute value of log2FoldChange) ###
res.lfcShrink_05_lfc <- res.lfcShrink_05[order(-abs(res.lfcShrink_05$log2FoldChange)),]


### Extract genes upregulated in adults relative to larvae ###
res.lfcShrink_05_adult_up <- res.lfcShrink_05 %>% 
  filter(log2FoldChange > 0)


### Extract genes upregulated in larvae relative to adults ###
res.lfcShrink_05_larva_up <- res.lfcShrink_05 %>% 
  filter(log2FoldChange < 0)
```

## Make a volcano plot of the expression data:

```
### Create dataframe with color column denoting differentially expressed genes (red =  FDR < 0.05; black = FDR < 0.05) ###
dat <- res.lfcShrink.ordered %>% 
  rownames_to_column() %>% 
  mutate(color = ifelse(res.lfcShrink.ordered$padj < 0.05, 
                        "red2", 
                        "black"))


### Convert any padj=0 instances to smallest possible number for plotting (avoids Inf values; -log10(0)=Inf) ###
dat$padj <- ifelse(dat$padj == 0, 
                    .Machine$double.xmin, 
                    dat$padj)


### Make volcano plot ###
vol <- ggplot(dat, 
              aes(x = log2FoldChange, 
                  y = -log10(padj), 
                  color=color)) +
  geom_point() + 
  scale_color_identity() +
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  labs(x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR")) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())  +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  coord_cartesian(xlim = c(-17.5, 17.5), 
                  ylim = c(-10, -log10(min(dat$padj, na.rm = T))*1.05),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(-15, 15, by = 5)) +
  scale_y_continuous(breaks = seq(0, 300, by = 50)) 


### Highlight differentially expressed genes of interest ###
gene_names <- read_csv("Genes_of_interest.csv") %>%
  select(SeqName)
gene_names2 <- subset(dat, dat$rowname %in% gene_names) %>% 
  filter(padj < 0.05)
gene_names2$color <- "dodgerblue"
vol2 <- vol + 
        geom_point(gene_names2, 
                   mapping = aes(x = log2FoldChange, 
                                 y = -log10(padj), 
                                 color=color))
```
