# B_anynana_head_RNA-seq

### This repository contains code for differential expression analysis of larval and adult *B. anynana* heads, including:

1. Aligning reads to the reference genome
2. DESeq2 analysis using the model *y* ~ *family* + *sex* + *stage* to analyze differental gene expression between stages while controlling for lineage and sex
3. DESeq2 analysis using the model *y* ~ *group* to analyze sex-specific differences between stages
4. Run BLASTX to annotate the *B. anynana* reference genome
