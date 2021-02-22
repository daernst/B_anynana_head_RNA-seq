## Create index for *B. anynana* genome and align reads:

```
#!/bin/bash

### ADD SCHEDULER HEADER HERE ###


### Load required modules ###
module load cufflinks


### Convert genome gff annotation file to gtf using Cufflinks gffread ###
gunzip /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2.gff3.gz

/share/apps/bioinformatics/cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread \
    /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2.gff3 \
    -T \
    -o /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2.gtf


### Create genome index ###
/home/daernst/Bioinformatics_programs/STAR-2.7.1a/bin/Linux_x86_64_static/STAR \
    --runThreadN 2 \
    --runMode genomeGenerate \
    --genomeDir /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome \
    --genomeFastaFiles /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2_-_scaffolds.fa \
    --sjdbGTFfile /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2.gtf \
    --sjdbOverhang 49


### Align reads to genome ###

for SAMPLE in F1-Ad-2_S9 F1-Ad-6_S6 F1-L-1_S1 F2-Ad-1_S2 F2-Ad-2_S11 F2-L-1_S5 F3-Ad-3_S4 F3-L-2_S12 F3-L-7_S7 F4-Ad-1_S8 F4-L-1_S10 F4-L-4_S3

do

    /home/daernst/Bioinformatics_programs/STAR-2.7.1a/bin/Linux_x86_64_static/STAR \
        --runThreadN 12 \
        --runMode alignReads \
        --genomeDir /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome \
        --readFilesCommand zcat \
        --readFilesIn /home/daernst/B_anynana_Larva_Adult_Head_Study/adapter-trimmed_fastq_files/EW-DE-12S-${SAMPLE}_L003_R1_001_adapt-trimmed.fastq.gz \
        --outFileNamePrefix /home/daernst/B_anynana_Larva_Adult_Head_Study/STAR_alignment_output/${SAMPLE}_trimmed_ \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic

done
```
