## Run BLASTX to annotate the *B. anynana* reference genome (v1.2)

#### This first code chunk splits a fasta containing all *B. anynana* genes into individual fasta files:

```
#!/bin/bash

### Load required modules ###
module load perl


### Split gene sequences into individual fasta files ###
./fasta-splitter.pl \
    --part-size 1 \
    --measure count \
    --line-length 0 \
    --nopad \
    --out-dir /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2_genes_only.fa_split/single_seqs \
    /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2_genes_only.fa
```

#### The following code submits a blast job array, submitting batches of 50 sequences per job:

```
#!/bin/bash

### ADD SCHEDULER HEADER HERE ### 


### Load required module ###
module load blast/2.6.0+


### Setup array job to submit batch of 50 fasta files ###
e=$(( $PBS_ARRAYID * 50 ))
s=$(( $e - 49 ))


### Blast sequences ###
for ((f=$s;f<=$e;f++))

do

    echo "Blasting sequence $f"
    blastx \
        -query /home/daernst/B_anynana_Larva_Adult_Head_Study/B_anynana_genome/Bicyclus_anynana_v1.2_genes_only.fa_split/single_seqs/Bicyclus_anynana_v1.2_genes_only.part-${f}.fa \
        -out /home/daernst/B_anynana_Larva_Adult_Head_Study/Blastx/XML_files/single_seqs_Razor/Bicyclus_anynana_v1.2_genes_only.part-${f}.xml \
        -db /scratch/nr/nr \
        -evalue 1e-3 \
        -outfmt 5 \
        -num_threads 8 \
        -max_target_seqs 10
    echo "Finished sequence $f"

done
```
