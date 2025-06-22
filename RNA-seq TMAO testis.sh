##########################################################################
# RNA-seq Analysis Pipeline (Paired-End)
#
# Steps:
#   1. Quality Control: FastQC + MultiQC
#   2. Adapter Trimming: Trim Galore!
#   3. Alignment: STAR or HISAT2
#   4. Sorting and Indexing: samtools
#   5. Quantification: FeatureCounts or StringTie
#   6. Expression Matrix Generation: prepDE, stringtie_expression_matrix.pl
#
# Tool Versions:
#   - FastQC v0.11.9
#   - MultiQC v1.9
#   - Trim Galore! v0.6.7
#   - STAR v2.5.3a
#   - HISAT2 v2.2.1
#   - samtools v1.10+
#   - FeatureCounts (Subread v2.0.1)
#   - StringTie v2.1.4
#
# Created on 2023-07-17 by Shuo Wang
# Updated on 2024-01-03
# Contact: mywanuo@163.com
##########################################################################

mkdir testistmao && cd testistmao

#=======================================================================================
#                                       1. Data Preparation
#=======================================================================================
mkdir rawdata && cd rawdata
# (Put your data download commands here)
cd ..

#=======================================================================================
#                               2. FastQC Quality Control
#=======================================================================================
# FastQC v0.12.1
mkdir fastqc_out

vim fastqc.sh
#!/bin/bash
fastqc -t 6 ./rawdata/*.gz -o ./fastqc_out/ > ./fastqc_out/fastqc.out.log 2>&1
multiqc -n fastqc.before.multi ./fastqc_out

#=======================================================================================
#                               3. Trim Galore Adapter Trimming
#=======================================================================================
# version 0.6.10
mkdir trim.out && cd trim.out
ln -s ../rawdata/*fq.gz ./

vim trim.p.sh
#!/bin/bash
parallel -j 10 --xapply \
'trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --fastqc --paired \
-o ./ {1}_1.fq.gz {1}_2.fq.gz' \
::: $(ls *_1.fq.gz | cut -d "_" -f 1 | sort -u)
cd ../

#=======================================================================================
#                                   4. STAR Alignment
#=======================================================================================
# STAR version: 2.7.10b

# Index building
vim starindex.sh
#!/bin/bash
cd testistmao
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir ./ref/mouseGRCm39 \
--genomeFastaFiles ./ref/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile ./ref/gencode.vM32.annotation.gtf \
--sjdbOverhang 150

mkdir star && cd star
ln -s ../trim.out/*val_1.fq.gz ./
ln -s ../trim.out/*val_2.fq.gz ./

vim star_gencode_align.sh
#!/bin/bash
cd testistmao/star
parallel -j 5 --xapply \
'STAR --genomeDir ../ref/star_index \
--runThreadN 4 \
--twopassMode Basic \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outSAMattributes All \
--outFilterIntronMotifs RemoveNoncanonical \
--readFilesIn {1}_1_val_1.fq.gz {1}_2_val_2.fq.gz \
--sjdbGTFfile ../ref/gencode.vM32.annotation.gtf \
--outFileNamePrefix ./{1}_ > ./{1}STARalign_out.log 2>&1' \
::: $(ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u)


vim bam.sort.parallel.sh
#!/bin/bash
cd testistmao/star
parallel -j 2 --xapply \
'samtools sort -@ 8 -O bam {1}_Aligned.out.bam -o {1}_sorted.bam && \
 samtools index -@ 8 {1}_sorted.bam > {1}samtools_out.log 2>&1' \
::: $(ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u)

cd ../..

#=======================================================================================
#                      5. featureCounts Quantification (v2.0.1)
#=======================================================================================
mkdir featureCounts && cd featureCounts
ln -s ../star/*sorted.bam ./

vim featurecount.sh
#!/bin/bash
cd testistmao/5.featureCounts
sorted_bam=$(ls *sorted.bam)
for each in ${sorted_bam}; do
    name=${each##*/}
    SRR=${name%%_*}
    featureCounts -a ../ref/annotation.gtf \
     --countReadPairs -p -T 6 \
     -o ${SRR}_featurecounts.txt -t exon -g gene_id $each > ${SRR}_count.log 2>&1
done

# Generate final count matrix
for n in *featurecounts.txt; do
    cut -f 1,7 $n | grep -v ^# > ${n}.readscount
done
paste *readscount | awk '{printf $1"\\t"; for(i=2;i<=NF;i+=2){printf $i"\\t"}; print ""}' > testistmao.counts.txt







