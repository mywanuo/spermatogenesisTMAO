```shell
##########################################################################
#This script takes a directory of Metagenomics,
##
# ref:https://github.com/YongxinLiu/EasyMetagenome
##
#Created by Shuo Wang on 2024-07-10.
#Issue report on mywanuo@163.com
#Copyright (c) 2023. All rights reserved.
##
#######################################################################

## Database prepare and software install
[Database prepare and software install link](https://github.com/YongxinLiu/EasyMetagenome/blob/master/0Install.sh)

## 1.1 Prepare

# Environment variable setup (must be run before starting each analysis)

# Set database, software, and working directories
```shell
db=./db/metag
wd=./meige
mkdir -p $wd && cd $wd
```

```shell
mkdir -p seq temp result
mv metadata.txt result/metadata.txt
sed -i 's/\r//' result/metadata.txt
```

# Configuration for NCBI_nt nucleotide database
```shell
lftp ftp://download.nmdc.cn/tools/meta/NCBI
mirror nt/
```

```shell
# Taxonomic classification database
mkdir NCBI_tax
cd NCBI_tax
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```

```shell
fastqc --version # FastQC v0.11.9
```

```shell
mkdir fastqc_out && cd fastqc_out

vim fastqc.sh
#!bin/bash
fastqc ../seq/*.fq.gz -t 12

multiqc --version # version 1.14
multiqc -d seq/ -o result/qc
```

## Quality control
```shell
cd seq
vim pigz.decompress.sh
#!/bin/bash
parallel -j 15 --xapply \
'pigz -p 4 -d -c {1}.gz > {1}' \
::: `ls *.fq.gz | sed 's/\.gz$//'`

mkdir -p temp/qc
```

### Quality control and host read removal using KneadData

```shell
kneaddata --version #kneaddata v0.12.0
trimmomatic -version # 0.39
bowtie2 --version # version 2.4.5
```

#### Parallel quality control for multiple samples

```shell
vim kneaddata.sh
#!/bin/bash
time parallel -j 12 --xapply \
  'kneaddata -i1 seq/{1}.R1.fq -i2 seq/{1}.R2.fq \
  -o temp/qc -v -t 10 --max-memory 506384m --remove-intermediate-output \
  --trimmomatic /tools/Trimmomatic-0.38/ \
  --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50 ILLUMINACLIP:/tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10" \
  --trf ./envs/kneaddata/bin/ \
  --reorder --bowtie2-options "--very-sensitive --dovetail" \
  -db ./db/metag/kneaddata/human_genome/ ' \
  ::: `tail -n+2 result/metadata.txt|cut -f1`
```

#### Rename QC results, remove temporary files, and summarize statistics
```shell
rm -rf temp/qc/*contam* temp/qc/*unmatched*  temp/qc/*.fq
ls -l temp/qc/
```
```shell
    awk '{system("mv temp/qc/"$1"_1_kneaddata_paired_1.fastq temp/qc/"$1"_1.fastq")}' <(tail -n+2 result/metadata.txt)
    awk '{system("mv temp/qc/"$1"_1_kneaddata_paired_2.fastq temp/qc/"$1"_2.fastq")}' <(tail -n+2 result/metadata.txt)
    ls -l temp/qc/
```

#### Quality assessment after QC

```shell
vim fastqcTrim.sh
#!bin/bash
fastqc temp/qc/*.fastq -t 12

multiqc -d temp/qc/ -o result/qc/
```

## Read-based analysis using HUMAnN2

```shell
mkdir -p temp/concat
# Merge paired-end reads into a single file

vim concat.combine.sh
#!/bin/bash
for i in `tail -n+2 result/metadata.txt|cut -f1`;do 
zcat 1.Kneaddata_Clean/clean_data/${i}.kneaddata_paired_?.fastq.gz \
> temp/concat/${i}.fq; done
```

### HUMAnN2 computation of species and functional composition

Check whether the database configuration is correct
```shell
humann --version # humann v3.6.1
metaphlan --version # MetaPhlAn version 4.0.6
metaphlan --install
db=./db/metag/
cd ${db}
mkdir -p ${db}/humann3_databases
humann_databases --download utility_mapping full ${db}/humann3_databases
humann_databases --download chocophlan full ${db}/humann3_databases
humann_databases --download uniref uniref90_diamond ${db}/humann3_databases
mkdir -p temp/humann

wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar
wget -c http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
tar xvf mpa_vOct22_CHOCOPhlAnSGB_202212.tar
tar xvf mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
```

```shell
vim humann3.sh

#!/bin/bash
time parallel -j 10 --xapply \
  'zcat temp/qc/{1}.R1_kneaddata_paired_?.fastq.gz > temp/concat/{1}.fq && \
  humann --input temp/concat/{1}.fq  \
  --diamond ./envs/metapy37/bin/ \
  --search-mode uniref90 \
  --threads 8 \
  --metaphlan /home/wangshuo/miniconda3/envs/metapy37/bin \
  --metaphlan-options "--bowtie2db ./metaphlan/metaphlan_databases/ -t rel_ab_w_read_stats" \
  --output temp/humann/ && \
  pigz -p 8 temp/concat/{1}.fq ' \
  ::: `tail -n+2 result/metadata.txt|cut -f1`

# Link important files to the humann directory
for i in `tail -n+2 result/metadata.txt|cut -f1`;do 
    cp -f temp/humann/${i}_humann_temp/${i}_metaphlan_bugs_list.tsv temp/humann/
done    
# Delete temporary files to free up disk space
rm -rf temp/concat/* temp/humann/*_humann_temp
```

### Species composition table
```shell
mkdir -p result/metaphlan4

merge_metaphlan_tables.py temp/humann/*_metaphlan_bugs_list.tsv | \
    sed 's/_metaphlan_bugs_list//g' > result/metaphlan4/taxonomy.tsv

merge_metaphlan_tables.py temp/humann/*_metaphlan_bugs_list.tsv | \
  sed 's/_metaphlan_bugs_list//g' | tail -n+2 | sed '1 s/clade_name/ID/' | sed '2i #metaphlan4'> result/metaphlan4/taxonomy.tsv

merge_metaphlan_tables.py temp/humann/*_metaphlan_bugs_list.tsv > result/metaphlan4/merged_abundance_table.txt

# Filter species-level taxa
grep -E "s__|clade" result/metaphlan4/merged_abundance_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed 's/_metaphlan_bugs_list//g' \
> result/metaphlan4/merged_abundance_table_Species.txt

# Filter genus-level taxa
grep -E "g__|clade" result/metaphlan4/merged_abundance_table.txt \
| grep -v "s__" \
| sed "s/^.*|//g" \
| sed 's/_metaphlan_bugs_list//g' \
> result/metaphlan4/merged_abundance_table_Genus.txt
```

### Merge, normalize, and stratify functional composition
```shell
mkdir -p result/humann
humann_join_tables --input temp/humann \
  --file_name pathabundance \
  --output result/humann/pathabundance.tsv
```


### genefamilies
```shell
humann_join_tables \
-i temp/humann \
-o result/humann/genefamilies.tsv \
--file_name genefamilies
```

### pathcoverage
```shell
humann_join_tables \
-i temp/humann \
-o result/humann/pathcoverage.tsv \
--file_name pathcoverage
```

```shell
# kegg
humann_regroup_table \
--input result/humann/genefamilies_cpm.tsv \
--groups uniref90_ko \
--output result/humann/genefamilies_uniref90_ko_Function_cpm.tsv
``