##########################################################################
#This script takes a directory of PE fastq RNA-seq reads,
#1. Fastq quality control using [fastqc]
#2. Trims adapters using [Trmmomatic]
#3. Mapping using [Hisat2/STAR]
#4. Count using [Featurecount/stringtie]

#Created by Shuo Wang on 2023-07-17.
#Issue report on mywanuo@163.com
#All rights reserved.
#First RUN in biotrain.vip server
##########################################################################

mkdir nhzy9 && cd nhzy9

vim downlist.txt
#=======================================================================================
#                                       1.数据下载
#=======================================================================================
mkdir rawdata && cd rawdata

#=======================================================================================
#								2.fastqc质控
#=======================================================================================

mkdir fastqc_out
vim fastqc.sh
#!/bin/bash
fastqc -t 6 *.gz -o ./fastqc_out/ > ./fastqc_out/fastqc.out.log 2>&1
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N fz -cwd -l vf=50G,p=15 -q all.q@icloud-lnode01 fastqc.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
multiqc -n fastqc.before.multi ./

#=======================================================================================
#								3.trim_galore
#=======================================================================================
mkdir 3.trim.out && cd 3.trim.out

ln -s ../rawdata/*fq.gz ./

vim trim.p.sh
#!/bin/bash
time parallel -j 10 --xapply \
'trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --fastqc --paired \
-o /home/data/ssy310/workspace/reproduction/nhzy9/3.trim.out {1}_1.fq.gz {1}_2.fq.gz' \
::: `ls *_1.fq.gz | cut -d "_" -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N trim -cwd -l vf=20G,p=100 -q all.q@icloud-hnode01 trim.p.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
multiqc -n fastqc.after.multi ./

#=======================================================================================
							###4.STAR
#=======================================================================================

# STAR version=STAR_2.5.3a
#构建ensembl index
vim starindex.sh
#!/bin/bash
#!/bin/bash
STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir /home/data/ssy310/db/gencode/star/mouseGRCm39 \
--genomeFastaFiles /home/data/ssy310/db/gencode/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile /home/data/ssy310/db/gencode/gencode.vM32.annotation.gtf \
--sjdbOverhang 149
#--runThreadN 12使用12核
#--runMode genomeGenerate选genomeGenerate时是构建index
#--genomeDir /data/home/vip467/meng/STAR_index创建index的路径
#--genomeFastaFiles /data/home/vip467/meng/STAR_index/ref_hg38.fa需要的reference位置
#--dbGTFfile /data/home/vip467/meng/reference/gtf/hg38_refseq_from_ucsc.gtf
#GTF文件的位置
#--sjdbOverhang 149 如果分析的read length长度是100就写100，长度是150就写149，本次数据是100，目前的illumina都是150
#运行错误因为没有足够容量，实际本次两条染色体index是4.4G

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N starindex -cwd -l vf=60G,p=15 -q all.q@icloud-lnode01 starindex.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#比对mapping 
mkdir 4.star && cd 4.star

ln -s /home/data/ssy310/workspace/reproduction/nhzy9/3.trim.out/*val_2.fq.gz /home/data/ssy310/workspace/reproduction/nhzy9/4.star
ln -s /home/data/ssy310/workspace/reproduction/nhzy9/3.trim.out/*val_1.fq.gz /home/data/ssy310/workspace/reproduction/nhzy9/4.star

vim star_gencode_align.sh
#!/bin/bash

parallel -j 5 --xapply \
'STAR --genomeDir /home/data/ssy310/db/gencode/star/mouseGRCm39 \
    --runThreadN 4 \
    --twopassMode Basic \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outSAMattributes All \
    --outFilterIntronMotifs RemoveNoncanonical \
    --readFilesIn {1}_1_val_1.fq.gz  {1}_2_val_2.fq.gz \
    --sjdbGTFfile /home/data/ssy310/db/gencode/gencode.vM32.annotation.gtf \
--outFileNamePrefix /home/data/ssy310/workspace/reproduction/nhzy9/4.star/{1}_ > /home/data/ssy310/workspace/reproduction/nhzy9/4.star/{1}STARalign_out.log 2>&1 ' \
::: `ls *_1.fq.gz | cut -d "_" -f 1 | sort -u`

# --twopassMode Basic \#先按索引进行第一次比对，而后把第一次比对发现的新剪切位点信息加入到索引中进行第二次比对。这个参数可以保证更精准的比对情况，但是费时也费内存。
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N staralign -cwd -l vf=300G,p=80 -q all.q@icloud-hnode01 star_gencode_align.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

vim bam.sort.parallel.sh
#!/bin/bash
time parallel -j 2 --xapply \
'samtools sort -@ 8 -O bam {1}_Aligned.out.bam -o {1}_sorted.bam ; samtools index -@ 8 {1}_sorted.bam > ./{1}samtools_out.log 2>&1' \
::: `ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u`

# -G参数指定参考基因组的gtf文件，-o指定输出的文件，格式也为gtf, -b指定ballgown的输出结果目录，这个参数是为了方便下游进行ballgown差异分析，-e参数要求软件只输出已知转录本的定量结果。
# 其中 -G 参数指定基因组注释文件， -o 输出的 gtf 路径， -e 参数限定分析基因组注释文件转录本，也就是说只分析已知转录本，不分析新转录本。
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N samto -cwd -l vf=100G,p=30 -q all.q@icloud-hnode01 bam.sort.parallel.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#=======================================================================================
###							5.featureCounts RNA-Seq的定量与标准化
#=======================================================================================
mkdir 5.featureCounts && cd 5.featureCounts

ln -s ../4.star/*sorted.bam ./

#featureCounts in subread

vim featurecount.sh
#!/bin/bash
sorted_bam=`ls *sorted.bam`
for each in ${sorted_bam} ;
do
    echo $each
    name=${each##*/}
    echo $name
    SRR=${name%%_*}
    echo $SRR
    featureCounts -a /home/data/ssy310/db/gencode/gencode.vM32.annotation.gtf \
     --countReadPairs -p -T 6 \
    -o ./${SRR}_featurecounts.txt -t exon -g gene_id $each > ./${SRR}_count.log 2>&1
done

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N feat -cwd -l vf=100G,p=30 -q all.q@icloud-lnode04 featurecount.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

for n in *txt ; do cat $n | cut -f 1,7 | grep -v ^# > ${n}.readscount ; done
paste *readscount | awk '{printf $1"\t"; for(i=2;i<=NF;i+=2){printf $i"\t"};print ""}' > nhzy9.counts.txt



######################################################################################################
######################################################################################################
######################################################################################################
##########################################################################
#This script takes a directory of PE fastq RNA-seq reads,
#1. Fastq quality control using [fastqc]
#2. Trims adapters using [Trmmomatic]
#3. Mapping using [Hisat2/STAR]
#4. Count using [Featurecount/stringtie]

#Created by Shuo Wang on 2024-01-03.
#Issue report on mywanuo@163.com
#All rights reserved.
#First RUN in zhw
##########################################################################

mkdir nhzy9 && cd nhzy9

#=======================================================================================
#                                       1.数据下载
#=======================================================================================
mkdir download && cd download

# 下载诺和致源下载工具lnd
# 然后
./linuxnd/lnd login -u xxx -p xxx
./linuxnd/lnd list oss://目录名称 :列举目录下的所有文件
./linuxnd/lnd cp oss:// 目录/文件 本地目录 : 下载文件 到 本地, 下载目录命令不能下载根目录
./linuxnd/lnd cp -d oss://CP2022012100093 一个目录 本地目录 : 下载一个目录到本地


#=======================================================================================
#								2.fastqc质控
#=======================================================================================

mkdir fastqc_out && cd fastqc_out

ln -s ../download/*fq.gz ./

vim fastqc.sh
#!/bin/bash
fastqc -t 4 *.gz -o ./ 
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N fz -cwd -l vf=50G,p=15 -q all.q@icloud-hnode01 fastqc.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
multiqc -n fastqc.before.multi ./


#=======================================================================================
#								3.trim_galore
#=======================================================================================
mkdir trim.out && cd trim.out

ln -s ../download/*fq.gz ./

vim trim.p.sh
#!/bin/bash
time parallel -j 20 --xapply \
'trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --fastqc --paired \
-o /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out {1}_1.fq.gz {1}_2.fq.gz' \
::: `ls *_1.fq.gz | cut -d "_" -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N trim -cwd -l vf=20G,p=40 -q all.q@icloud-lnode03 trim.p.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
multiqc -n fastqc.after.multi ./

#=======================================================================================
#                                   4.Hisat2
#=======================================================================================
# hisat2版本2.2.1
mkdir mm39_104_hisat2

# https://mp.weixin.qq.com/s/rWk_Su9nX9w3VzBIfhJhjQ
# 提取剪接位点和外显子信息
extract_splice_sites.py Mus_musculus.GRCm39.104.gtf > Mus_musculus.GRCm39.104.ss
extract_exons.py Mus_musculus.GRCm39.104.gtf > Mus_musculus.GRCm39.104.exon
# 建立索引
# 官网直接下载http://daehwankimlab.github.io/hisat2/download/#m-musculus
# hisat2-build --ss Mus_musculus.ss \
#                --exon Mus_musculus.exon \
#                Mus_musculus.GRCm39.dna.primary_assembly.fa \
#                Mus_musculus.GRCm39_tran

vim hisat2index.mm39.sh
#!/bin/bash
fa=/home/wangshuo/db/ensembl/genome/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa
fa_baseName=GRCm39.104.dna
index_address=/home/wangshuo/db/ensembl/index/hisat2/hisat2.2.1/mm39_104_hisat2
gtf_address=/home/wangshuo/db/ensembl/genome/GRCm39
hisat2-build -p 12 -f ${fa} \
    --ss ${gtf_address}/Mus_musculus.GRCm39.104.ss \
    --exon ${gtf_address}/Mus_musculus.GRCm39.104.exon \
    ${index_address}/${fa_baseName}

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hsmm39 -cwd -l vf=50G,p=15 -q all.q@icloud-lnode05 hisat2index.mm39.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# onlyfa
vim hisat2index.mm39onlyfa.sh
#!/bin/bash
fa=/home/wangshuo/db/ensembl/genome/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa
fa_baseName=GRCm39.104.dna
index_address=/home/wangshuo/db/ensembl/index/hisat2/hisat2.2.1/mm39_104_hisat2onlyfa
hisat2-build -p 12 -f ${fa} ${index_address}/${fa_baseName}

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hsmm39onlyfa -cwd -l vf=50G,p=15 -q all.q@icloud-lnode06 hisat2index.mm39onlyfa.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 多个样本批量进行比对，排序，建索引
# Hisat.sh内容： 注意命令中的-，表示占位符，表示|管道符前面的输出。
## 此处索引直接使用服务器上已经构建好的进行练习

mkdir hisat2  && cd hisat2
vim hisat2.sh

#!/bin/bash
index=/home/wangshuo/db/ensembl/index/hisat2/hisat2.2.1/mm39_104_hisat2/GRCm39.109.dna
inputdir=/data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out
outdir=/data3/Group2/wangshuo3/rnaseq/prjna449550hfsex/4.hisat2

ls /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/*.gz |cut -d "_" -f 1 | sort -u | while read id
do
    ls ${id}_1_val_1.fq.gz  ${id}_2_val_2.fq.gz;
    SRR=${id##*/};
    echo $SRR;
    hisat2 -p 6 -x ${index} -1 ${inputdir}/${SRR}_1_val_1.fq.gz -2 ${inputdir}/${SRR}_2_val_2.fq.gz 2>${SRR}.log  | samtools sort -@ 6 -o ${outdir}/${SRR}.Hisat_aln.sorted.bam - &&  samtools index ${outdir}/${SRR}.Hisat_aln.sorted.bam ${outdir}/${SRR}.Hisat_aln.sorted.bam.bai
done
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hfsali -cwd -l vf=50G,p=15 -q all.q@icloud-mnode01 hisat2.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#基于--ss --exon建立index
ln -s ../trim.out/*_1_val_1.fq.gz ./
vim hisat2parallel.sh
#!/bin/bash
time parallel -j 10 --xapply \
'hisat2 -p 6 \
-x /home/wangshuo/db/ensembl/index/hisat2/hisat2.2.1/mm39_104_hisat2/GRCm39.104.dna \
--summary-file /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2/{1}_summary.txt \
-1 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_1_val_1.fq.gz \
-2 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_2_val_2.fq.gz  | \
samtools sort -@ 6 -o /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2/{1}.Hisat_aln.sorted.bam - &&  \
samtools index /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2/{1}.Hisat_aln.sorted.bam \
/data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2/{1}.Hisat_aln.sorted.bam.bai ' \
::: `ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hfsali -cwd -l vf=50G,p=30 -q all.q@icloud-hnode01 hisat2parallel.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 【参考】https://mp.weixin.qq.com/s/rWk_Su9nX9w3VzBIfhJhjQ
# 创建一个脚本文件然后输入批量比对的脚本
# $ vi hisat2_map.sh
# #!/bin/bash
# batch map to genome
# for i in test1 test2
# do
#     for j in rep1 rep2
#     do
#         hisat2 -p 10 --dta -x 2.index/grcm38_tran/genome_tran \
#                --summary-file 3.map-data/${i}_${j}_summary.txt \
#                -1 1.fastq-data/${i}_R1_${j}.fq.gz \
#                -2 1.fastq-data/${i}_R2_${j}.fq.gz \
#                | samtools sort -@ 10 -o 3.map-data/${i}_${j}.sorted.bam
#     done
# done
# 保存后挂后台运行
# $ nohup ./map.sh &

ln -s ../trim.out/*_1_val_1.fq.gz ./
vim hisat2parallel.sh
#!/bin/bash
time parallel -j 10 --xapply \
'hisat2 -p 6 --dta \
-x /home/wangshuo/db/ensembl/index/hisat2/hisat2.2.1/mm39_104_hisat2onlyfa/GRCm39.104.dna \
--summary-file /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa/{1}_summary.txt \
-1 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_1_val_1.fq.gz \
-2 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_2_val_2.fq.gz  | \
samtools sort -@ 6 -o /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa/{1}.Hisat2_aln.sorted.bam - &&  \
samtools index /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa/{1}.Hisat2_aln.sorted.bam \
/data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa/{1}.Hisat2_aln.sorted.bam.bai ' \
::: `ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hsof -cwd -l vf=50G,p=30 -q all.q@icloud-hnode01 hisat2parallel.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# -x跟索引名前缀，-1，-2跟双端测序文件，-U跟单端测序文件，-S输出为sam格式的文件,-p线程数量
# 我们直接输出为排序好的bam文件
# 单个样本比对，--dta输出为转录本组装的reads，-@为samtools的线程数，--summary-file输出比对信息




#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 官网下载index
mkdir hisat2official  && cd hisat2official

ln -s ../trim.out/*_1_val_1.fq.gz ./
vim hisat2parallel.sh
#!/bin/bash
time parallel -j 10 --xapply \
'hisat2 -p 6 --dta \
-x /home/wangshuo/db/ensembl/index/hisat2/hisat2offical/grcm38_tran/genome_tran \
--summary-file /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2offical/{1}_summary.txt \
-1 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_1_val_1.fq.gz \
-2 /data3/Group2/wangshuo3/rnaseq/nhzy9/trim.out/{1}_2_val_2.fq.gz  | \
samtools sort -@ 6 -o /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2offical/{1}.Hisat2_aln.sorted.bam - &&  \
samtools index /data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2offical/{1}.Hisat2_aln.sorted.bam \
/data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2offical/{1}.Hisat2_aln.sorted.bam.bai ' \
::: `ls *_1_val_1.fq.gz | cut -d "_" -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hsof -cwd -l vf=50G,p=30 -q all.q@icloud-hnode02 hisat2parallel.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#=======================================================================================
###							5.stringtie RNA-Seq的定量与标准化
#=======================================================================================
mkdir hisat_stringtie && cd hisat_stringtie

ln -s ../hisat2/*sorted.bam ./

vim stringtie.sh
#!/bin/bash
sorted_bam=`ls *sorted.bam`
for each in ${sorted_bam} ;
do
    echo $each
    name=${each##*/}
    echo $name
    SRR=${name%%_*}
    echo $SRR
    stringtie  -p 6 -e -B \
    -G /home/wangshuo/db/ensembl/genome/GRCm39/Mus_musculus.GRCm39.104.gtf \
    -A ./${SRR}_stringtie_abundance.txt \
    -o ./${SRR}_stringtie.gtf $each
done


# -G参数指定参考基因组的gtf文件，-o指定输出的文件，格式也为gtf, -b指定ballgown的输出结果目录，这个参数是为了方便下游进行ballgown差异分析，-e参数要求软件只输出已知转录本的定量结果。
# 其中 -G 参数指定基因组注释文件， -o 输出的 gtf 路径， -e 参数限定分析基因组注释文件转录本，也就是说只分析已知转录本，不分析新转录本。
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N hsc -cwd -l vf=50G,p=30 -q all.q@icloud-lnode06 stringtie.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


find . -name "*.gtf"
ls *.gtf >  mergelist_simple.txt # 将所有组装的转录本文件名合并到一个文
# 合并转录本
# 建立一个 gtf 的 list 文件，里面为上一步输出的文件的路径
ls -l ./*.gtf | awk '{print $9}' > mergelist.txt

vim stringtiemerge.sh
#!/bin/bash
annotation_file=/home/wangshuo/db/ensembl/genome/GRCm39/Mus_musculus.GRCm39.104.gtf
stringtie --merge  -p 8 -G $annotation_file -o nhzy9_hisat2onlyfa_stringtie_merged.gtf mergelist.txt
#这一步是用--merge指令将所有转录本合并输出到merge.gtf文件中

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N strmerge -cwd -l vf=50G,p=30 -q all.q@icloud-lnode05 stringtiemerge.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#########################################################################################################
# 基于merge的定量
#########################################################################################################

# 建立定量文件夹
mkdir quantity && cd quantity
# 定量，单个样本，-e评估转录本表达丰度，-A评估基因表达丰度并输出，-G跟合并的gtf文件
# $ stringtie  -p 10 -e -G ./stringtie_merged.gtf \
#              -o quantity/test1_rep1/test1_rep1.gtf \
#              -A 5.quantity/test1_rep1/gene_abundances.tsv \
#              3.map-data/test1_rep1.sorted.bam

vim stringtiequantity.sh
#!/bin/bash
sorted_bam=`ls *sorted.bam`
for each in ${sorted_bam} ;
do
    echo $each
    name=${each##*/}
    echo $name
    SRR=${name%%.*}
    echo $SRR
    stringtie  -p 10 -e \
    -G ./nhzy9_hisat2onlyfa_stringtie_merged.gtf \
    -o quantity/${SRR}/${SRR}_quantity.gtf \
    -A quantity/${SRR}/${SRR}gene_abundances.tsv $each
done

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N sq -cwd -l vf=50G,p=30 -q all.q@icloud-lnode05 stringtiequantity.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# sample_list.txt文件格式如下所示：两个处理的名称（s1，s2）+对应gtf文件的地址。注意第二行末尾不要有换行符号（\n）！！！
# 处理名称（s1，s2）和gtf文件地址之间应该有一个空格。
sample_list.txt
# s1 /mnt/1223/NGS2021/differential_expression/s1_out.gtf
# s2 /mnt/1223/NGS2021/differential_expression/s2_out.gtf
#...
ls *.Hisat2_stringtie.gtf | awk -F '\\.' '{print $1, "/data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa_stringtie/quantity/" $1 "/" $1 "_quantity.gtf"}' > sample_list.txt

# 需要使用一个 prepDE.py 脚本，在 python2 环境中使用，下载地址如下：
# prepDE.py (python2) : http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
# prepDE.py (python3) : http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
# prepDE.py 需要一个 sample_list，第一列为样本名，第二列为 gtf 文件路径

python prepDE.py3 -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv
#-i # 输入文件，就是前面做的sample_list.txt
#-g # 自定义基因组表达矩阵名字，默认也是gene_count_matrix.csv
#-t # 自定义转录本表达矩阵名字，默认也是transcript_count_matrix.csv
#其中 -l 是设置平均 reads 长度。prepDE.py 计算 read counts 代码如下，read_len 就是 -l 参数。
#可见这个参数需要按实际实验去设置，默认值是75.
# 运行后得到 gene_count_matrix.csv 和 transcript_count_matrix.csv 两个文件，
# 其中 gene_count_matrix.csv 就是需要的基因 read counts 矩阵。查看一下。


#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 提取 FPKM/TPM 或 coverage 结果
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 获取提取脚本：https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/stringtie_expression_matrix.pl



#########################################################################################################
# 基于不需要发现新的转录本和基因,直接定量
#########################################################################################################
# 如果不需要发现新的转录本和基因的话，可直接基于 bam 文件走定量步骤

# 1、HTseq 定量（使用参考：使用 htseq-count 进行定量分析[11]）
# 2、Subread 包中的 featureCounts 定量（使用参考：转录组定量-featureCounts[12]）
# 3、stringtie 定量
# [这里直接使用 stringtie 定量]

# 创建一个脚本文件然后输入批量定量的脚本
# $ vi quantity.sh
# #!/bin/bash
# quantity
# for i in test1 test2
# do
#         for j in rep1 rep2
#         do
#             stringtie  -p 10 -e -G ./Mus_musculus.GRCm38.102.gtf \
#                         -o 5.quantity-data/${i}_${j}/${i}_${j}.gtf \ #该处${i}_${j}文件夹下必须命名为${i}_${j}.gtf
#                         -A 5.quantity-data/${i}_${j}/gene_abundances.tsv \ #该处${i}_${j}文件夹下必须命名为gene_abundances.tsv
#                         3.map-data/${i}_${j}.sorted.bam
#         done
# done

# 建立定量文件夹
mkdir quantitybam && cd quantitybam

vim stringbamtiequantityparallel.sh
#!/bin/bash
time parallel -j 10 --xapply \
'stringtie -p 4 -e \
    -G /home/wangshuo/db/ensembl/genome/GRCm39/Mus_musculus.GRCm39.104.gtf \
    -o quantitybam/{1}/transcripts.gtf \
    -A quantitybam/{1}/gene_abundances.tsv {1}.Hisat_aln.sorted.bam ' \
::: `ls *sorted.bam | cut -d "." -f 1 | sort -u`

#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
conda activate dnaws
qsub -V -N sq -cwd -l vf=50G,p=30 -q all.q@icloud-lnode06 stringbamtiequantityparallel.sh
#-------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# sample_list.txt文件格式如下所示：两个处理的名称（s1，s2）+对应gtf文件的地址。注意第二行末尾不要有换行符号（\n）！！！
# 处理名称（s1，s2）和gtf文件地址之间应该有一个空格。
sample_list.txt
# s1 /mnt/1223/NGS2021/differential_expression/s1_out.gtf
# s2 /mnt/1223/NGS2021/differential_expression/s2_out.gtf
#...
ls *.Hisat2_stringtie.gtf | awk -F '\\.' '{print $1, "/data3/Group2/wangshuo3/rnaseq/nhzy9/hisat2onlyfa_stringtie/quantitybam/" $1 "/transcripts.gtf"}' > sample_list.txt

# 需要使用一个 prepDE.py 脚本，在 python2 环境中使用，下载地址如下：
# prepDE.py (python2) : http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
# prepDE.py (python3) : http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
# prepDE.py 需要一个 sample_list，第一列为样本名，第二列为 gtf 文件路径

python prepDE.py3 -i sample_list.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv
#-i # 输入文件，就是前面做的sample_list.txt
#-g # 自定义基因组表达矩阵名字，默认也是gene_count_matrix.csv
#-t # 自定义转录本表达矩阵名字，默认也是transcript_count_matrix.csv
#其中 -l 是设置平均 reads 长度。prepDE.py 计算 read counts 代码如下，read_len 就是 -l 参数。
#可见这个参数需要按实际实验去设置，默认值是75.
# 运行后得到 gene_count_matrix.csv 和 transcript_count_matrix.csv 两个文件，
# 其中 gene_count_matrix.csv 就是需要的基因 read counts 矩阵。查看一下。



#########################################################################################################
# 提取 FPKM/TPM 或 coverage 结果
#########################################################################################################
# 获取提取脚本：https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/stringtie_expression_matrix.pl

# 查看用法，--result_dirs跟上含有每个样本gtf的文件夹名称
./stringtie_expression_matrix.pl
Required parameters missing
Usage:  ./stringtie_expression_matrix.pl --expression_metric=TPM  --result_dirs='HBR_Rep1,HBR_Rep2,HBR_Rep3,UHR_Rep1,UHR_Rep2,UHR_Rep3' --transcript_matrix_file=transcript_tpms_all_samples.tsv --gene_matrix_file=gene_tpms_all_samples.tsv


ls *.Hisat2_stringtie.gtf | awk -F '\\.' '{print $1}' | tr '\n' ',' > sample_result_dirs.txt
cd quantitybam

# 提取TPM
./stringtie_expression_matrix.pl --expression_metric=TPM \
                                   --result_dirs='T0-1,T0-2,T0-3,T0-4,T1-1,T1-2,T1-3,T1-4,T2-1,T2-2,T2-3,T2-4,T4-1,T4-2,T4-3,T4-4,T5-1,T5-2,T5-3,T5-4' \
                                   --transcript_matrix_file=transcript_tpms_all_samples.tsv \
                                   --gene_matrix_file=gene_tpms_all_samples.tsv

# 提取FPKM
./stringtie_expression_matrix.pl --expression_metric=FPKM \
                                   --result_dirs='T0-1,T0-2,T0-3,T0-4,T1-1,T1-2,T1-3,T1-4,T2-1,T2-2,T2-3,T2-4,T4-1,T4-2,T4-3,T4-4,T5-1,T5-2,T5-3,T5-4' \
                                   --transcript_matrix_file=transcript_fpkms_all_samples.tsv \
                                   --gene_matrix_file=gene_fpkms_all_samples.tsv

# 提取coverage
./stringtie_expression_matrix.pl --expression_metric=coverage \
                                   --result_dirs='T0-1,T0-2,T0-3,T0-4,T1-1,T1-2,T1-3,T1-4,T2-1,T2-2,T2-3,T2-4,T4-1,T4-2,T4-3,T4-4,T5-1,T5-2,T5-3,T5-4' \
                                   --transcript_matrix_file=transcript_coverage_all_samples.tsv \
                                   --gene_matrix_file=gene_coverage_all_samples.tsv
# 在当前目录就会生成相应的基因和转录本的tpm、fpkm、coverage 结果


