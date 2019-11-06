* [基于shell的demo](#基于shell的demo)
* [1)、准备条件](#1准备条件)
    * [①、HPC集群](#hpc集群)
    * [②、脚本及配置文件](#脚本及配置文件)
* [2)、脚本文件](#2脚本文件)
* [3)、运行方法](#3运行方法)

***

### 基于shell的demo
补充材料，可自行测试
#### 1)、准备条件
##### ①、HPC集群
```shell
pcluster create XXX #创建SGE集群
```

##### ②、脚本及配置文件
#### 2)、脚本文件
保存以下代码为指定文件名，需要与后续运行命令相匹配。
FileName：***command.sh***
```shell
#bwa

mkdir -p $1/01.bwa
cd $1/01.bwa
bwa mem -t 4 -R "@RG\tID:group_1\tLB:library_1\tPL:illumina\tSM:sample_1" /genomes/reference/hg38/v0/Homo_sapiens_assembly38.fasta /genomes/testdata/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq /genomes/testdata/wgs_fastq/NA12878_20k/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq > $1/01.bwa/sample.sam


#samtools
samtools view -b $1/01.bwa/sample.sam > $1/01.bwa/sample.bam
samtools sort $1/01.bwa/sample.bam > $1/01.bwa/sample_sort.bam
samtools index $1/01.bwa/sample_sort.bam


#gatk
mkdir -p $1/02.gatk
cd $1/02.gatk
gatk HaplotypeCaller -R /genomes/reference/hg38/v0/Homo_sapiens_assembly38.fasta -I $1/01.bwa/sample_sort.bam -O $1/02.gatk/sample.vcf
gatk HaplotypeCaller -R /genomes/reference/hg38/v0/Homo_sapiens_assembly38.fasta -I $1/01.bwa/sample_sort.bam -O $1/02.gatk/sample.g.vcf.gz -L  /genomes/testdata/intervals/wgs_calling_regions.hg38.interval_list -ERC GVCF -dbsnp /genomes/reference/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
gatk MergeVcfs --INPUT=$1/02.gatk/sample.g.vcf.gz --OUTPUT=$1/02.gatk/sample.vcf

```

FileName：***run.sh***
```shell
sh command.sh /genomes/temp/1
```

#### 3)、运行方法
```shell
#pbs
echo "sh run.sh" | qsub  -l nodes=4,walltime=2:00:00,mem=8gb -q batch

#sge
echo "sh run.sh" | qsub -l vf=8G,s_core=4 -q all.q

#slurm
sbatch -n 4 run.sh
```
