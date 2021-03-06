* [基于workflow工具调度的WES分析](#基于workflow工具调度的wes分析)
* [1)、准备条件](#1准备条件)
    * [①、HPC集群](#hpc集群)
    * [②、nextflow调度软件](#nextflow调度软件)
    * [③、脚本及配置文件](#脚本及配置文件)
* [2)、脚本文件](#2脚本文件)
* [3)、运行方法](#3运行方法)

***

### 基于workflow工具调度的WES分析
补充材料，可自行测试
#### 1)、准备条件
##### ①、HPC集群
```shell
pcluster create XXX #创建SGE集群
```

##### ②、nextflow调度软件
```shell
curl -s https://get.nextflow.io | bash
./nextflow run hello  #测试是否正常
```

##### ③、脚本及配置文件
#### 2)、脚本文件
保存以下代码为指定文件名，需要与后续运行命令相匹配。
FileName：***custom.conf***
```nextflow
process{
    executor="sge"
    clusterOptions='-S /bin/bash'
}
aws {
  accessKey = 'xxx'
  secretKey = 'xxx'
  region = 'cn-north-1'
}
```

FileName：***genome.nf***:

```nextflow
#!/usr/bin/env nextflow
/* wget -qO- https://get.nextflow.io | bash*/

VERSION="0.2"

log.info "===================================================================="
log.info "GATK4 Best Practice Nextflow Pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}


fastq1 = file("$params.fastq1")
fastq2 = file("$params.fastq2")
params.outdir = "./Results"
params.samplename = fastq1.baseName
params.rg = fastq1.baseName
reference = file("/genomes/reference/hg19/v0/Homo_sapiens_assembly19.fasta")
dbsnp = file("/genomes/reference/hg19/v0/Homo_sapiens_assembly19.dbsnp138.vcf")
golden_indel = file("/genomes/reference/hg19/v0/Homo_sapiens_assembly19.known_indels_20120518.vcf")
hapmap = file("/genomes/reference/hg19/v0/hapmap_3.3.b37.vcf.gz")
omni = file("/genomes/reference/hg19/v0/1000G_omni2.5.b37.vcf.gz")

process BWA {
    publishDir "${params.outdir}/MappedRead"

    output:
    file 'aln-pe.sam' into samfile
    
    """
    bwa mem -M -t 8 -R '@RG\\tID:${params.rg}\\tSM:${params.samplename}\\tPL:Illumina' $reference $fastq1 $fastq2 > aln-pe.sam
    """
        
}

process BWA_sort {
    publishDir "${params.outdir}/MappedRead"
    
    input:
    file samfile

    output:
    file 'aln-pe-sorted.bam' into bam_sort

    """
    samtools sort -@ 8 -o aln-pe-sorted.bam -O BAM $samfile
    """

}

process MarkDuplicates {
    publishDir "${params.outdir}/MappedRead"
    
    input:
    file bam_sort

    output:
    file 'aln-pe_MarkDup.bam' into bam_markdup

    """
    /Pipeline/01.software/gatk/build/install/gatk/bin/gatk MarkDuplicates -I $bam_sort -M metrics.txt -O aln-pe_MarkDup.bam    
    """

}

process BaseRecalibrator {
    publishDir "${params.outdir}/BaseRecalibrator"
    
    input:
    file bam_markdup

    output:
    file 'recal_data.table' into BaseRecalibrator_table

    """
    /Pipeline/01.software/gatk/build/install/gatk/bin/gatk BaseRecalibrator \
    -I $bam_markdup \
    --known-sites $dbsnp \
    --known-sites $golden_indel \
    -O recal_data.table \
    -R $reference
    """
}

process ApplyBQSR {
    publishDir "${params.outdir}/BaseRecalibrator"
    
    input:
    file BaseRecalibrator_table
    file bam_markdup

    output:
    file 'aln-pe_bqsr.bam' into bam_bqsr
    
    script:
    """
    /Pipeline/01.software/gatk/build/install/gatk/bin/gatk ApplyBQSR -I $bam_markdup -bqsr $BaseRecalibrator_table -O aln-pe_bqsr.bam
    """
}

process HaplotypeCaller {
    publishDir "${params.outdir}/HaplotypeCaller"
    
    input:
    file bam_bqsr

    output:
    file 'haplotypecaller.g.vcf' into haplotypecaller_gvcf
    
    script:
    """
    /Pipeline/01.software/gatk/build/install/gatk/bin/gatk HaplotypeCaller -I $bam_bqsr -O haplotypecaller.g.vcf --emit-ref-confidence GVCF -R $reference
    """
}

process GenotypeGVCFs {
    publishDir "${params.outdir}/HaplotypeCaller"
    
    input:
    file haplotypecaller_gvcf

    output:
    file 'haplotypecaller.vcf' into haplotypecaller_vcf
    
    script:
    """
    /Pipeline/01.software/gatk/build/install/gatk/bin/gatk GenotypeGVCFs  --dbsnp $dbsnp --variant haplotypecaller.g.vcf -R $reference -O haplotypecaller.vcf
    """
}

process finish {
    publishDir "${params.outdir}/Report", mode: 'copy'

    input:
    file haplotypecaller_vcf

    output:
    file "${params.samplename}.vcf"
    
    """
    mv $haplotypecaller_vcf ${params.samplename}.vcf
    """

}
```
#### 3)、运行方法
```shell
#nextflow run -c ***[custom.conf]*** ***[nextflow script]*** --fastq1 ***[fastq1]*** --fastq2 ***[fastq2] -with-report xxx -with-timeline xxx -with-dag xxx ***-resume******
nextflow run -c custom.conf genome.nf --fastq1 /genomes/project/nf/SRR622461_1.fastq.gz --fastq2 /genomes/project/nf/SRR622461_1.fastq.gz ***-with-report xxx -with-timeline xxx -with-dag xxx ***-resume******
```
