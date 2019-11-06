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
