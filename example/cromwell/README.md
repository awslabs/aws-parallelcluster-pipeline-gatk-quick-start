* [基于cromwell工具调度的demo](#基于cromwell工具调度的demo)
* [1)、准备条件](#1准备条件)
    * [①、HPC集群](#hpc集群)
    * [②、cromwell调度软件](#nextflow调度软件)
    * [③、脚本及配置文件](#脚本及配置文件)
* [2)、脚本文件](#2脚本文件)
* [3)、运行方法](#3运行方法)

***

### 基于cromwell工具调度的demo
补充材料，可自行测试
#### 1)、准备条件
##### ①、HPC集群
```shell
pcluster create XXX #创建SGE集群
```

##### ②、cromwell调度软件
```shell
wget https://parallelcluster-gatk.s3.cn-north-1.amazonaws.com.cn/01.software/cromwell-39.jar
```

##### ③、脚本及配置文件
#### 2)、脚本文件
保存以下代码为指定文件名，需要与后续运行命令相匹配。
FileName：***test.wdl***
```wdl
## Copyright Broad Institute, 2017
##
## This WDL workflow runs HaplotypeCaller from GATK4 in GVCF mode on a single sample
## according to the GATK Best Practices (June 2016), scattered across intervals.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support
## - Successfully tested on v29
## - Does not work on versions < v23 due to output syntax
##
## IMPORTANT NOTE: HaplotypeCaller in GATK4 is still in evaluation phase and should not
## be used in production until it has been fully vetted. In the meantime, use the GATK3
## version for any production needs.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow HaplotypeCallerGvcf_GATK4 {
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File scattered_calling_intervals_list

  String gatk_docker

  String gatk_path

  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  String sample_basename = basename(input_bam, ".bam")

  String gvcf_name = sample_basename + ".g.vcf.gz"
  String gvcf_index = sample_basename + ".g.vcf.gz.tbi"

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval_list = interval_file,
        gvcf_name = gvcf_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path
    }

  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      vcf_name = gvcf_name,
      vcf_index = gvcf_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_merged_gvcf = MergeGVCFs.output_vcf
    File output_merged_gvcf_index = MergeGVCFs.output_vcf_index
  }
}

# TASK DEFINITIONS

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String gvcf_name
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File interval_list
  Int? interval_padding
  Float? contamination
  Int? max_alt_alleles

#   Int preemptible_tries
#   Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    ${gatk_path} --java-options ${java_opt} \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${gvcf_name} \
      -L ${interval_list} \
      -ip ${default=100 interval_padding} \
      -contamination ${default=0 contamination} \
      --max-alternate-alleles ${default=3 max_alt_alleles} \
      -ERC GVCF
  }
  output {
    File output_gvcf = "${gvcf_name}"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  Array [File] input_vcfs
  String vcf_name
  String vcf_index

#   Int preemptible_tries
#   Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    ${gatk_path} --java-options ${java_opt} \
      MergeVcfs \
      --INPUT=${sep=' --INPUT=' input_vcfs} \
      --OUTPUT=${vcf_name}
  }
  output {
    File output_vcf = "${vcf_name}"
    File output_vcf_index = "${vcf_index}"
  }
}
```

FileName：***sge.conf***
```conf
#include the application.conf file.
include required(classpath("application"))
webservice {
  port = 8000
}

aws {

  application-name = "cromwell"

  auths = [
    {
      name = "default"
      scheme = "default"
    }
  ]

  region = "cn-north-1"
}

engine {
  filesystems {
    s3 {
      auth = "default"
    }
  }
}

backend {
  # Switch the default backend to "SGE"
  default = "SGE"
  providers {

    # Configure the SGE backend
    SGE {
      # Use the config backend factory
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

      config {
        # Limits the number of concurrent jobs
        concurrent-job-limit = 500

        # Define runtime attributes for the SGE backend.
        # memory_gb is a special runtime attribute. See the cromwell README for more info.
        runtime-attributes = """
        Int cpu = 1
        Float? memory_gb
        String? sge_queue
        String? sge_project
        """

        # Script for submitting a job to SGE, using runtime attributess.
        # See the cromwell README for more info.
        submit = """
        qsub \
        -terse \
        -V \
        -b n \
        -N ${job_name} \
        -wd ${cwd} \
        -o ${out} \
        -e ${err} \
        -pe smp ${cpu} \
        ${"-l h_vmem=" + memory_gb + "g"} \
        ${"-q " + sge_queue} \
        ${"-P " + sge_project} \
        ${script}
        """

        # command for killing/aborting
        kill = "qdel ${job_id}"

        # Command used at restart to check if a job is alive
        check-alive = "qstat -j ${job_id}"

        # How to search the submit output for a job_id
        job-id-regex = "(\\d+)"
        filesystems {
          s3 {
            auth = "default"
          }
        }
      }
    }
  }
}
```


FileName：***input.json***
```json
{
    "##_COMMENT1": "INPUT BAM",
    "HaplotypeCallerGvcf_GATK4.input_bam": "s3://parallelcluster-gatk/99.testdata/wgs_bam/NA12878_24RG_b37/NA12878_24RG_small.b37.bam",
    "HaplotypeCallerGvcf_GATK4.input_bam_index": "s3://parallelcluster-gatk/99.testdata/wgs_bam/NA12878_24RG_b37/NA12878_24RG_small.b37.bai",
  
    "##_COMMENT2": "REFERENCE FILES",
    "HaplotypeCallerGvcf_GATK4.ref_dict": "/genomes/reference/hg19/v0/Homo_sapiens_assembly19.dict",
    "HaplotypeCallerGvcf_GATK4.ref_fasta": "/genomes/reference/hg19/v0/Homo_sapiens_assembly19.fasta",
    "HaplotypeCallerGvcf_GATK4.ref_fasta_index": "/genomes/reference/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
  
    "##_COMMENT3": "INTERVALS",
    "HaplotypeCallerGvcf_GATK4.scattered_calling_intervals_list": "/genomes/testdata/intervals/b37_wgs_scattered_calling_intervals.txt",
    "HaplotypeCallerGvcf_GATK4.HaplotypeCaller.interval_padding": 100,
  
    "##_COMMENT4": "DOCKERS",
    "HaplotypeCallerGvcf_GATK4.gatk_docker": "broadinstitute/gatk:4.0.0.0",
  
    "##_COMMENT5": "PATHS",
    "HaplotypeCallerGvcf_GATK4.gatk_path": "/Pipeline/01.software/gatk/gatk",
  
    "##_COMMENT6": "JAVA OPTIONS",
    "HaplotypeCallerGvcf_GATK4.HaplotypeCaller.java_opt": "-Xms4000m",
    "HaplotypeCallerGvcf_GATK4.MergeGVCFs.java_opt": "-Xms4000m",
  
    "##_COMMENT7": "MEMORY ALLOCATION",
    "HaplotypeCallerGvcf_GATK4.HaplotypeCaller.mem_size": "4 GB",
    "HaplotypeCallerGvcf_GATK4.MergeGVCFs.mem_size": "4 GB"
}
```

#### 3)、运行方法
```shell
#java -Dconfig.file=[*path to custom.**conf]* -jar [*cromwell jar package]* run [*wdl script]* -i *[input json]*
#example:
java -Dconfig.file=./conf/sge.conf -jar /genomes/project/wdl/simple-pipeline-test/cromwell-39.jar run /genomes/project/wdl/simple-pipeline-test/test.wdl -i /genomes/project/wdl/simple-pipeline-test/input.json
```


#### 4)、结果输出

测试标志：控制台以json格式输出结果文件路径，如：
```
#...other output
 [2019-05-31 06:41:55,48] [info] WorkflowExecutionActor-621bcdec-19a9-411d-8670-3899fa9e8509 [621bcdec]: Workflow HaplotypeCallerGvcf_GATK4 complete. Final Outputs:
 {
  "HaplotypeCallerGvcf_GATK4.output_merged_gvcf_index": "/genomes/project/wdl/simple-pipeline-test/cromwell-executions/HaplotypeCallerGvcf_GATK4/621bcdec-19a9-411d-8670-3899fa9e8509/call-MergeGVCFs/execution/NA12878_24RG_small.b37.g.vcf.gz.tbi",
  "HaplotypeCallerGvcf_GATK4.output_merged_gvcf": "/genomes/project/wdl/simple-pipeline-test/cromwell-executions/HaplotypeCallerGvcf_GATK4/621bcdec-19a9-411d-8670-3899fa9e8509/call-MergeGVCFs/execution/NA12878_24RG_small.b37.g.vcf.gz"
}
[2019-05-31 06:41:55,51] [info] WorkflowManagerActor WorkflowActor-621bcdec-19a9-411d-8670-3899fa9e8509 is in a terminal state: WorkflowSucceededState
[2019-05-31 06:42:20,69] [info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "HaplotypeCallerGvcf_GATK4.output_merged_gvcf": "/genomes/project/wdl/simple-pipeline-test/cromwell-executions/HaplotypeCallerGvcf_GATK4/621bcdec-19a9-411d-8670-3899fa9e8509/call-MergeGVCFs/execution/NA12878_24RG_small.b37.g.vcf.gz",
    "HaplotypeCallerGvcf_GATK4.output_merged_gvcf_index": "/genomes/project/wdl/simple-pipeline-test/cromwell-executions/HaplotypeCallerGvcf_GATK4/621bcdec-19a9-411d-8670-3899fa9e8509/call-MergeGVCFs/execution/NA12878_24RG_small.b37.g.vcf.gz.tbi"
  },
  "id": "621bcdec-19a9-411d-8670-3899fa9e8509"
}
```