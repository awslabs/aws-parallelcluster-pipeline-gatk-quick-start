# GATK Best Practice on AWS

+ [Chinese](./README.md)
+ [English](./README-en.md)

* [一、方案简介](#一方案简介)
    * [1、现状及方案概述](#1现状及方案概述)
    * [2、解决的需求和问题](#2解决的需求和问题)
      * [1)、生产与测试系统未隔离，各产品线共用一套集群](#1生产与测试系统未隔离各产品线共用一套集群)
      * [2)、本地集群无容灾，一旦受损（宕机、环境变动）、停止（如物业计划维护、断网断电等）将影响所有产品业务](#2本地集群无容灾一旦受损宕机环境变动停止如物业计划维护断网断电等将影响所有产品业务)
      * [3)、计算存储资源量估算困难，部署周期长，业务线变动较大](#3计算存储资源量估算困难部署周期长业务线变动较大)
      * [4)、开发迭代混乱，集群环境大家共同维护，难以做到版本控制和环境复现!](#4开发迭代混乱集群环境大家共同维护难以做到版本控制和环境复现)
      * [5)、行业其他问题：](#5行业其他问题)
          * [①、项目管理靠人工跟踪，无信息化系统或信息化程度低](#项目管理靠人工跟踪无信息化系统或信息化程度低)
          * [②、信息环节无法做到成本核算](#信息环节无法做到成本核算)
          * [③、数据量大且数据复杂，难以做到生命周期管理](#数据量大且数据复杂难以做到生命周期管理)
* [二、详细方案说明](#二详细方案说明)
    * [1、集群解决方案](#1集群解决方案)
    * [2、端到端的基因大数据分析、归档、交付方案](#2端到端的基因大数据分析归档交付方案)
* [三、迁移部署流程](#三迁移部署流程)
* [四、测试文档](#四测试文档)
    * [1、10分钟集群部署](#110分钟集群部署)
      * [1)、 安装pip及awscli并配置必要信息](#1-安装pip及awscli并配置必要信息)
      * [2)、安装pcluster](#2安装pcluster)
      * [3)、pcluster安装及配置](#3pcluster安装及配置)
          * [①、相关准备](#相关准备)
          * [②、配置pcluster config(可参考官方博客)](#配置pcluster-config可参考官方博客)
      * [4)、启动集群](#4启动集群)
      * [5)、登陆集群master节点](#5登陆集群master节点)
      * [6)、投递任务](#6投递任务)
    * [2、AMI](#2ami)
    * [3、DEMO](#3demo)
* [五、参考资料：](#四参考资料)
* [六、FAQ](#五faq)


## 一、方案简介

***相关基本介绍请参考[Parallelcluster官方博客](https://amazonaws-china.com/cn/blogs/china/aws-parallelcluster/).***

### 1、现状及方案概述
当前基因行业的分析平台大致分为三类，分别是单机、HPC集群、K8S集群，占比最的的仍然是HPC集群这种形式。
对于本地HPC环境来说可能会遇到以下问题：

+ 系统部署周期长(从采购到部署长达数个月)，难以跟上业务发展的速度
+ 集群运维难度较大，对于业务时效性要求较高的场景，难以做到高可用
+ 数据的生命周期管理有很大的挑战，且成本难以控制
+ 开发人员被各种集群问题所困扰，运维人员的精力被日常清理工作所消耗，难以集中做集群的优化工作

对于传统云上HPC环境来说可能会遇到以下问题：
+ 资源难以精确控制，计算节点存在资源浪费
+ 运维部署周期长，一般需要1～2周，运维难度较大，环境不便复制迁移


除此之外还有很多问题会困扰大家,为了帮助克服、解决这些问题，我们提供了一套简单易用的方案，可一键创建完整的HPC集群，除此之外通过参数配置的调整，还可以以用户习惯的架构来精细化调整HPC环境，通过模版的复制，也可以最小化的成本实现本地的测试、迁移以及云上的集群环境复制，以此来帮助生命科学领域的用户去构建安全、可靠、高效、低成本的HPC集群，将开发、运维人员的精力从琐事中解放，专注在更有创造力的事情上，借助云的优势，可以打造更灵活、高可用的系统架构。

本方案的核心是parallelcluster,它自带的jobwatcher和配套服务，可以每分钟监控SGE、Slurm或Torque作业情况，以决定节点的何时进行弹性伸缩，这样可以直接带来30%左右成本的节约，具体工作原理参见 ***[Parallelcluster用户文档——作业处理流程](https://docs.aws.amazon.com/zh_cn/parallelcluster/latest/ug/processes.html).***


***

### 2、解决的需求和问题
#### 1)、生产与测试系统未隔离，各产品线共用一套集群
#### 2)、本地集群无容灾，一旦受损（宕机、环境变动）、停止（如物业计划维护、断网断电等）将影响所有产品业务
![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/1.png)

#### 3)、计算存储资源量估算困难，部署周期长，业务线变动较大
![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/2.png)

#### 4)、开发迭代混乱，集群环境大家共同维护，难以做到版本控制和环境复现!
![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/3.png)

#### 5)、行业其他问题：
##### ①、项目管理靠人工跟踪，无信息化系统或信息化程度低
##### ②、信息环节无法做到成本核算
##### ③、数据量大且数据复杂，难以做到生命周期管理

***

## 二、详细方案说明
### 1、集群解决方案
![集群解决方案](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/4.png)
### 2、端到端的基因大数据分析、归档、交付方案
![端到端的基因大数据分析、归档、交付方案](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/5.png)
## 三、迁移部署流程
[建议迁移部署方法](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/doc/01.migrate/README.md)
## 四、测试文档
### 1、10分钟集群部署
下述文档示例会启动一个完整的HPC集群，包括主节点、计算节点、共享存储以及预装SGE作业调度系统，AMI为预装GATK相关软件的镜像，包括bwa，Samtools，gatk4等，镜像snapshot为GATK公开数据集，包括数据库及测试文件，启动后挂载到/genomics目录下。

***[注：测试以宁夏区为例，本测试需要账号有多个资源的创建权限(初期测试可用admin access)](https://docs.aws.amazon.com/zh_cn/parallelcluster/latest/ug/iam.html)***

#### 1)、 安装pip及awscli并配置必要信息

+ [awscli安装参考](https://docs.aws.amazon.com/zh_cn/cli/latest/userguide/cli-chap-install.html)
+ aws_access_key_id及aws_secret_access_key


**请登录console，并点击*我的安全凭证***

**创建访问密钥并记录aws_access_key_id及aws_secret_access_key**

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/6.png)

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/7.png)

<br>
<br>

```shell
#pip安装
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py

#awscli安装
sudo pip install awscli

#aws配置
#根据具体情况设置AK,SK，所在区域及输出格式
aws configure 
```
#### 2)、安装pcluster
```shell
sudo pip install aws-parallelcluster
```

#### 3)、pcluster安装及配置
##### ①、相关准备

+ VPC及子网

**请搜索VPC服务，并选择你要启动集群的VPC**

**记录vpc_id，master_subnet_id**

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/8.png)

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/9.png)

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/10.png)

<br>
<br>

+ EC2访问密钥

**搜索EC2服务，并进入EC2服务页面选择密钥对，创建新密钥并下载密钥文件**

**记录key_name为创建密钥对输入的名字**

<br>

![](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/blob/master/images/11.png)

##### ②、配置pcluster config([可参考官方博客](https://docs.aws.amazon.com/zh_cn/parallelcluster/latest/ug/configuration.html))

```shell
#创建配置模版,会提示无配置，请忽略错误信息
pcluster create new

#编辑配置文档
vim ~/.parallelcluster/config

#复制下述配置信息，粘贴到配置文档~/.parallelcluster/config
[aws]
aws_region_name = cn-northwest-1

[global]
update_check = true
sanity_check = true
cluster_template = GATK-pipeline

[aliases]
ssh = ssh {CFN_USER}@{MASTER_IP} {ARGS}

[cluster GATK-pipeline]
base_os = alinux
custom_ami = ami-056169db492793d02 #根据需要修改
vpc_settings = public
scheduler = slurm
key_name = ZHY_key  #需要修改
compute_instance_type = m5.xlarge
master_instance_type = m5.xlarge
compute_root_volume_size = 50
master_root_volume_size = 50
ebs_settings = genomes
scaling_settings = GATK-ASG
initial_queue_size = 1
max_queue_size = 4
maintain_initial_size = false
extra_json = { "cluster" : { "cfn_scheduler_slots" : "cores" } }

[vpc public]
vpc_id = vpc-a817aac5  #需要修改
master_subnet_id = subnet-26fcc86cd  #需要修改

[ebs genomes]
shared_dir = genomes
ebs_snapshot_id = snap-040c71fd2bb5d4236 #根据需要修改
volume_type = gp2
volume_size =  1024

[scaling GATK-ASG]
scaledown_idletime = 5
```
    
#### 4)、启动集群
```shell
pcluster create GATK-pipeline
```
    
#### 5)、登陆集群master节点
```shell
#根据集群启动后的反馈信息输入
#ssh -i <private key_name> <username>@<public ip>
ssh -i <private key_name> ec2-user@master-public-ip #alinux
ssh -i <private key_name> ubuntu@master-public-ip #ubuntu
ssh -i <private key_name> centos@master-public-ip #centos
```
    
#### 6)、投递任务
默认预装SGE作业调度系统，示例sge调度系统投递命令参考如下：
```
echo "sleep 180" | qsub
echo "sh run.sh" | qsub -l vf=2G,s_core=1 -q all.q
for((i=1;i<=10;i++));do echo "sh /genomes/temp/run.sh $i" | qsub -cwd -S /bin/bash -l vf=2G,s_core=1 -q all.q;done
```    

示例slurm调度系统投递命令参考如下：
```shell
sbatch -n 4 run.sh  #4核，可根据需要修改
squeue #查看队列情况
sinfo #查看节点情况
scancel jobid #取消任务
```

示例pbs调度系统投递命令参考如下：
```shell
echo "sleep 180" | qsub
echo "sh run.sh" | qsub -l nodes=1,walltime=2:00:00,mem=2gb -q batch
for((i=1;i<=10;i++));do echo "sh /genomes/temp/run.sh $i" | qsub -l nodes=1,walltime=2:00:00,mem=2gb -q batch;done
```

### 2、AMI
***制作自定义ami请参考官方文档***

+ [alinux](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/AMI/alinux.sh)
+ [ubuntu](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/AMI/ubuntu.sh)

### 3、DEMO

+ [基于nextflow工具调度的demo](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/example/nextflow/README.md)
+ [基于cromwell工具调度的demo](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/example/cromwell/README.md)
+ [基于shell的demo](https://github.com/awslabs/aws-parallelcluster-pipeline-gatk-quick-start/tree/master/example/shell/README.md)


## 五、参考资料：

+ [Parallelcluster官方博客.](https://amazonaws-china.com/cn/blogs/china/aws-parallelcluster/)
+ [parallelcluster文档.](https://docs.aws.amazon.com/zh_cn/parallelcluster/latest/ug/what-is-aws-parallelcluster.html)
+ [aws-parallelcluster GitHub 存储库.](https://github.com/aws/aws-parallelcluster)
+ 当前版本
  
|系统	|版本号	|pcluster版本	|AMI ID	|更新描述	|地域	|是否公开	|可用性	|备注	|
|---	|---	|---	|---	|---	|---	|---	|---	|---	|
|alinux	|0.2	|2.4.1	|ami-ami-005d4f6437dca2d6d	|基础软件环境AMI + Golang环境 + goofys	|BJS	|是	|是	|	|
|alinux	|0.2	|2.4.1	|ami-056169db492793d02	|基础软件环境AMI + Golang环境 + goofys	|ZHY	|是	|是	|	|

+ 镜像版本迭代

|系统	|版本号	|pcluster版本	|AMI ID	|更新描述	|地域	|是否公开	|可用性	|备注	|
|---	|---	|---	|---	|---	|---	|---	|---	|---	|
|alinux-base	|	|2.3.1	|ami-0e58e06d5b958ccb6	|基础镜像	|BJS	|是	|是	|	|
|ubuntu-base	|16.04	|2.3.1	|ami-0a9c1879e6583621e	|基础镜像	|BJS	|是	|是	|	|
|alinux	|0.1	|2.3.1	|ami-0997595bce93c6e7b	|基础软件环境AMI	|BJS	|是	|是	|	|
|alinux	|0.2	|2.3.1	|ami-0cad4e9d804bd9c15	|基础软件环境AMI + Golang环境 + goofys；修复pip问题并安装awscli；修复goofys无法挂载问题，安装fuse依赖	|BJS	|是	|是	|	|
|alinux	|0.2	|2.4.0	|ami-0b876120ec98b9a7c	|基础软件环境AMI	|BJS	|是	|是	|	|
|ubuntu	|0.1	|2.3.1	|ami-097d3bf901991372e	|基础软件环境AMI	|BJS	|是	|是	|	|
|ubuntu	|0.2	|2.3.1	|ami-041e4a3bce09385b9	|修改ubuntu默认shell(dash)为bash	|BJS	|是	|是	|不再更新	|
|ubuntu	|0.2-a	|2.3.1	|ami-026882b56146cdc1b	|基础软件环境AMI + Golang环境 + goofys	|BJS	|是	|是	|	|
|alinux	|0.1	|2.3.1	|ami-007f6ed61542ae017	|基础软件环境AMI	|ZHY	|是	|是	|	|
|alinux	|0.2	|2.4.0	|ami-005db8a58ebd4e9a4	|基础软件环境AMI	|ZHY	|是	|是	|	|
|ubuntu	|0.1	|2.3.1	|ami-0a1d99c2c70e3f86c	|基础软件环境AMI	|ZHY	|是	|否	|	|
|ubuntu	|0.2	|2.3.1	|ami-071aa7a2927cc02a8	|修改ubuntu默认shell(dash)为bash	|ZHY	|是	|否	|不再更新	|
|ubuntu	|0.2-a	|2.3.1	|ami-015f3a018cc98b6cc	|基础软件环境AMI + Golang环境 + goofys	|ZHY	|是	|否	|	|

+ 参考文件EBS快照迭代：

|名称	|版本号	|snap ID	|大小	|更新描述	|地域	|
|---	|---	|---	|---	|---	|---	|
|gatk-reference v0.1	|0.1	|snap-09c16ac9809cf4359	|100G	|基础环境快照，包含hg19参考基因组	|BJS	|
|gatk-reference v0.2	|0.2	|snap-06f5e874571e44510	|100G	|增加hg38及GATK数据集	|BJS	|
|gatk-reference-v0.3	|0.3	|snap-08a4b975a2f40736f	|1T	|增加测试文件及GATK-TEST-DATA	|BJS	|
|gatk-reference-v0.3	|0.3	|snap-040c71fd2bb5d4236	|1T	|增加测试文件及GATK-TEST-DATA	|ZHY	|
|	|	|	|	|	|	|

## 六、FAQ

