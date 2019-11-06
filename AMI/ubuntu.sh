sudo -s
ln -sf /bin/bash /bin/sh
#library
apt-get update && apt install zlib1g-dev gcc libncurses5-dev libncursesw5-dev libboost-all-dev libbz2-dev liblzma-dev python-pip unzip python-software-properties software-properties-common openjdk-8-jdk-headless
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && apt-get install git-lfs

mkdir -p /Pipeline/01.software
mkdir -p /Pipeline/02.resource

#BWA
cd /Pipeline/01.software && git clone https://github.com/lh3/bwa
cd bwa && make
ln -s /Pipeline/01.software/bwa/bwa /usr/bin/

#Samtools
cd /Pipeline/01.software && git clone https://github.com/samtools/htslib.git
cd /Pipeline/01.software && git clone https://github.com/samtools/samtools.git
cd /Pipeline/01.software/samtools && make
ln -s /Pipeline/01.software/samtools/samtools /usr/bin/

#GATK
cd /Pipeline/01.software && git clone https://github.com/broadinstitute/gatk.git
cd /Pipeline/01.software/gatk/ && ./gradlew installAll
ln -s /Pipeline/01.software/gatk/gatk /usr/bin/

#Picard
cd /Pipeline/01.software && git clone https://github.com/broadinstitute/picard.git
cd /Pipeline/01.software/picard/ && ./gradlew shadowJar


#Golang
mkdir /Pipeline/01.software/go
cd /Pipeline/01.software/go && echo 'export GOPATH='`pwd` >> ~/.bashrc && source ~/.bashrc 
wget https://dl.google.com/go/go1.12.5.linux-amd64.tar.gz && tar -C `pwd` -xzf go1.12.5.linux-amd64.tar.gz
echo 'export PATH=$PATH:/Pipeline/01.software/go/go/bin/' >> ~/.bashrc && source ~/.bashrc

#goofys
mkdir -p $GOPATH/src/golang.org/x/ && cd $GOPATH/src/golang.org/x/
git clone https://github.com/golang/net.git $GOPATH/src/golang.org/x/net
git clone https://github.com/golang/text.git $GOPATH/src/golang.org/x/text
git clone https://github.com/golang/sys.git $GOPATH/src/golang.org/x/sys
git clone https://github.com/golang/crypto.git $GOPATH/src/golang.org/x/crypto
git clone https://github.com/grpc/grpc-go.git $GOPATH/src/google.golang.org/grpc
go get -u github.com/golang/protobuf/{proto,protoc-gen-go}
git clone https://github.com/google/go-genproto.git $GOPATH/src/google.golang.org/genproto
cd $GOPATH/src/
go install google.golang.org/grpc
git clone https://github.com/golang/sync.git $GOPATH/src/golang.org/x/sync
git clone https://github.com/golang/oauth2.git $GOPATH/src/golang.org/x/oauth2
git clone https://github.com/googleapis/google-api-go-client.git $GOPATH/src/google.golang.org/api
go install google.golang.org/api
go install net
go get github.com/kahing/goofys && go install github.com/kahing/goofys
echo 'export PATH=$PATH:/Pipeline/01.software/go/bin/' >> ~/.bashrc && source ~/.bashrc