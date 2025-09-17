RACKJ (Read Analysis & Comparison Kit in Java) is a set of Java programs that analyze and compare RNA-seq data made by NGS (Next-Generation Sequencing) technologies. In addition to RPKM (Reads Per Kbp per Million reads) values, RACKJ computes read counts for exons and splicing events. In so doing, it is feasible to compare two samples and identify genes with most significant difference in exon(splicing)-level.

This repository is currently for [rackj](https://sourceforge.net/projects/rackj/) at SourceForge. We would setup manual pages here first, and then programs.

## Installation

It would be most convenient to use the [docker image](https://hub.docker.com/r/wdlin/rackj). Programs and packages required in our example walkthroughs should be contained in the image file.

<details>
  <summary>Manual install instruction</summary>
For manual installation, the following steps should work for most of our programs in Ubuntu 20.04.

```
# Assume username is ubuntu and at home directory
# Create a directory as you like and enter it
mkdir Tools
cd Tools/

# Download rackJ (check https://sourceforge.net/projects/rackj/files/ for any updated versions) and extract
wget https://downloads.sourceforge.net/project/rackj/0.99a/rackJ.tar.gz
tar -zxvf rackJ.tar.gz

# Down related package, picard
wget https://downloads.sourceforge.net/project/picard/sam-jdk/1.89/sam-1.89.jar

# make sure R installed and available from command line
# example steps from R offical website for ubuntu
sudo apt update -qq
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt install --no-install-recommends r-base

# Install bioperl, necessary perl libraries, make, and java11
sudo apt install bioperl make openjdk-11-jdk-headless
sudo cpan Statistics::Distributions
sudo cpan Statistics::R

# The MappingBowtie.pl will call bowtie2. Install if you need it.
sudo apt install unzip
sudo apt install python
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip
unzip bowtie2-2.4.4-linux-x86_64.zip

# The MappbinBlat.pl will call blat. Install if you need it.
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod 755 blat

# Set PATH environment into ~/.profile, adjust any path if necessary
echo "PATH=\"/home/ubuntu/Tools/rackJ/scripts:\$PATH\"" >> ~/.profile
echo "PATH=\"/home/ubuntu/Tools/bowtie2-2.4.4-linux-x86_64:\$PATH\"" >> ~/.profile
echo "PATH=\"/home/ubuntu/Tools:\$PATH\"" >> ~/.profile
echo "INC=\"/home/ubuntu/Tools/rackJ/scripts/PerlLib\"" >> ~/.profile
echo "CLASSPATH=\"/home/ubuntu/Tools/rackJ/rackj.jar:/home/ubuntu/Tools/sam-1.89.jar\"" >> ~/.profile
echo "export PATH" >> ~/.profile
echo "export INC" >> ~/.profile
echo "export CLASSPATH" >> ~/.profile

# Source if needed
source ~/.profile
```
</details>
