# Comparing isoform expression levels and ratios (against gene expression level) using StringTie

This page will go through the following items using short read pair-ended datasets of three control samples and three treatment samples:
1. isoform read count estimation and gene read count computation
2. isoform expression level comparison
3. isofrom ratio (isoform against gene) comparison

We will use ExampleData.zip and the [docker image](https://hub.docker.com/r/wdlin/rackj) for all the programs. In this walkthrough, we will use Singularity to run the docker image. Usage example:

```
wdlin@comp10:/RAID2/R418/20250930_AS$ unzip ExampleData.zip
Archive:  ExampleData.zip
   creating: ExampleData/
 extracting: ExampleData/control_rep1.merged.bam
 extracting: ExampleData/control_rep2.merged.bam
 extracting: ExampleData/control_rep4.merged.bam
 extracting: ExampleData/README.txt
   creating: ExampleData/src/
 extracting: ExampleData/src/control_rep1_R1.fq.gz
 extracting: ExampleData/src/control_rep1_R2.fq.gz
 extracting: ExampleData/src/control_rep2_R1.fq.gz
 extracting: ExampleData/src/control_rep2_R2.fq.gz
 extracting: ExampleData/src/control_rep4_R1.fq.gz
 extracting: ExampleData/src/control_rep4_R2.fq.gz
 extracting: ExampleData/src/treatment_rep5_R1.fq.gz
 extracting: ExampleData/src/treatment_rep5_R2.fq.gz
 extracting: ExampleData/src/treatment_rep7_R1.fq.gz
 extracting: ExampleData/src/treatment_rep7_R2.fq.gz
 extracting: ExampleData/src/treatment_rep9_R1.fq.gz
 extracting: ExampleData/src/treatment_rep9_R2.fq.gz
  inflating: ExampleData/tair10.strand.cgff
  inflating: ExampleData/tair10.strand.model
  inflating: ExampleData/TAIR10_chr_all.fas
  inflating: ExampleData/TAIR10_GFF3_genes_transposons.gff
 extracting: ExampleData/treatment_rep5.merged.bam
 extracting: ExampleData/treatment_rep7.merged.bam
 extracting: ExampleData/treatment_rep9.merged.bam

# this is to make sure the folder writable
wdlin@comp10:/RAID2/R418/20250930_AS$ chmod -R 755 ExampleData

wdlin@comp10:/RAID2/R418/20250930_AS$ cd ExampleData/

# this binds the ExampleData folder to /mnt in the container
wdlin@comp10:/RAID2/R418/20250930_AS/ExampleData$ singularity run --bind "$PWD:/mnt" docker://wdlin/rackj
INFO:    Using cached SIF image

Singularity> cd /mnt/
Singularity> ls
README.txt                         TAIR10_chr_all.fas       control_rep2.merged.bam  src                 tair10.strand.model        treatment_rep7.merged.bam
TAIR10_GFF3_genes_transposons.gff  control_rep1.merged.bam  control_rep4.merged.bam  tair10.strand.cgff  treatment_rep5.merged.bam  treatment_rep9.merged.bam
```

Note that the `ExampleData` folder was bounded as `/mnt` in the container. All necessary programs should be available so no need to do any installation. Also note the raw reads in this dataset contains only very small part of adapters so the adapter removal was not applied.

## Mapping using TopHat2

Since we are going to map reads, existing BAM files are not needed.
```
rm *.bam
```

Build Bowtie2 genome index.
```
bowtie2-build TAIR10_chr_all.fas tair10.genome
```

Build TopHat2 transcriptome index. This is strongly suggested although TopHat2 would automatically build the index. Prebuilding the index would save time and avoid possible race condition. 
```
tophat2 -G TAIR10_GFF3_genes_transposons.gff --transcriptome-index=tair10.transcriptome/known tair10.genome
```

Since we have pair-ended raw read files `src/*.fq.gz`,
```
Singularity> ls -l src/*.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21454332 Oct  4  2021 src/control_rep1_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 22753585 Oct  4  2021 src/control_rep1_R2.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21983578 Oct  4  2021 src/control_rep2_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 23494113 Oct  4  2021 src/control_rep2_R2.fq.gz
-rwxr-xr-x+ 1 wdlin R418 20786103 Oct  4  2021 src/control_rep4_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21921178 Oct  4  2021 src/control_rep4_R2.fq.gz
-rwxr-xr-x+ 1 wdlin R418 20703897 Oct  4  2021 src/treatment_rep5_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21943326 Oct  4  2021 src/treatment_rep5_R2.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21909066 Oct  4  2021 src/treatment_rep7_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 23300279 Oct  4  2021 src/treatment_rep7_R2.fq.gz
-rwxr-xr-x+ 1 wdlin R418 21281263 Oct  4  2021 src/treatment_rep9_R1.fq.gz
-rwxr-xr-x+ 1 wdlin R418 22591761 Oct  4  2021 src/treatment_rep9_R2.fq.gz
```

the following perl one-liner can help us to form 6 commands of running TopHat2. Note that option `--transcriptome-index=tair10.transcriptome/known` was provided for the transcriptome index and `tair10.genome` was for the genome index. Option `-p 4` was for using 4 threads for mapping, adjust it if needed.
```
Singularity> ls src/*.fq.gz | sort | perl -ne 'chomp; /.+\/(.+)_R\d\./; push @{$hash{$1}},$_; if(eof){ for $k (sort keys %hash){ $cmd="tophat2 -o $k"."_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome @{$hash{$k}}"; print "\nCMD: $cmd\n"; } }'

CMD: tophat2 -o control_rep1_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/control_rep1_R1.fq.gz src/control_rep1_R2.fq.gz

CMD: tophat2 -o control_rep2_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/control_rep2_R1.fq.gz src/control_rep2_R2.fq.gz

CMD: tophat2 -o control_rep4_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/control_rep4_R1.fq.gz src/control_rep4_R2.fq.gz

CMD: tophat2 -o treatment_rep5_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/treatment_rep5_R1.fq.gz src/treatment_rep5_R2.fq.gz

CMD: tophat2 -o treatment_rep7_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/treatment_rep7_R1.fq.gz src/treatment_rep7_R2.fq.gz

CMD: tophat2 -o treatment_rep9_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome src/treatment_rep9_R1.fq.gz src/treatment_rep9_R2.fq.gz
```

The commands seem OK. So we add `system $cmd` for actual executing them.
```
ls src/*.fq.gz | sort | perl -ne 'chomp; /.+\/(.+)_R\d\./; push @{$hash{$1}},$_; if(eof){ for $k (sort keys %hash){ $cmd="tophat2 -o $k"."_tophat2 -p 4 --transcriptome-index=tair10.transcriptome/known tair10.genome @{$hash{$k}}"; print "\nCMD: $cmd\n"; system $cmd } }'
```
