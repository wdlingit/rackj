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

## 1. Mapping using TopHat2

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

Output BAM files would be `*_tophat2/accepted_hits.bam`.
```
Singularity> ls -l *_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 40533905 Sep 18 14:43 control_rep1_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 41634756 Sep 18 14:45 control_rep2_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 39169407 Sep 18 14:47 control_rep4_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 39090724 Sep 18 14:49 treatment_rep5_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 41104212 Sep 18 14:50 treatment_rep7_tophat2/accepted_hits.bam
-rwxrwxrwx+ 1 wdlin R418 40041940 Sep 18 14:52 treatment_rep9_tophat2/accepted_hits.bam
```

## 2. Guided assembly by StringTie

Firstly we translate the `TAIR10_GFF3_genes_transposons.gff` into GTF file `tair10.gtf` using `gffread`. Note that the `2>/dev/null` part is to direct warning messages of `transposon_fragment` to `/dev/null`. We omit them because we know they are harmless.
```
gffread -T -o tair10.gtf TAIR10_GFF3_genes_transposons.gff 2>/dev/null
```

StringTie-guided assembly has two steps. The first step is to generate a GTF file for each of the BAM file.
```
ls *_tophat2/accepted_hits.bam | perl -ne 'chomp; /(.+?)_tophat2/; $cmd="stringtie $_ -o $1.gtf -p 4 -G tair10.gtf"; print "\nCMD: $cmd\n"; system $cmd'
```

The above perl one-liner should give us the following 6 commands. Each of them reads the read alignment BAM file and does guided assembly with regards to the existing genome annotation (`tair10.gtf`).
```
CMD: stringtie control_rep1_tophat2/accepted_hits.bam -o control_rep1.gtf -p 4 -G tair10.gtf

CMD: stringtie control_rep2_tophat2/accepted_hits.bam -o control_rep2.gtf -p 4 -G tair10.gtf

CMD: stringtie control_rep4_tophat2/accepted_hits.bam -o control_rep4.gtf -p 4 -G tair10.gtf

CMD: stringtie treatment_rep5_tophat2/accepted_hits.bam -o treatment_rep5.gtf -p 4 -G tair10.gtf

CMD: stringtie treatment_rep7_tophat2/accepted_hits.bam -o treatment_rep7.gtf -p 4 -G tair10.gtf

CMD: stringtie treatment_rep9_tophat2/accepted_hits.bam -o treatment_rep9.gtf -p 4 -G tair10.gtf
```

The second step is to merge all GTF file for the 6 samples. The merged GTF file is named `merged.gtf`.
```
stringtie --merge -G tair10.gtf -o merged.gtf control_rep1.gtf control_rep2.gtf control_rep4.gtf treatment_rep5.gtf treatment_rep7.gtf treatment_rep9.gtf
```

## 2A. Not applying guided assembly by StringTie

The guided assembly step by StringTie would produce many novel transcripts and novel gene loci. In case that you don't want those novel annotation or you have a genome annotation file made by high confident source (ex: ISOseq), you can just replace `merged.gtf` with your GTF file. For example, we can replace `merged.gtf` with `tair10.gtf` if we don't want to do the guided assembly.

## 3. Generate count matrixes using StringTie

The following step should generate count files for each of the samples.
```
ls *_tophat2/accepted_hits.bam | perl -ne 'chomp; /(.+?)_tophat2/; $cmd="stringtie $_ -eB -o $1/$1.gtf -p 4 -G merged.gtf"; print "\nCMD: $cmd\n"; system $cmd'
```

The above perl one-liner should give us the following 6 commands. Each of them means "generate transcript read counts and gene read counts using the genome annotation file `merged.gtf`".
```
CMD: stringtie control_rep1_tophat2/accepted_hits.bam -eB -o control_rep1/control_rep1.gtf -p 4 -G merged.gtf

CMD: stringtie control_rep2_tophat2/accepted_hits.bam -eB -o control_rep2/control_rep2.gtf -p 4 -G merged.gtf

CMD: stringtie control_rep4_tophat2/accepted_hits.bam -eB -o control_rep4/control_rep4.gtf -p 4 -G merged.gtf

CMD: stringtie treatment_rep5_tophat2/accepted_hits.bam -eB -o treatment_rep5/treatment_rep5.gtf -p 4 -G merged.gtf

CMD: stringtie treatment_rep7_tophat2/accepted_hits.bam -eB -o treatment_rep7/treatment_rep7.gtf -p 4 -G merged.gtf

CMD: stringtie treatment_rep9_tophat2/accepted_hits.bam -eB -o treatment_rep9/treatment_rep9.gtf -p 4 -G merged.gtf
```

StringTie offeres a convenient command `prepDE.py` to generate count matrixes for transcripts and genes, given that all above commands were correctly done and there are 6 files by `ls */*.gtf`.
```
Singularity> ls */*.gtf
control_rep1/control_rep1.gtf  control_rep4/control_rep4.gtf      treatment_rep7/treatment_rep7.gtf
control_rep2/control_rep2.gtf  treatment_rep5/treatment_rep5.gtf  treatment_rep9/treatment_rep9.gtf

Singularity> prepDE.py

Singularity> head transcript_count_matrix.csv
transcript_id,control_rep1,control_rep2,control_rep4,treatment_rep5,treatment_rep7,treatment_rep9
AT4G04480.1,0,0,0,0,0,0
AT1G07730.2,0,0,0,0,0,2
AT1G38430.1,0,0,0,0,0,0
AT1G03340.1,7,4,4,4,4,11
AT2G25040.1,0,0,0,0,0,0
AT1G04440.1,52,68,88,106,85,113
AT5G13090.1,34,24,23,19,29,19
AT2G30190.1,0,0,0,0,0,0
AT1G31390.1,0,0,0,0,0,0

Singularity> head gene_count_matrix.csv
gene_id,control_rep1,control_rep2,control_rep4,treatment_rep5,treatment_rep7,treatment_rep9
AT5G59170,0,0,0,0,0,0
AT1G38440,0,0,0,0,0,0
AT2G46660,0,0,0,0,0,0
AT5G25130,7,11,7,3,4,4
AT5G54062,0,0,0,0,0,0
AT5G61100,0,0,0,0,0,0
AT5G06000,4,0,4,0,0,0
AT4G25040,0,0,0,0,0,0
MSTRG.10381,44,24,25,27,22,45
```

Note that `prepDE.py` is for Python2. There is another script named `prepDE.py3` for Python3. Also note that `MSTRG.10381` is a novel gene locus predicted by StringTie in case that the guided assembly was done as described in above step 2.

## 4. Isoform expression level comparison

Since we have three control samples and three treatment samples where their transcript read counts are stored in file `transcript_count_matrix.csv`, the following R code should work for comparing the 3 treatment samples against the 3 control samples using `DESeq2`.

```
Singularity> R

> library(DESeq2)

> countData <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))
> head(countData)
            control_rep1 control_rep2 control_rep4 treatment_rep5
AT4G04480.1            0            0            0              0
AT1G07730.2            0            0            0              0
AT1G38430.1            0            0            0              0
AT1G03340.1            7            4            4              4
AT2G25040.1            0            0            0              0
AT1G04440.1           52           68           88            106
            treatment_rep7 treatment_rep9
AT4G04480.1              0              0
AT1G07730.2              0              2
AT1G38430.1              0              0
AT1G03340.1              4             11
AT2G25040.1              0              0
AT1G04440.1             85            113

# We do this because we know the first 3 columns are for the 3 control samples
# where "A" will be the reference
> condition= c("A","A","A","B","B","B")
> df = data.frame(condition,row.names=colnames(countData))

> dds <- DESeqDataSetFromMatrix(countData,colData=df,design=~condition)
> dds <- DESeq(dds)

> resultsNames(dds)
[1] "Intercept"        "condition_B_vs_A"

> write.csv(as.data.frame(results(dds,name="condition_B_vs_A")),file="desqOut.csv")

> q()
Save workspace image? [y/n/c]: n

Singularity> head desqOut.csv
"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"AT4G04480.1",0,NA,NA,NA,NA,NA
"AT1G07730.2",0.322550149365598,1.8429268299077,4.03846570623945,0.45634331549736,0.648143120481551,NA
"AT1G38430.1",0,NA,NA,NA,NA,NA
"AT1G03340.1",5.52317051161533,0.331267683593407,1.30268193524514,0.254296674138703,0.799266369024993,0.999939960474882
"AT2G25040.1",0,NA,NA,NA,NA,NA
"AT1G04440.1",83.4751224366111,0.53923757786864,0.280024189372634,1.92568213152138,0.0541440763572674,0.555418282168473
"AT5G13090.1",24.0349600628139,-0.285104934557586,0.517964774480711,-0.55043305762138,0.582022380056347,0.999939960474882
"AT2G30190.1",0,NA,NA,NA,NA,NA
"AT1G31390.1",0,NA,NA,NA,NA,NA
```
