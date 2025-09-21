# Comparing alternative-splicing events using RackJ

This page will go through the following items using short read pair-ended datasets of three control samples and three treatment samples:
1. alternative-splicing event comparison based on read counts
2. alternative-splicing evnet comparison based on ratios with respects to biological replicates

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

For Docker and Docker Desktop user:
1. extract ExampleData.zip and remember the path to the ExampleData folder
2. (only for Docker Desktop) open _terminal_ (buttom-right corner of Docker Desktop)
3. run command `docker run -it --rm --mount type=bind,src=<path_to_ExampleData>,dst=/mnt wdlin/rackj`
4. just remember that the ExampleData folder is bound to `/mnt`

The maximum memory usage is about 2GB for this workthrough so it seems not needed to adjust any resource limits.

## 1. Mapping using BLAT

This is a optional step, you may adopt `*.merged.bam` in ExampleData.zip directly. 

In case that we are going to map reads, existing BAM files are not needed.
```
rm *.bam
```

Build indexes for `BioPerl` and `samtools`. This is strongly suggested for avoiding race condition in case you are going to submit mapping commands to a job scheduler like slurm.
```
preindex.pl TAIR10_chr_all.fas
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

the following perl one-liner can help us to form 12 commands of running BLAT via rackj scripts `Mapping.pl` and `MappingBlat.pl`. Note that option option `-split 4` was for using 4 processes for mapping, adjust it if needed.
```
Singularity> ls src/*.fq.gz | perl -ne 'chomp; /.+\/(.+?)\./; $cmd="gzip -dc $_ > $1.fq; Mapping.pl -split 4 x $1.fq $1.blat.bam MappingBlat.pl -target TAIR10_chr_all.fas -t=dna -q=rna; rm $1.fq"; print "\nCMD: $cmd\n";'

CMD: gzip -dc src/control_rep1_R1.fq.gz > control_rep1_R1.fq; Mapping.pl -split 4 x control_rep1_R1
.fq control_rep1_R1.blat.bam MappingBlat.pl -target TAIR10_chr_all.fas -t=dna -q=rna; rm control_re
p1_R1.fq

CMD: gzip -dc src/control_rep1_R2.fq.gz > control_rep1_R2.fq; Mapping.pl -split 4 x control_rep1_R2
.fq control_rep1_R2.blat.bam MappingBlat.pl -target TAIR10_chr_all.fas -t=dna -q=rna; rm control_re
p1_R2.fq

CMD: gzip -dc src/control_rep2_R1.fq.gz > control_rep2_R1.fq; Mapping.pl -split 4 x control_rep2_R1
.fq control_rep2_R1.blat.bam MappingBlat.pl -target TAIR10_chr_all.fas -t=dna -q=rna; rm control_re
p2_R1.fq
(... deleted)
```

TODO: need parameter description here.

The commands seem OK. So we add `system $cmd` for actual executing them.
```
ls src/*.fq.gz | perl -ne 'chomp; /.+\/(.+?)\./; $cmd="gzip -dc $_ > $1.fq; Mapping.pl -split 4 x $1.fq $1.blat.bam MappingBlat.pl -target TAIR10_chr_all.fas -t=dna -q=rna; rm $1.fq"; print "\nCMD: $cmd\n"; system $cmd'
```

Output BAM files of the last command would be `*.blat.bam`.
```
Singularity> ls -l *.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20462399 Sep 21 13:54 control_rep1_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20337560 Sep 21 13:54 control_rep1_R2.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20980407 Sep 21 13:54 control_rep2_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20874369 Sep 21 13:55 control_rep2_R2.blat.bam
-rwxrwxrwx+ 1 wdlin R418 19785876 Sep 21 13:55 control_rep4_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 19655968 Sep 21 13:55 control_rep4_R2.blat.bam
-rwxrwxrwx+ 1 wdlin R418 19713120 Sep 21 13:55 treatment_rep5_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 19576079 Sep 21 13:56 treatment_rep5_R2.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20835627 Sep 21 13:56 treatment_rep7_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20736819 Sep 21 13:56 treatment_rep7_R2.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20219439 Sep 21 13:57 treatment_rep9_R1.blat.bam
-rwxrwxrwx+ 1 wdlin R418 20108811 Sep 21 13:57 treatment_rep9_R2.blat.bam
```

As they are read alignment files of read1 and read2 separately, (not finished)
```
Singularity> ls *.blat.bam | perl -ne 'chomp; /(.+?)_R\d\./; push @{$hash{$1}},$_; if(eof){ for $key (sort keys %hash){ $cmd="samtools merge -fn /dev/stdout @{$hash{$key}} | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo $key.merged.bam -T TAIR10_chr_all.fas /dev/stdin"; print "\nCMD: $cmd\n"; system $cmd } }'

CMD: samtools merge -fn /dev/stdout control_rep1_R1.blat.bam control_rep1_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo control_rep1.merged.bam -T TAIR10_chr_all.fas /dev/stdin

CMD: samtools merge -fn /dev/stdout control_rep2_R1.blat.bam control_rep2_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo control_rep2.merged.bam -T TAIR10_chr_all.fas /dev/stdin

CMD: samtools merge -fn /dev/stdout control_rep4_R1.blat.bam control_rep4_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo control_rep4.merged.bam -T TAIR10_chr_all.fas /dev/stdin

CMD: samtools merge -fn /dev/stdout treatment_rep5_R1.blat.bam treatment_rep5_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo treatment_rep5.merged.bam -T TAIR10_chr_all.fas /dev/stdin

CMD: samtools merge -fn /dev/stdout treatment_rep7_R1.blat.bam treatment_rep7_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo treatment_rep7.merged.bam -T TAIR10_chr_all.fas /dev/stdin

CMD: samtools merge -fn /dev/stdout treatment_rep9_R1.blat.bam treatment_rep9_R2.blat.bam | samtools view /dev/stdin | SamReverse.pl _1 | samtools view -Sbo treatment_rep9.merged.bam -T TAIR10_chr_all.fas /dev/stdin
```

```
Singularity> ls -l *.merged.bam
-rwxrwxrwx+ 1 wdlin R418 37814246 Sep 21 14:40 control_rep1.merged.bam
-rwxrwxrwx+ 1 wdlin R418 38837864 Sep 21 14:41 control_rep2.merged.bam
-rwxrwxrwx+ 1 wdlin R418 36550185 Sep 21 14:41 control_rep4.merged.bam
-rwxrwxrwx+ 1 wdlin R418 36394270 Sep 21 14:41 treatment_rep5.merged.bam
-rwxrwxrwx+ 1 wdlin R418 38590467 Sep 21 14:41 treatment_rep7.merged.bam
-rwxrwxrwx+ 1 wdlin R418 37492552 Sep 21 14:41 treatment_rep9.merged.bam
```

## 2. Extract gene-exon coordinates

## 3. Compute basic numbers, separate biological replicates

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+?)\./; push @{$hash{$1}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.RPKMComputer -GFF tair10.strand.cgff -model tair10.strand.model -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -O $key -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.RPKMComputer -GFF tair10.strand.cgff -model tair10.strand.model -M SAMun control_rep1.merged.bam -O control_rep1 -ID 0.90 -direction true
(deleted...)
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+?)\./; push @{$hash{$1}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.ExonCounter -GFF tair10.strand.cgff -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -O $key -intronic true -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.ExonCounter -GFF tair10.strand.cgff -M SAMun control_rep1.merged.bam -O control_rep1 -intronic true -ID 0.90 -direction true
(deleted...)
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+?)\./; push @{$hash{$1}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.GeneCoverageArray -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -GFF tair10.strand.cgff -exon false -O $key -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.GeneCoverageArray -M SAMun control_rep1.merged.bam -GFF tair10.strand.cgff -exon false -O control_rep1 -ID 0.90 -direction true
(deleted...)
```

```
Singularity> ls *.geneCoverage | perl -ne 'chomp; $cmd="intronExonAvgDepth.pl tair10.strand.cgff $_"; print "\nCMD: $cmd\n"; system $cmd'

CMD: intronExonAvgDepth.pl tair10.strand.cgff control_rep1.geneCoverage
(deleted...)
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+?)\./; push @{$hash{$1}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.FineSpliceCounter -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -GFF tair10.strand.cgff -model tair10.strand.model -O $key -ID 0.90 -J 2 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.FineSpliceCounter -M SAMun control_rep1.merged.bam -GFF tair10.strand.cgff -model tair10.strand.model -O control_rep1 -ID 0.90 -J 2 -direction true
(deleted...)
```

## 3. Compute basic numbers, merged biological replicates

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+)_rep/; $key=$1; push @{$hash{$key}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.RPKMComputer -GFF tair10.strand.cgff -model tair10.strand.model -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -O $key -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.RPKMComputer -GFF tair10.strand.cgff -model tair10.strand.model -M SAMun control_rep1.merged.bam -M SAMun control_rep2.merged.bam -M SAMun control_rep4.merged.bam -O control -ID 0.90 -direction true

CMD: java rnaseq.RPKMComputer -GFF tair10.strand.cgff -model tair10.strand.model -M SAMun treatment_rep5.merged.bam -M SAMun treatment_rep7.merged.bam -M SAMun treatment_rep9.merged.bam -O treatment -ID 0.90 -direction true
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+)_rep/; $key=$1; push @{$hash{$key}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.ExonCounter -GFF tair10.strand.cgff -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -O $key -intronic true -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.ExonCounter -GFF tair10.strand.cgff -M SAMun control_rep1.merged.bam -M SAMun control_rep2.merged.bam -M SAMun control_rep4.merged.bam -O control -intronic true -ID 0.90 -direction true

CMD: java rnaseq.ExonCounter -GFF tair10.strand.cgff -M SAMun treatment_rep5.merged.bam -M SAMun treatment_rep7.merged.bam -M SAMun treatment_rep9.merged.bam -O treatment -intronic true -ID 0.90 -direction true
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+)_rep/; $key=$1; push @{$hash{$key}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.GeneCoverageArray -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -GFF tair10.strand.cgff -exon false -O $key -ID 0.90 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.GeneCoverageArray -M SAMun control_rep1.merged.bam -M SAMun control_rep2.merged.bam -M SAMun control_rep4.merged.bam -GFF tair10.strand.cgff -exon false -O control -ID 0.90 -direction true

CMD: java rnaseq.GeneCoverageArray -M SAMun treatment_rep5.merged.bam -M SAMun treatment_rep7.merged.bam -M SAMun treatment_rep9.merged.bam -GFF tair10.strand.cgff -exon false -O treatment -ID 0.90 -direction true
```

```
Singularity> ls *.geneCoverage | perl -ne 'print if not /.+_.+\./' | perl -ne 'chomp; $cmd="intronExonAvgDepth.pl tair10.strand.cgff $_"; print "\nCMD: $cmd\n"; system $cmd'

CMD: intronExonAvgDepth.pl tair10.strand.cgff control.geneCoverage

CMD: intronExonAvgDepth.pl tair10.strand.cgff treatment.geneCoverage
```

```
Singularity> ls *.merged.bam | perl -ne 'chomp; /(.+)_rep/; $key=$1; push @{$hash{$key}},$_; if(eof){ for $key (sort keys %hash){ $cmd="java rnaseq.FineSpliceCounter -M SAMun ".join(" -M SAMun ",@{$hash{$key}})." -GFF tair10.strand.cgff -model tair10.strand.model -O $key -ID 0.90 -J 2 -direction true"; print "\nCMD: $cmd\n"; system "$cmd >> $key.log" } }'

CMD: java rnaseq.FineSpliceCounter -M SAMun control_rep1.merged.bam -M SAMun control_rep2.merged.bam -M SAMun control_rep4.merged.bam -GFF tair10.strand.cgff -model tair10.strand.model -O control -ID 0.90 -J 2 -direction true

CMD: java rnaseq.FineSpliceCounter -M SAMun treatment_rep5.merged.bam -M SAMun treatment_rep7.merged.bam -M SAMun treatment_rep9.merged.bam -GFF tair10.strand.cgff -model tair10.strand.model -O treatment -ID 0.90 -J 2 -direction true
```

## 4. Alternative-splicing event comparison based on read counts of merged samples

## 5. Alternative-splicing evnet comparison based on ratios with respects to biological replicates
