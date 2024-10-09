# Bioinformatics Exome Analysis Assignment

This repository contains scripts and workflow for downloading, processing, aligning, and analyzing paired-end exome sequencing data from SRA. Steps include quality control, alignment to the human genome, variant calling using GATK, and variant annotation.


##  Download Paired-end FASTQ Files

[Sample 1 was download from this link](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR30834366
)

## Pre-process the Paired-end FASTQ Files

SRA Toolkit was used to download paired-end FASTQ files from the Sequence Read Archive (SRA) and split the forward and reverse read from the paired-end sequence. Here is the link to download the tool:
[SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
Here while executing the code replace 'SRRxxxxxx' with the SRR id of the respective sample.

```bash
prefetch SRRxxxxxx
fastq-dump --split-files SRRxxxxxx
```

After downloading the files, rename them as sample1_R1.fastq, sample1_R2.fastq, sample2_R1.fastq, sample2_R2.fastq

## 1. To calculate the fastq statistics: 
For this step seqkit tool was used. Here is the link to download the tool:
[SeqKit](https://bioinf.shenwei.me/seqkit/download/)
```bash
seqkit stats sample1_R1.fastq sample1_R2.fastq > sample1_stats.txt
seqkit stats sample2_R1.fastq sample2_R2.fastq > sample2_stats.txt

```

## 2. To perform basic QC checks on both files use FastQC:
For this step fastqc tool was used.  Here is the link to download the tool:
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
fastqc sample1_R1.fastq sample1_R2.fastq sample2_R1.fastq sample2_R2.fastq -o qc_results/

```

## 3. Now Download human genome fasta file from UCSC:
[Human GENOME FASTA FILE](https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED) and unzip it and rename it as 'hg38.fasta'



## 4. To align both sets of fastq files onto the human genome BWA is used:
For this step BWA-MEM tool was used.  Here is the link to download the tool:
[BWA-MEM](https://github.com/lh3/bwa)
First the human genome fasta file is  indexed.
```bash
bwa index hg38.fasta

```
And then to align both the files:
```bash
bwa mem -t 8 hg38.fasta sample1_R1.fastq sample1_R2.fastq > sample1_aln.sam
bwa mem -t 8 hg38.fasta sample2_R1.fastq sample2_R2.fastq > sample2_aln.sam

```

## 5. To calculate alignment statistics and store them in a file:
For this step SAMtools was used. Here is the link to download the tool:
[SAMtools](https://www.htslib.org/download/)
Use SAMtools to convert SAM to BAM and calculate alignment stats:
```bash
samtools view -S -b sample1_aln.sam > sample1_aln.bam
samtools flagstat sample1_aln.bam > sample1_alignment_stats.txt
samtools view -S -b sample2_aln.sam > sample2_aln.bam
samtools flagstat sample2_aln.bam > sample2_alignment_stats.txt


```

## 6. To perform mark duplicates and to store mark duplicate statistics:
For this step use Picard tools to mark duplicates. Here is the link to download the tool:
[Picard](https://www.htslib.org/download/)
```bash
picard MarkDuplicates I=sample1_aln.bam O=sample1_dedup.bam M=sample1_dedup_metrics.txt
picard MarkDuplicates I=sample2_aln.bam O=sample2_dedup.bam M=sample2_dedup_metrics.txt



```

## 7. To call variants using standard GATK best practices:
For this step use GATK to call variants. Here is the link to download the tool:
[GATK](https://github.com/broadinstitute/gatk/releases)
```bash
gatk HaplotypeCaller -R hg38.fasta -I sample1_dedup.bam -O sample1_variants.vcf
gatk HaplotypeCaller -R hg38.fasta -I sample2_dedup.bam -O sample2_variants.vcf


```

## 8. To separate out SNPs and Indels into individual vcf files:
For this step use GATK to call variants. Here is the link to download the tool:
[GATK](https://github.com/broadinstitute/gatk/releases)
```bash
gatk SelectVariants -R hg38.fasta -V sample1_variants.vcf --select-type-to-include SNP -O sample1_snps.vcf
gatk SelectVariants -R hg38.fasta -V sample1_variants.vcf --select-type-to-include INDEL -O sample1_indels.vcf

gatk SelectVariants -R hg38.fasta -V sample2_variants.vcf --select-type-to-include SNP -O sample2_snps.vcf
gatk SelectVariants -R hg38.fasta -V sample2_variants.vcf --select-type-to-include INDEL -O sample2_indels.vcf


```

