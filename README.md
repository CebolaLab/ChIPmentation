# Analysis pipeline for ChIPmentaion data

A step-by-step analysis pipeline for ChIP-seq data using the ChIPmentaion protocol from the [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/).

Correspondence: hannah.maude12@imperial.ac.uk

## ** \*UNDER CONSTRUCTION\* **

The resources and references used to build this tutorial are found at the bottom, in the [resources](#resources) section.

## Table of Contents

- [Pre-alignment quality control (QC)](#pre-alignment-qc) 
- [Alignment](#alignment) 
- [Post-alignment QC](#post-alignment-qc) - filter, check library complexity and format for peak calling
- [Peak Calling](#peak-calling)
- [Peak Calling QC and differential accessibility (DA) analysis](#peak-QC-and-DA)
- [Visualisation](#visualisation)
- Functional analysis & Motif Discovery


### Pre-alignment QC

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. As described by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```
fastqc <sample>_R1.fastq.gz -d . -o .

fastqc <sample>_R2.fastq.gz -d . -o .
```

### Adapter trimming 

Adapters and low quality reads/bases should be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). Adapter contamination can be seen in the fastqc report:

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/adapters.png" width="800">

For this pipeline, fastp is used to remove adapter sequences. The minimum fragment length is set at 20.

```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 20 -j <sample>.fastp.json -h <sample>.fastp.html
```

The output of fastp includes a html report, part of which is shown below. This presents the total number of reads before and after filtering, including the % of high quality (Q30) bases. The report also shows the main causes of read removal. In the example below, 1.9% of reads were removed because they were shorter than the minimum read length specified above by the -l argument (35bp).

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastp-html.png" width="600">

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This pipeline will use bowtie2 to align reads to the hg19 reference genome. If the user is aligning to the more recent GRCh38 release, it is recommended to remove alternative contigs, otherwise reads may not map uniquely and will be assigned a low quality score. Suggested guidelines for preparing the GRCh38 genome are discussed in [this tutorial](https://www.biostars.org/p/342482/). If the user selects an alternative alignment tool, such as bwa, they are referred to [this blog post](https://www.acgt.me/?offset=1426809676847) which discusses the resulting differences in alignment quality scores.

#### Bowtie2 alignment

The `local` parameter is used to 'soft clip' the end of reads to allow the best possible alignment, including any remaining adapter sequences (e.g. 1 or 2bp).  By using the `--no-mixed` and `--no-discordant` parameters, reads will only be aligned if both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired, which is a common post-alignment QC step). The `-I 35` and `-X 700` require fragments to be greater than 35bp and less than 700bp. The user can adjust these depending on the experiment/sequencing protocol (see the fastp html report for a plot of the estimated insert sizes). The maximum fragment length of 700bp prevents reads from aligning incorrectly outside the expected fragment range. 

***Which reference genome to use?*** See [this](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) discussion on which reference genome to use. The recommended downloads for both hg19/b37 and GRCh38 are shown below.


```bash
#to download the build 37 (hg19) reference genome
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

#to download the more recent GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Bowtie2 should be used to create the reference genome index files (see the bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)). After the index files have been generated, align the trimmed fastq files to the genome (here using hg19/b37). If a sample has been sequenced across multiple lanes, assuming the samples are balanced and there are no batch effects, bowtie2 can access multiple files as a comma seperated list:

```bash
#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx=/path/to/reference-genome

#Run the bowtie2 alignment and output a bam alignment file
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 20 -X 700 -x $bt2idx/human_g1k_v37.fasta -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam

#If your sample was sequenced across multiple lanes, in this case lane 1 (L002) and lane 2 (L003):
#bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 35 -X 700 -x $bt2idx/human_g1k_v37.fasta -1 <sample>_L001_R1.trimmed.fastq.gz,<sample>_L002_R1.trimmed.fastq.gz -2 <sample>_L001_R2.trimmed.fastq.gz,<sample>_L002_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam
```

The output `bam` file should be sorted and indexed prior to the downstream analysis:

```bash
#Sort the output bam file by coordinate
samtools sort <sample>.bam -o <sample>_sorted.bam 

#Generate an index file
samtools index <sample>_sorted.bam
```

## Post-alignment QC

The post-alignment QC involves several steps:

- [Remove duplicates & low-quality alignments](#remove-duplicates-&-low-quality-alignments) (including non-uniquely mapped reads)
- [Calculate library complexity and QC](#calculate-library-complexity-and-QC)
- [Remove ENCODE blacklist regions](#remove-encode-blacklist-regions)
- [Shift read coordinates](#shift-read-coordinates)

### Removed unmapped reads, multi-mapped reads and duplicates and estimate library complexity

*A note on sam file flags:* the output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which can map to more than one position are assigned a low quality scores. The user can assess the proportion of uniquely mapped reads (prior to the filtering step above) using `samtools -view -q 30 -c <sample_sorted.bam` (divide this number of reads by 2 to calculate the number of DNA fragments). In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). *A low % of uniquely mapped reads* may result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). 

The following command will remove unmapped reads and secondary alignments `samtools fixmate -r` and mark and remove duplicate reads `samtools markdup -r`. The `samtools markdup` will also estimate library complexity, with the output saved in the `<sample>.markdup.stats` file.

```bash
samtools sort -o <sample>_sorted.bam - | samtools fixmate -rcm - - | samtools sort - | samtools markdup -d 100 -r -f <sample>.markdup.stats - <sample>.rmdup.bam

samtools index <sample>.rmdup.bam
```

*A note on multi-mapping*: here, reads which align to more than one position have been removed (through the `samtools fixmate -r` option). Some users may opt to retain these 'multi-mapped reads', especially if single-end data is beign used. Removing multi-mapped reads can result in the loss of biologically informative reads (false negatives), but retaining them can lead to false potivies. The choice will depend on the study design. It is recommended by the Harvard hbctraining tutorial to *remove* multi-mapped reads to increase confidence and reproducibility:

*Sam file flags*: the read identity as a PCR duplicate, or uniquely mapped read is stored in the sam/bam file 'flag'. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevented the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.


### Remove ENCODE blacklist regions

The [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist/), most recently reported by [Amemiya et al. (2019)](https://www.nature.com/articles/s41598-019-45839-z) are defined as 'a comprehensive set of regions in the human, mouse, worm, and fly genomes that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment.' These problematic regions should be removed before further analysis. Download the blacklist files for your chosen reference genome from the [Boyle Lab github repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). Details regarding the identification of blacklist regions are reported [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist-README.pdf).

```bash
bedtools intersect -nonamecheck -v -abam <sample>.rmdup.bam -b hg19-blacklist.v2.bed > <sample>.blacklist-filtered.bam
```

### Shift read coordinates

An optional step in analysing data generated using the Tn5 transposase (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion, introducted as repair of the transposase-induced nick introduces a 9bp insertion. Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand should be offset by -5bp. For references, see the first ATAC-seq paper by [Buenrostro et al., (2013)](https://www.nature.com/articles/nmeth.2688) and the analysis by [Adey et al., (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119) which showed this insertion bias. Shifting coordinates is only really important if single-base resolution is required, for example in the analysis of transcription factor motifs in ATAC-seq peak footprints. Be aware that some tools do this shifting themselves (so double check manuals!).
 

We can use the `bedtools` command  `alignmentSieve`.


```bash
#The user can set the preferred number of processors 
alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName hg19-blacklist.v2.bed --bam <sample>.blacklist-filtered.bam -o <sample>.tmp.bam

#Sort and index the bam file
#Set the number of preferred threads with the -@ option
samtools sort -@ 8 -O bam -o <sample>.shifted.bam <sample>.tmp.bam
samtools index -@ 8 <sample>.shifted.bam

rm <sample>.tmp.bam
```

## Peak calling

The ChIP-seq peaks, either of histone marks or protein binding, will be called using the [MACS2](https://pypi.org/project/MACS2/) algorithm. 

```bash
macs2 callpeak -t <sample>.shifted.bam
```


MACSs - broad for histone marks, narrow for transcription factors. 

The ChIPmentaion paper uses MACS2 with:

- Peak-calling with MACS2. 
	Independently for biological replicates.
	Bandwideth of 200bp and matched IgG control as background. 

For both ChIP-seq and ChIPmentation data, MACS2 was run independently for biological replicates using a bandwidth of 200 bp and the matched IgG control as background. For broad histone marks (H3K27me3, H3K36me3) the “--broad”, “--nomodel”, “--extsize 73” and “--pvalue 1e-3” flags and arguments were provided. After ensuring consistency among replicates, downstream analysis was performed on peaks called from merged biological replicates in the same way as described. 


## Visualisation

The following code can be used to generate log<sub>10</sub> p-value tracks from the output of MACS peak calling. With ChIP-seq data, each sample should have a control input, which the data is normalised to. 


- Genome browser tracks - genomeCoverageBed command in BEDTools and bedGraphToBigWig tool (UCSC) was used to produce a bigWig file

## Peak Quality Control

To assess the quality of our peaks, we will use the *R* package ChIPQC as described in this [online tutorial](https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lessons/06_combine_chipQC_and_metrics.md) by the Harvard hbctraining. 

## Differential binding

[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)



## ChIPmentation analysis steps 

The original ChIPmentaion paper 'ChIPmentation: fast, robust, low-input ChIP-seq for histones and transcription factors' by Schmidl et al. 2015 can be found [here](https://www.nature.com/articles/nmeth.3542).

- Read trimming (Skewer20)
- Bowtie2 - hg19/GRCh37 assembly of the human genome using Bowtie2 (ref. 21) with the “--very-sensitive” and allowing for multimapper reads according to the aligner's default. 
- Duplicate reads were marked and removed using Picard. 
- Read shifting 
- Genome browser tracks - genomeCoverageBed command in BEDTools and bedGraphToBigWig tool (UCSC) was used to produce a bigWig file
- Peak-calling with MACS2. 
	Independently for biological replicates.
	Bandwideth of 200bp and matched IgG control as background. 

For both ChIP-seq and ChIPmentation data, MACS2 was run independently for biological replicates using a bandwidth of 200 bp and the matched IgG control as background. For broad histone marks (H3K27me3, H3K36me3) the “--broad”, “--nomodel”, “--extsize 73” and “--pvalue 1e-3” flags and arguments were provided. After ensuring consistency among replicates, downstream analysis was performed on peaks called from merged biological replicates in the same way as described. 



