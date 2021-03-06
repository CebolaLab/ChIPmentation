# Analysis pipeline for ChIPmentation data for histones

A step-by-step analysis pipeline for ChIP-seq data of histone modifications using the ChIPmentaion protocol [(Schmidl et al. 2015)](https://www.nature.com/articles/nmeth.3542) from the [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/).

Correspondence: hannah.maude12@imperial.ac.uk

## ** \*UNDER CONSTRUCTION\* **

The resources and references used to build this tutorial are found at the bottom, in the [resources](#resources) section.

## Experiment guidelines

According to [ENCODE](https://www.encodeproject.org/chip-seq/histone/), the current standards of ChIP-seq experiments are as follows:

- Each ChIP-seq experiment should have a corresponding input control experiment with matching run type, read length, and replicate structure. 
- For narrow-peak histone experiments, each replicate should have 20 million usable fragments.
- For broad-peak histone experiments, each replicate should have 45 million usable fragments.
- Quality control metrics are collected to determine library complexity, read depth, FRiP score, and reproducibility.

For uniformity of experiments, it is recommended that:

- The read length should be a minimum of 50 base pairs, though longer read lengths are encouraged; the pipeline can process read lengths as low as 25 base pairs. Sequencing may be paired- or single-ended.
- The sequencing platform used should be indicated.
- Replicates should match in terms of read length and run type. 
- Pipeline files are mapped to either the GRCh38 or mm10 sequences.

This GitHub repository contains an [excel spreadsheet](https://github.com/CebolaLab/ChIPmentation/blob/main/QC-template-ChIP-seq.xlsx) with the QC measures that should be generated and the recommended values from [ENCODE](https://www.encodeproject.org/chip-seq/histone/). The user will be prompted when to fill in values obtained during this pipeline by red boxes ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+). Blue boxes ![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) indicate where the expected output files are generated. 

The table below, lifted from [ENCODE](https://www.encodeproject.org/chip-seq/histone/), detailed the outputs:

<img src="https://github.com/CebolaLab/ChIPmentation/blob/main/Figures/outputs-ENCODE.png" width="800">


The main output files will therefore be:

For **individual** replicates:

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bed` peaks file (used mainly to assess reproducibility and as an input to call pooled peaks)

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bigWig` tracks for fold-enrichment

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bigWig` tracks for -log<sub>10</sub> *p*-value

For **pooled** replicates:

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bed` replicated peaks file pooled replicates

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bigWig` tracks for fold-enrichment

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) `bigWig` tracks for -log<sub>10</sub> *p*-value

In addition, this pipeline will cover differential binding analysis, functional analysis and motif discovery. 

## Table of Contents

- [Pre-alignment quality control (QC)](#pre-alignment-qc): trim reads
- [Alignment](#alignment): align to the reference genome
- [Post-alignment QC](#post-alignment-qc): filter, check library complexity and format for peak calling
- [Alignment visualisation](#alignment-visualisation): generate BigWig tracks from the aligned bam file
- [Peak calling](#peak-calling): call relaxed (per replicate) and replicated (pooled replicates) peaks and generate -log<sub>10</sub> *p*-value and fold-enrichment bigWig trakcs
- [Peak quality control](#peak-quality-control): calculate QC scores including fraction of reads in peaks (FRiP)
- [Differential binding analysis](#differntial-binding)
- Functional analysis & Motif Discovery

### Pre-alignment QC

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. As described by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end. The GC content may be expected to vary depending on the ChIP-seq target. For paired-end reads, run fastqc on both files, with the results output to the current directory. These fastQC reports can be combined into one summary report using [multiQC](https://multiqc.info/). 

```bash
#Generate fastQC reports for the forward and reverse reads
fastqc <sample>_R1.fastq.gz -d . -o .
fastqc <sample>_R2.fastq.gz -d . -o .

#Combine reports across samples using multiQC
multiqc *.html
```

The total number of DNA reads is given in the fastQC report under 'Total Sequences'. If using paired-end reads, the total sequences in one file (R1 *or* R2) is the number of DNA fragments. Sum the total reads in the read 1 (R1) *and* read 2 (R2) files to calculate the **total number of reads**. (Note that for paired-end reads, each DNA fragment has two reads since it is sequenced from both ends). If a sample has been sequenced across multiple lanes, sum together all of the 'total sequences' to calculate the total number of reads. 

To automatically extract the total number of sequences, run:

```bash
totalreads=$(unzip -c <sample>_fastqc.zip <sample>_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2)

echo $totalreads
#This number will be used again later so is saved as a variable 'totalreads'
```

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value:** input the total number of reads into the QC spreadsheet. 

### Adapter trimming 

Adapters and low quality reads/bases should be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). Adapter contamination can be seen in the fastqc report:

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/adapters.png" width="800">

For this pipeline, fastp is used to remove adapter sequences. The minimum fragment length is set at 20.

```bash
#Assuming paired-end reads:
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 20 -j <sample>.fastp.json -h <sample>.fastp.html

#If using single-end data, refer to the fastp manual. The --detect_adapter_for_pe argument should be removed and the adapter sequence should be provided.
```

The output of fastp includes a html report, part of which is shown below. This presents the total number of reads before and after filtering, including the % of high quality (Q30) bases. The report also shows the main causes of read removal. In the example below, 1.9% of reads were removed because they were shorter than the minimum read length specified above by the -l argument (35bp).

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastp-html.png" width="600">


To automatically calculate the % of reads removed during trimming:

```bash
printf %.2f `echo "(1-$(echo $(grep "total_reads" <sample>.fastp.json | head -n2 | cut -d : -f2 | cut -d , -f1 | paste - - | awk '{ print $2/$1 }')))*100" | bc -l`
```

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the % of reads removed during trimming to the QC spreadsheet. 

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This pipeline will use bowtie2 to align reads to the GRCh38 reference genome. If the user is aligning to the more recent GRCh38 release, it is recommended to remove alternative contigs, otherwise reads may not map uniquely and will be assigned a low quality score. Suggested guidelines for preparing the GRCh38 genome are discussed in [this tutorial](https://www.biostars.org/p/342482/). If the user selects an alternative alignment tool, such as bwa, they are referred to [this blog post](https://www.acgt.me/?offset=1426809676847) which discusses the resulting differences in alignment quality scores.

#### Bowtie2 alignment

The `local` parameter is used to 'soft clip' the end of reads to allow the best possible alignment, including any remaining adapter sequences (e.g. 1 or 2bp).  By using the `--no-mixed` and `--no-discordant` parameters, reads will only be aligned if both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired, which is a common post-alignment QC step). The `-I 35` and `-X 700` require fragments to be greater than 35bp and less than 700bp. The user can adjust these depending on the experiment/sequencing protocol (see the fastp html report for a plot of the estimated insert sizes). The maximum fragment length of 700bp prevents reads from aligning incorrectly outside the expected fragment range. 

***Which reference genome to use?*** See [this](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) discussion on which reference genome to use. The recommended downloads for both hg19/b37 and GRCh38 are shown below.


```bash
#to download the build 37 (hg19) reference genome
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

#to download the more recent GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Bowtie2 should be used to create the reference genome index files (see the bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)). After the index files have been generated, align the trimmed fastq files to the genome (here using GRCh38/hg38). If a sample has been sequenced across multiple lanes, assuming the samples are balanced and there are no batch effects, bowtie2 can access multiple files as a comma seperated list:

```bash
#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx=/path/to/reference-genome

#Run the bowtie2 alignment and output a bam alignment file
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 20 -X 700 -x $bt2idx/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam

#If your sample was sequenced across multiple lanes, in this case lane 1 (L002) and lane 2 (L003):
#bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 35 -X 700 -x $bt2idx/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -1 <sample>_L001_R1.trimmed.fastq.gz,<sample>_L002_R1.trimmed.fastq.gz -2 <sample>_L001_R2.trimmed.fastq.gz,<sample>_L002_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam
```

The output `bam` file should be sorted and indexed prior to the downstream analysis:

```bash
#Sort the output bam file by coordinate
samtools sort <sample>.bam -o <sample>_sorted.bam 

#Generate an index file
samtools index <sample>_sorted.bam

#Calculate the total number of mapped reads
samtools flagstat <sample>_sorted.bam > <sample>_sorted.flagstat
```

## Post-alignment QC

The post-alignment QC involves several steps:

- [Remove unmapped, multi-mapped and duplicates reads](#removed-unmapped-multi-mapped-and-duplicates-reads-and-estimate-library-complexity)
- [Remove ENCODE blacklist regions](#remove-encode-blacklist-regions)
- [Shift read coordinates](#shift-read-coordinates-optional)

### Remove unmapped, multi-mapped and duplicates reads

*A note on sam file flags:* the output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which can map to more than one position are assigned a low quality scores. The user can assess the proportion of uniquely mapped reads (prior to the filtering step above) using `samtools -view -q 30 -c <sample_sorted.bam` (divide this number of reads by 2 to calculate the number of DNA fragments). In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). A low % of uniquely mapped reads may result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). 

The following command will remove unmapped reads and secondary alignments `samtools fixmate -r` and mark and remove duplicate reads `samtools markdup -r`. The `samtools markdup` will also estimate library complexity, with the output saved in the `<sample>.markdup.stats` file.

```bash
#Sort by order, fixmate and remove unmapped/multi-mapped reads, sort by coordinate, mark and remove duplicates and estimate library complexity
samtools sort -n <sample>_sorted.bam | samtools fixmate -rcm - - | samtools sort - | samtools markdup -d 100 -r -f <sample>.markdup.stats - <sample>.rmdup.bam

#Index the resulting bam file
samtools index <sample>.rmdup.bam
```

From the `<sample>.markdup.stats` file:

Calculate the % of duplicates removed by dividing the **DUPLICATE TOTAL** by the number **READ**. The number READ is the total number of uniquely mapped reads, before duplicate removal. Note that multi-mapped and unmapped reads were removed prior in the `samtools fixmate -rcm` step, leaving only uniquely mapped reads (reads which map to one location in the genome).

**Calculate the % of (duplicate) reads removed:** 

```bash
printf %.2f $(echo $(echo -e `grep -e 'DUPLICATE TOTAL' -e 'READ' <sample>.markdup.stats | cut -d ':' -f2` | paste - - | awk '{ print $2/$1 }')*100 | bc -l) 
```
![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the % of duplicates into the QC spreadsheet. 

**Exctract the total number of uniquely mapped reads after duplicate removal:** 

```bash
totalunique=$(grep 'WRITTEN' <sample>.markdup.stats | cut -d ' ' -f 2)

echo $totalunique
```
![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the total number of uniquely mapped, non-duplicated reads into the QC spreadsheet.

**Extract the estimated library complexity:** 

```bash
grep 'ESTIMATED_LIBRARY_SIZE' <sample>.markdup.stats | cut -d ' ' -f 2
```

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the estimated library size to the QC spreadsheet.

**Calculate the non-redundant fraction (NRF)**

The NRF score is calculated as the number of distinct uniquely mapping reads (`$totalunique`) / total number of reads (`$totalreads`).

```bash
printf %.2f $(echo $totalunique/$totalreads | bc -l)
```

<img src="https://github.com/CebolaLab/ChIPmentation/blob/main/Figures/NRF-guidelines.png" width="400">


![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the NRF score into the QC spreadsheet.

*A note on multi-mapping*: here, reads which align to more than one position have been removed (through the `samtools fixmate -r` option). Some users may opt to retain these 'multi-mapped reads', especially if single-end data is beign used. Removing multi-mapped reads can result in the loss of biologically informative reads (false negatives), but retaining them can lead to false potivies. The choice will depend on the study design. It is recommended by the Harvard hbctraining tutorial to *remove* multi-mapped reads to increase confidence and reproducibility:

*Sam file flags*: the read identity as a PCR duplicate, or uniquely mapped read is stored in the sam/bam file 'flag'. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevented the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.


### Remove ENCODE blacklist regions

The [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist/), most recently reported by [Amemiya et al. (2019)](https://www.nature.com/articles/s41598-019-45839-z) are defined as 'a comprehensive set of regions in the human, mouse, worm, and fly genomes that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment.' These problematic regions should be removed before further analysis. Download the blacklist files for your chosen reference genome from the [Boyle Lab github repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). Details regarding the identification of blacklist regions are reported [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists).

```bash
#Remove reads within blacklist regions
bedtools intersect -nonamecheck -v -abam <sample>.rmdup.bam -b hg38-blacklist.v2.bed > <sample>.blacklist-filtered.bam
```

### Shift read coordinates [optional]

An optional step in analysing data generated using the Tn5 transposase (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion, introducted as repair of the transposase-induced nick introduces a 9bp insertion. Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand should be offset by -5bp [(Buenrostro et al., 2013; ](https://www.nature.com/articles/nmeth.2688)[Adey et al., 2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119). Shifting coordinates is only really important if single-base resolution is required, for example in the analysis of transcription factor motifs in ATAC-seq peak footprints. Be aware that some tools do this shifting themselves (so double check manuals!).
 

We can use the `bedtools` command  `alignmentSieve`.

```bash
#The user can set the preferred number of processors 
alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName hg38-blacklist.v2.bed --bam <sample>.blacklist-filtered.bam -o <sample>.tmp.bam

#Sort and index the bam file
#Set the number of preferred threads with the -@ option
samtools sort -O bam -o <sample>.shifted.bam <sample>.tmp.bam
samtools index <sample>.shifted.bam

rm <sample>.tmp.bam
```

## Alignment visualisation

The `<sample>.blacklist-filtered.bam` or `<sample>.shifted.bam` file can be converted to a `BigWig` file to visualise the alignment as a track in a genome browser, such as UCSC. For ChIP-seq, each sample is expected to have a control 'input' sample, which the aligned `bam` file may be normalised to. 

The input to `bamCoverage` (see below) requires the effective genome size to be estimated; a table is provided [at this link](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html). Select the appropriate value depending on the read length and reference genome. 

To visualise the input and data tracks:

```bash
#RPGC is reads per genome coverage
#2862010578 is the effective genome size for GRCh38 when using 150bp reads and including only regions which are uniquely mappable
#Edit <sample>.shifted.bam to <sample>.blacklist-filtered.bam if you did not run the coordinate shifting
bamCoverage --bam <sample>.shifted.bam -o <sample>.SeqDepthNorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2862010578 --ignoreForNormalization chrX --extendReads --blackListFileName hg38-blacklist.v2.bed

#Repeat for the input control
bamCoverage --bam <sample>.input.bam -o <sample>.input.SeqDepthNorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2862010578 --ignoreForNormalization chrX --extendReads --blackListFileName hg38-blacklist.v2.bed
```

## Peak calling

The ChIP-seq peaks, either of histone marks or protein binding, will be called using the [MACS2](https://pypi.org/project/MACS2/) algorithm. It is important to first know whether you want to call **broad** or **narrow** peaks. Typically, transcription factors form narrow peaks, although there are exceptions such as PolII which binds across the gene body and thus forms 'broad' peaks of binding. For histone marks, examples of narrow peaks include marks enriched at transcription state sites, wherease marks which mark heterochromatin may cover extensive regions and therefore form broad peaks. The table below, lifted from [ENCODE](https://www.encodeproject.org/chip-seq/histone/), shows the categories of histone-peak types. H3K9me3 is an exception as it is enriched in repetitive regions of the genome. MACS2 is typically more reliable for calling narrow peaks.

<img src="https://github.com/CebolaLab/ChIPmentation/blob/main/Figures/broad-vs-narrow-histones-ENCODE.png" width="600">

According to the [ENCODE guidelines](https://www.encodeproject.org/chip-seq/histone/), **narrow-peak** histone experiments should have **at least 20 million usable fragments**, while **broad-peak** histone experiments should have at least **45 million usable fragments**. 

A useful tutorial on how MACS2 calls peaks is provided [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Detailed information on the subcommands used within `macs2 callpeak` are provided [at this link](https://github.com/macs3-project/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands).

The following sections cover:

- [Call peaks for individual replicates](#call-peaks-for-individual-replicates) and generate bigWig tracks
- [Call peaks for pooled replicates](#call-peaks-for-pooled-replicates) and generate bigWig tracks
- [Extract replicated peaks](#extract-replicated-peaks)


### Call peaks for individual replicates

To call **narrow** peaks for individual replicates:

```bash
#Call peaks
macs2 callpeak -t <sample>.shifted.bam -c <input>.bam -f BAM -B --broad --keep-dup all --cutoff-analysis -g 2862010578 -n <sample> --outdir <sample>.macs2 2> <sample>_macs2.log
```

Note, for broad histone marks (H3K27me3, H3K36me3) the parameters used in the original ChIPmentation paper by [Schmidl et al. (2015)](https://www.nature.com/articles/nmeth.3542) are `--broad --nomodel --extsize 73 --pvalue 1e-3`.

The output files:

- `_peaks.narrowPeak`: a BED6+4 file detailing the peak locations, along with the peak summits, *p*-value and *q*-values 
- `_peaks.xls`: a tabular file containing addition information, such as pileup and fold-enrichment.
- `_summits.bed`: the locations of the summits for all peaks. 
- `_model.R`: an R script used to plot a PDF model based on your data and cross-correlation plot
- `_control_lambda.bdg`: bedGraph format for the input sample
- `_treat_pileup.bdg`: bedGraph format for the treatment sample

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: the `<sample>_peaks.narrowPeak` is a `bed` file containing the peak information for the INDIVIDUAL replicate\*.

\*The `<sample>_peaks.narrowPeak` can be uploaded and visualised via a genome browser such as UCSC. The `bed` file of peak calls is referred to at this stage as 'relaxed' peak calls, since they are called for individual replicates. This `bed` file is therefore most useful as an input to the next stage: calling replicated peaks.

The total number of peaks can be obtained using `wc -l <sample>_peaks.narrowPeak`. 

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the total number of peaks into the QC spreadsheet.


From these output files, we will generate:

1. A `bigWig` track of the fold-enrichment (treatment over the background)
2. A `bigWig` track of the -log<sub>10</sub> *p*-value (treatment over the background)

**1. Fold-enrichment bigWig**

The following commands require an input file detailing the chromosome sizes. Use the UCSC tool `fetchChromSizes` (install via [conda](https://anaconda.org/bioconda/ucsc-fetchchromsizes)): `fetchChromSizes hg38 > hg38.chrom.sizes`

```bash
#Generate the fold-change bedGraph
macs2 bdgcmp -t <sample>_treat_pileup.bdg -c <sample>_control_lambda.bdg -m FE -o <sample>_FE.bdg 

#Sort the bedGraph file and convert to bigWig
sort -k1,1 -k2,2n <sample>_FE.bdg > <sample>_FE.sorted.bdg

bedGraphToBigWig <sample>_FE.sorted.bdg hg38.chrom.sizes <sample>_macs2_FE.bw
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_macs2_FE.bw`

**2. -log<sub>10</sub> *p*-value bigWig**

```bash
#Generate the p-value bedGraph
macs2 bdgcmp -t <sample>_treat_pileup.bdg -c <sample>_control_lambda.bdg -m ppois -o <sample>_ppois.bdg

#Sort the bedGraph file and convert it to bigWig using the hg38 chromosome sizes
sort -k1,1 -k2,2n <sample>_ppois.bdg > <sample>_ppois.sorted.bdg

bedGraphToBigWig <sample>_ppois.sorted.bdg hg38.chrom.sizes <sample>_macs2_pval.bw
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_macs2_pval.bw`

### Call peaks for pooled replicates

The step assumes that the ChIP-seq expriment includes *biological replicates* for each treated condition. Best practise requires a combined set of peaks for the pooled replicates to be called. Replicated peaks are defined as the peaks in the pooled replicates which are observed in at least two replicates. 

The following steps assume that there are two biological replicates, `rep1` and `rep2`. Any number of replicates can be added. First, check the correlation between the replicates using the UCSC tool wigCorrelate:

```bash
wigCorrelate <sample>_rep1_macs2_FE.bw <sample>_rep2_macs2_FE.bw
```
Assuming there is a satisfactory correlation, peaks should be called using the pooled replicates by including all `bam` files in the `macs2 callpeak` command:

```bash
macs2 callpeak -q 0.01 -t <sample>_rep1_shifted.bam <sample>_rep2_shifted.bam -c <input>_rep1.bam <input>_rep2.bam -f BAM -B --broad --keep-dup all --cutoff-analysis -g 2862010578 -n <sample>_pooled -outdir <sample>_pooled.macs2 2> <sample>_pooled.macs2/<sample>_pooled_macs2.log
```

The outut `<sample>_pooled_peaks.narrowPeak` file can be used to define the replicated peaks. `intersectBed` from [bedTools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) will be used.

#### Generate bigWigs

Generate the Fold-enrichment and -log<sub>10</sub> *p*-value bigWig files for the pooled replicates:

```bash
#Generate the fold-change bedGraph
macs2 bdgcmp -t <sample>_pooled_treat_pileup.bdg -c <sample>_pooled_control_lambda.bdg -m FE -o <sample>_pooled_FE.bdg 
sort -k1,1 -k2,2n <sample>_pooled_FE.bdg > <sample>_pooled_FE.sorted.bdg
bedGraphToBigWig <sample>_pooled_FE.sorted.bdg hg38.chrom.sizes <sample>_pooled_macs2_FE.bw

#Generate the p-value bedGraph
macs2 bdgcmp -t <sample>_pooled_treat_pileup.bdg -c <sample>_pooled_control_lambda.bdg -m ppois -o <sample>_pooled_ppois.bdg
sort -k1,1 -k2,2n <sample>_ppois.bdg > <sample>_ppois.sorted.bdg
bedGraphToBigWig <sample>_pooled_ppois.sorted.bdg hg38.chrom.sizes <sample>_pooled_macs2_pval.bw
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_macs2_FE.bw`

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_pooled_macs2_pval.bw`

### Extract replicated peaks

*The following code is adapted from the ENCODE pipeline.* As done by ENCODE, overlapping peaks should overlap by at least 50% for either of the two peaks. First, the pooled peaks will be subsetted for those which overlap replicate 1, then further subsetted for those which also overlap replicate 2. 

**For narrow peaks:** 

```bash
#Identify peaks from the POOLED replicates which are in BOTH replicate 1 and replicate 2

#First extract pooled peaks which are in replicate 1
intersectBed -wo -a <sample>_pooled.narrowPeak -b <sample>_rep2_peaks.narrowPeak  | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$13-$12; if(($21/s1 > 0.5) || ($21/s2 > 0.5)) {print $0}}' | cut -f 1-10 > tmp.bed

#Next, take these peaks and extract the ones which overlap with replicate 2
intersectBed -wo -a tmp.bed -b <sample>_rep1_peaks.narrowPeak | awk_command | cut_command > tmp_pooled | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$13-$12; if(($21/s1 > 0.5) || ($21/s2 > 0.5)) {print $0}}' | cut -f 1-10 > replicated_narrowPeaks.bed
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `replicated_narrowPeaks.bed`

**For broad peaks:** 

```bash
#Identify peaks from the POOLED replicates which are in BOTH replicate 1 and replicate 2

#First extract pooled peaks which are in replicate 1
intersectBed -wa -a <sample>_pooled.broadPeak -b <sample>_rep2_peaks.broadPeak  | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$12-$11; if(($19/s1 > 0.5) || ($19/s2 > 0.5)) {print $0}}' | cut -f 1-9 > tmp.bed

#Next, take these peaks and extract the ones which overlap with replicate 2
intersectBed -wa -a tmp.bed -b <sample>_rep1_peaks.broadPeak | awk_command | cut_command > tmp_pooled | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$12-$11; if(($19/s1 > 0.5) || ($19/s2 > 0.5)) {print $0}}' | cut -f 1-9 > replicated_broadPeak.bed
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `replicated_broadPeak.bed`


## Peak quality control

To assess the quality of our peaks, we will use the *R* package ChIPQC as described in this [online tutorial](https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lessons/06_combine_chipQC_and_metrics.md) by the Harvard Chan Bioinformatics Core.


## Differential binding

Differential binding will be assessed using [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html). (See this [paper](https://www.nature.com/articles/s41598-020-66998-4)). DiffBind makes use of DESeq2 to assess differential binding. 

To create a conda environment with the required packages:

```bash
conda create -n DiffBind r-essentials r-base r-tidyverse
conda install -n DiffBind -c bioconda bioconductor-diffbind 
conda install -n DiffBind -c bioconda bioconductor-genomeinfodbdata
```

Or, run the following from R:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DiffBind")
BiocManager::install("GenomeInfoDbData")

install.packages("tidyverse") 
```

Create a csv file containing 



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

# Resources

- ENCODE define replicated peaks [script](https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/overlap_peaks/src/overlap_peaks.py)
- ENCODE [terms and definitions](https://www.encodeproject.org/data-standards/terms)
- ChIP-seq tutorials from the [Harvard Chan Bioinformatics Core (HBC)](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/)
