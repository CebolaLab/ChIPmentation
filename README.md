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

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```
fastqc <sample>_R1.fastq.gz -d . -o .

fastqc <sample>_R2.fastq.gz -d . -o .
```

### Adapter trimming 

Adapters and low quality reads/bases should be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). Adapter contamination can be seen in the fastqc report:

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/adapters.png" width="800">

For this pipeline, fastp is used to remove adapter sequences. The minimum fragment length is set at 35, since short ATAC-seq fragments can be observed if the transposase cuts adjacent nucleosome-free DNA. 

```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 35 -j <sample>.fastp.json -h <sample>.fastp.html
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

Bowtie2 should be used to create the reference genome index files (see the bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)). After the index files have been generated, align the trimmed ATAC-seq fastq files to the genome (here using hg19/b37). If a sample has been sequenced across multiple lanes, assuming the samples are balanced and there are no batch effects, bowtie2 can access multiple files as a comma seperated list:

```bash
#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx=/path/to/reference-genome

#Run the bowtie2 alignment and output a bam alignment file
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 35 -X 700 -x $bt2idx/human_g1k_v37.fasta -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam

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

### Mark duplicates 

To mark duplicate reads and view the % of duplicates:

```bash
picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

#View the % of duplicates
head -n 8 <sample>.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
```

### Remove duplicates & low-quality alignments 

The output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which are uniquely mapped are assigned a high alignment quality score and one genomic position. If reads can map to more than one location, Bowtie2 reports one position and assigns a low quality score. The proportion of uniquely mapped reads can be assessed. In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). Secondly, the 'flag' reports information such as whether the read is mapped as a pair or is a PCR duplicate. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevent the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.

**Multi-mapping:** the user should decide whether or not to retain reads which have multi-mapped, i.e. aligned to more than one position in the reference genome. When using paired-end data, it may be the case that one read aligns to a repetitive region (and therefore can map elsewhere), while the mate aligns to a unique sequence with a high quality. The bowtie2 parameters used above required reads to align within 50-700bp, so there should be no reads incorrectly aligned outside this distance. As such, the user may decide to keep multi-mapping reads on the assumption that they are likely to be mapped to the correct sequence, within the length of the DNA fragment. This may, however, cause incorrect alignments in extended repetitive regions where a read could map to multiple positions within the length of the DNA fragment. This should be minimised by the downstream removal of the [ENCODE Blacklisted regions](https://www.nature.com/articles/s41598-019-45839-z).

If a read is multi-mapped, it is assigned a low quality score by bowtie2. To view how many DNA reads align with a quality score >30, run the following (divide this number by 2 to calculate the # of DNA fragments):

```bash
samtools view -q 30 -c <sample>.marked.bam
```

A low % of uniquely mapped reads map result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). The following code uses the `sam/bam` flags to retain properly mapped pairs (`-f 2`) and to remove reads which fail the platform/vendor QC checks (`-F 512`), duplicate reads (`-F 1024`) and those which are unmapped (`-F 12`). The three flags to be removed can be combined into `-F 1548`, which will remove reads which meet any of the three individual flags

To ***retain*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 <sample>.rmChrM.bam | samtools sort -n -o <sample>.filtered.bam 
```

To ***remove*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 -q 30 <sample>.rmChrM.bam | samtools sort -o <sample>.filtered.bam
```

The output `bam` file, which is now sorted by name, should be indexed: 

```bash
samtools index <sample>.filtered.bam
```

### Remove ENCODE blacklist regions

The [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist/), most recently reported by [Amemiya et al. (2019)](https://www.nature.com/articles/s41598-019-45839-z) are defined as 'a comprehensive set of regions in the human, mouse, worm, and fly genomes that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment.' These problematic regions should be removed before further analysis. Download the blacklist files for your chosen reference genome from the [Boyle Lab github repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). Details regarding the identification of blacklist regions are reported [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist-README.pdf).

```bash
bedtools intersect -nonamecheck -v -abam <sample>.filtered.bam -b hg19-blacklist.v2.bed > <sample>.blacklist-filtered.bam
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


MACSs - broad for histone marks, narrow for transcription factors. 

The ChIPmentaion paper uses MACS2 with:

- Peak-calling with MACS2. 
	Independently for biological replicates.
	Bandwideth of 200bp and matched IgG control as background. 

For both ChIP-seq and ChIPmentation data, MACS2 was run independently for biological replicates using a bandwidth of 200 bp and the matched IgG control as background. For broad histone marks (H3K27me3, H3K36me3) the “--broad”, “--nomodel”, “--extsize 73” and “--pvalue 1e-3” flags and arguments were provided. After ensuring consistency among replicates, downstream analysis was performed on peaks called from merged biological replicates in the same way as described. 


## Visualisation

The following code can be used to generate <sub>10</sub> p-value tracks from the output of MACS peak calling. With ChIP-seq data, each sample should have a control input, with which the data is normalised to.

**Convert the bam file to a bed file**

```bash
bedtools bamtobed -i <chip>.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > <chip>.bed.gz
bedtools bamtobed -i <input>.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > <input>.bed.gz
```

Generate a pileup file, which contains the number of reads in each peak

```bash
estd=300
###  Extend ChIP sample to get ChIP coverage track                                                                                
macs2 pileup -i ${chip}.bed.gz -o ${chip}.pileup.bdg --extsize ${estd}
###  Extend the control read to both sides (-B option) using pileup function.                                                     
macs2 pileup -i ${input}.bed.gz -B --extsize 150 -o d_bg.bdg

###  Create a local background ( 1kb window)       
macs2 pileup -i ${input}.bed.gz -B --extsize 500 -o 1k_bg.bdg
### Normalize the 1kb noise by multiplying the values by d/slocal ( its 300/1000 )      
macs2 bdgopt -i 1k_bg.bdg -m multiply -p 0.3 -o 1k_bg_norm.bdg

###  Create a large local background ( 10 kb window)                                                                              
macs2 pileup -i ${input}.bed.gz -B --extsize 5000 -o 10k_bg.bdg
### Normalize the 10kb noise by multiplying the values by d/slocal ( its 300/10000 )                                              
macs2 bdgopt -i 10k_bg.bdg -m multiply -p 0.03 -o 10k_bg_norm.bdg

### Compute the maximum bias for each genomic location.                                                                           
macs2 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
macs2 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
ctrl_reads=`zcat ${input}.bed.gz | wc -l`

###  The whole genome background can be calculated as the number of control reads/genome size*fragment length                     
genome_background=`echo "123" | awk -v ctrl_reads=$ctrl_reads '{ print (ctrl_reads*300)/2700000000 }'`
### Compute the genome wide background                                                                                            
macs2 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p ${genome_background} -o local_bias_raw.bdg
chip_reads=`zcat ${chip}.bed.gz | wc -l `
###   Create a temp copy of chip and local bias sample                                                                            
cp ${chip}.pileup.bdg scaled_treat.pileup.bdg
cp local_bias_raw.bdg local_lambda.bdg

###  Scale the ChIP and control to the same sequencing depth.                                                                     
### If the control has more reads, it needs to be scaled down by multiplying the ratio of Chip/control or vice versa              
echo "foo" | awk -v ctrl=${ctrl_reads} -v chip=${chip_reads} -v chip_sample=${chip}.pileup.bdg '{                                 
        if ( ctrl > chip ) {                                                                                                      
                scale=chip/ctrl;                                                                                                  
                print "macs2 bdgopt -i local_bias_raw.bdg -m multiply -p " scale " -o local_lambda.bdg"                           
        }                                                                                                                         
        else {                                                                                                                    
                scale=ctrl/chip;                                                                                                  
                print "macs2 bdgopt -i " chip_sample " -m multiply -p " scale " -o scaled_treat.pileup.bdg"                       
        }                                                                                                                         
 }' | parallel

mv scaled_treat.pileup.bdg ${chip}.pileup.bdg
mv local_lambda.bdg local_bias_raw.bdg
macs2 bdgcmp -t ${chip}.pileup.bdg -c local_bias_raw.bdg -m ppois -o ${chip}.pileup_pvalue.bdg
macs2 bdgcmp -t ${chip}.pileup.bdg -c local_bias_raw.bdg -m qpois  -o ${chip}.pileup_qvalue.bdg
cp local_bias_raw.bdg ${chip}.pileup_pvalue.bdg ${chip}.pileup_qvalue.bdg ${chip}.pileup.bdg ${wd}
```


- Genome browser tracks - genomeCoverageBed command in BEDTools and bedGraphToBigWig tool (UCSC) was used to produce a bigWig file



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