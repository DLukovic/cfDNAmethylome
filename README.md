# cfDNAmethylome
Workflow for Twist targeted methylation sequencing of cfDNA 

md5 (Message Digest Algorithm 5)
==================================
+ md5 creates md5 hash
+ To get md5 hash for multiple files in a directory with .fa files, use a wildcard to match all files as follows:

		md5 /check_this_dir/* > files.md5 

md5sum
=============
+ use the md5sum command to verify MD5 hashes for multiple files

		$ md5sum * | tee md5sum.txt 

md5sum compute the hash values of the files in the current directory `*` and
`tee`is used to simultaneously display the outuput on stdout and write the otutput to `md5sum.txt` file


+ run md5sum -c command with the text file as argument

		$ md5sum -c mdsum.txt
		
-`c` option of `md5sum` can accept a text file as an input
Download reference sequences
============================
* Download Twist Human Methylome Panel hg38 from Twist Bioscience website 	
<https://www.twistbioscience.com/resources/data-files/twist-human-methylome-panel-target-bed-file>
* Download hg38/GRCh38 FASTA human reference genome via wget 
* Terminal: 

		$ brew install wget
		$ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
* Decompressed the zipped file

		$ gzip -d hg38.fa.gz

Trim adapters 	
=============================
Trim adapters using TrimGalore!

Dependencies: cutadapt - export to your PATH

search for cutadapt: 
		$ echo $PATH | tr ":" "\n" | grep "cutadapt" 
write into .zprofile file:

		# Setting PATH for CUTADAPT
		  export PATH=$PATH:$HOME/.local/bin
		 
run in terminal:

		$ source ~/.zprofile

Use Snakemake Workflow to trim adapter from paired fastqfiles.

Trimming individual sample:

 
	
	
 	~/TrimGalore-0.6.10/trim_galore --cores 8 --output_dir output --fastqc --2colour 20 --paired DE17NGSLABD100901_1.fq.gz DE17NGSLABD100901_2.fq.gz 


Generate reference genome index and dictionary file (bwa, samtools, Picard)
==============================================================================
* Make an index of the reference genome in FASTA format

		bwa index -a bwtsw hg38.fa
	
		
bwa needs to contstruct an index, -a for BWT index, bwtsv is optiom implemented in BWT-SW, whorks with human whole genome

* Create an index for your FASTA.file (must be unziped)

		samtools faidx hg38.fa > hg38.fa.fai
		

Reference Sequence Directory 
---------------------------

		java -Xmx4g -Xms4g -jar picard.jar CreateSequenceDictionary REFERENCE=/Users/gyongyosilab/Documents/ccfDNA\ Projekt/Reference\ Sequences/hg38.fa  	OUTPUT=/Users/gyongyosilab/Documents/ccfDNA\ Projekt/Reference\ Sequences/hg38.dict


Alignment
==========
install bwameth.py in conda environment

		conda install bioconda::bwameth

create environment for alignment

		conda create --name methylenv

create indices
		
   		bwameth.py index /Users/gyongyosilab/Documents/ccfDNA_Projekt/Reference_Sequences/hg38.fa 

Apply following options:

	@RG: Read group header.
	ID:foo: Identifier for the read group.
	SM:bar: Sample name. These annotations help downstream tools distinguish samples or batches of reads.
	
Execute bwameethalignment using *snakemake*  [workflow](Rules/bwameth/Snakefile) management system.

Filter BAM files
==================
Next steps requires presence of *--read-group* in BAM files. For filtering and sorting .bam files use sambamba package in following steps
	
		sambamba view -h -t 4 -p \
		--filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
		-f bam -l 0 \
		DE17NGSLABD100901_trimmed.bam \
		 -o filtered.bam
		 
Run in Snakemake [workflow](Rules/filtering_sorting_reads/Snakefile).


Mark duplicates and Picard metrics
==================================
Mark duplicates: run different parameters for a random flow cell or a patterned flow cell. The optical pixel
distance should be changed accordingly, random = 100 and patterned = 2500 based on GATK best practices.

		picard -Xmx4g -Xms4g MarkDuplicates \
		-I DE07NGSLABD100887_filtered.sorted.bam \
		-O DE07NGSLABD100887.markdup.bam
		--REFERENCE_SEQUENCE /Path/to/Reference_Sequences/hg38.fa\
		--METRICS_FILE 	DE07NGSLABD100887_picard_markdup_raw_metrics.txt \
		--CREATE_INDEX false \
		--MAX_RECORDS_IN_RAM 1000 \
	        --SORTING_COLLECTION_SIZE_RATIO 0.15 \
		--ASSUME_SORTED true \
		--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

Convert BED file to interval list:

		picard -Xmx4g -Xms4g BedToIntervalList \
		-I /Path/to/Reference_Sequences/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed \
		-O covered_targets_Twist_Methylome_hg38_annotated_collapsed_intervals \
		-SD /Users/gyongyosilab/Documents/ccfDNA_Projekt/Reference_Sequences/hg38.dict

Generate performance metrics:
Coverage cap setting is increased to 1000 for accurate metrics at higher sequencing depths Picard HsMetrics calculate important performance metrics , such as Fold-80 Base Penalty, HS Library Size, Percent Duplicates, and Percent Off Bait.

		picard -Xmx4g -Xms4g CollectHsMetrics \
		--INPUT DE94NGSLABD100873.markdup.bam \
		--OUTPUT DE94NGSLABD100873_hs_metrics.txt \
		-R /Path/to/Reference_Sequences/hg38.fa \
		--BAIT_INTERVALS  /Path/to/Reference_Sequences/covered_targets_Twist_Methylome_hg38_annotated_collapsed_intervals \
    		--TARGET_INTERVALS /Path/to/Reference_Sequences/covered_targets_Twist_Methylome_hg38_annotated_collapsed_intervals \
		--MINIMUM_MAPPING_QUALITY 20 \
		--COVERAGE_CAP 1000 \
		--PER_TARGET_COVERAGE DE94NGSLABD100873_hsmetrics_pertargetcoverage.txt \
		--NEAR_DISTANCE 500

Generate additional performance metrics:

		picard -Xmx4g -Xms4g CollectMultipleMetrics \ 
  		--INPUT markdupl/DE17NGSLABD100901.markdup.bam \
    		--OUTPUT DE17NGSLABD100901multiplemetrics \
		-R /Path/to/Reference_Sequences/hg38.fa \
  		--PROGRAM CollectGcBiasMetrics \
    		--PROGRAM CollectInsertSizeMetrics \
      		--PROGRAM CollectAlignmentSummaryMetrics
  		
	


Apply Snakemake [workflow](Rules/picard/picard_metrics) rules.

