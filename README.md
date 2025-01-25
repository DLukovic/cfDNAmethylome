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

Calling methylation (MethylDackel)
==================================
Check for bias and clean reads before calling CpGs. In an ideal experiment, we expect that the probability of observing a methylated C is constant across the length
of any given read. In practice, however, there are often increases/decreases in observed methylation rate at the
ends of reads and/or more global changes. These are termed methylation bias and including such regions in
the extracted methylation metrics will result in noisier and less accurate data.

Example: <img width="294" alt="image" src="https://github.com/user-attachments/assets/b5a4b9cf-a82e-4b63-a69e-995aa6e8d8ea" />

Run mbias:

		 MethylDackel mbias /Path/to//Reference_Sequences/hg38.fa \
                 DE07NGSLABD100887_markdup.bam DE07NGSLABD100887_mbias

Use Snakemake [rule](Rules/mbias/Snakefile).


Results
-------
Suggested inclusion options:

		Cardio		--OT 5,0,0,80 --OB 0,96,22,0
		CardioTox   	--OT 5,0,0,80 --OB 0,96,20,0
		Healthy		--OT 5,0,0,80 --OB 0,96,20,0	
		NonCardioTox	--OT 5,0,0,80 --OB 0,96,20,0
		Pulmo		--OT 5,0,0,80 --OB 0,96,20,0


Explanation:

Explanation: 

		--OT
Include reads for the Original Top (OT) strand starting from position 5 and ending at position 80.
Exclude positions 1–4 and anything beyond position 80.

		--OB
Include reads for the Original Bottom (OB) strand starting from position 2 (position 1 is excluded).
End inclusion at position 97.

Run MethylDackel:

		MethylDackel extract --minDepth 5 \
			             --maxVariantFrac 0.25  \
				     --OT 5,101,1,80 \
				     --OB 1,96,20,101 \
		 		     --mergeContext \
				     /Path/to/Reference_Sequences/hg38.fa \
				     DE07NGSLABD100887.markdup.bam \
				     -o DE07NGSLABD100887_methylome_biased_report.txt


Run MethylDackel with Snakemake [workflow](Rules/MethylDackel/Snakefile). Output is .bedGraphfile.

awk '{methylated+=$4; unmethylated+=$5} END {print "Methylated:", methylated, "Unmethylated:", unmethylated, "Methylation Level:", methylated/(methylated+unmethylated)}' DE07NGSLABD100880_methylome_biased_CpG.bedGraph

Exclude XY chromosomes from the MethylDackel output to eliminate sex-related biases:

		awk '$1 != "X" && $1 != "Y"' "$input_file" > "$output_file"							 
		awk '$1 != "X" && $1 != "Y"'  DE02NGSLABD100880_CpG.bedGraph > DE02NGSLABD100880_CpG_noXY.bedGraph"	

Run Snakemake [workflow](Rules/bedGraph_preprocess/Snakefile) to preprocess .bedGraph files and retrieve statistics.

Calling DMRs using Metilene
===========================

1.    .bedGraph files preprocessing:
	.bedGraphfiles containg track and 6 columns- fix it: (Use awk to remove the header and keep only the 					first 4 columns:)

		awk 'NR > 1 {print $1, $2, $3, $4}' OFS='\t' Sample_1.bedGraph > Sample_1_clean.bedGraph
		
2.    create Matrix of samples/.bedGraph files

Bedtools unionbedg combines multiple .bedGraph files into a single file such that one can directly compare coverage (and other text-values such as genotypes) across multiple sample. .bedGraph files have to be tab delimited, sorted, containing 4 columns (chromosomes, start, end, methylation proportion).

		
		
		bedtools unionbedg -i Sample_1_clean_sorted.bedGraph Sample_2_clean_sorted.bedGraph \
				   -header \ 
				   -names cardio_1 pulmo_2 > input_metilene.bed

Use Snakemake [rule](Rules/matrix_prep/Snakefile) to construct matrix of .bedGraph output from MethylDackel.
			
3.	To prepare input for metilene, it is necessary to remove the 3th column ("end) from the matrix.		
		
		awk -F'\t' 'BEGIN {OFS="\t"} { $3=""; $1=$1; print }' merged_matrix.bed > matrix.bed	
		

Metilene :
The first group is considered as a control group (g1) and the second (g2) group is relative to control group
		
		metilene -a Healthy -b Pulmo matrix.bed | sort -V -k1,1 -k2,2n > metilene_output.bed 
		metilene -a Healthy -b Pulmo_LC matrix.bed | sort -V -k1,1 -k2,2n > metilene_output_healthy_pulmo_lc.bed

To visualise data, use Perls script:

		perl /System/Volumes/Data/Users/gyongyosilab/metilene_v0.2-8/metilene_output.pl
		perl /Users/gyongyosilab/metilene_v0.2-8/metilene_output.pl -q metilene_output.bed -a Healthy -b Cardio
		perl /Users/gyongyosilab/metilene_v0.2-8/metilene_output.pl -q metilene_output.bed -a Healthy -b Pulmo

Metilene output:

*metilene_output.bed* the raw, unfiltered output, all identified DMRs, also outside of gived threshold (0.05)

*metilene_qval.0.05.out* filtered according to the parameters set for the pipeline run.

Metilene_output.bed (RAW output) :

		chr	start	stop	q-value	     meth.diff	CpGs	p (MWU)		o(2D KS)	mean g1	mean g2	
		chr1	10469	10576	1		-5.749990	10		0.24023		0.38309		71.25	77
		chr1	135159	135280	1		2.923485	10		5.5992e-05	0.0032001	91.642	88.718


Metilene_qval.0.05.out:

		chr	start	stop	q-value		met.diff	CpGs	mean g1 mean g1	
		chr1	1439669	1439938	0.0011631	4.702479	44		82.068	77.366
		chr1	1629043	1630187	0.00060914	3.722891	88		19.764	16.041

q-val (adjusted p-Vlaue; per default Bonferroni adjusted based on MWU-test p-values )

How many DMRs were identified ?

		wc -l metilene_qval.0.05.out

Healthy vs Cardio_LC : 
* 89754 metilene_output.bed
* 292 metilene_qval.0.05.out







