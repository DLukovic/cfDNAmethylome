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
		
Downsampling process using seqtk
=========================================
+	Inspect the number of reads and pasepairs in fastq file

		gzip -dc DE02NGSLABD100880_1.fq.gz | 
    	 awk 'NR%4==2{c++; l+=length($0)}
          	END{
                print "Number of reads: "c; 
                print "Number of bases in reads: "l
              }'
	`gzip -dc` inspect fasq file without having  to decompress the file 
	
	`awk` operates on each line of input from the decompressed fastq file focusing on the certain lines
	
	`NR%4==2` NR is build in variable in awk that tracks the record number of line number, in this case awk operates only on the second line of every 4 line group
	
	`{c++; l+=length($0)}` or accumulating counts
	
	`c++`	increments the read count (c) by 1 each time sequence line is encountered, providing the total count of reads
	
	`l += length($0)` adds the length of the current sequence line ($0 represents the entire line) to l, accumulating the total number of bases across all reads.
	
	`END` block fter all lines are processed, the END block is executed. Here, the total number of reads (c) and the total number of bases (l) are printed.
print "Number of reads: " c outputs the count of reads.
print "Number of bases in reads: " l outputs the cumulative number of bases.

Downsample to final total reads according to formula:

Total wanted basepairs: count of base pairs in merged probes x 250 raw coverage

Final Total reads = total wanted basepairs / (2*100[forward and reverse reads in 100 bp])

+ Count total number of base pairs in all samples

		for file in /input_path/*.gz; do
   		gzip -dc "$file" | 
       	 awk 'NR%4==2{c++; l+=length($0)}
             	END{
                   print "Number of reads: "c; 
                   print "Number of bases in reads: "l
                 }'
		   /output_path
		done

Reference Sequence Directory 
===============================

		java -Xmx4g -Xms4g -jar picard.jar CreateSequenceDictionary REFERENCE=/Users/gyongyosilab/Documents/ccfDNA\ Projekt/Reference\ Sequences/hg38.fa  	OUTPUT=/Users/gyongyosilab/Documents/ccfDNA\ Projekt/Reference\ Sequences/hg38.dict


Alignment
==========
install bwameth.py in conda environment

		conda install bioconda::bwameth

create environment for alignment

		conda create --name methylenv

create indices
		
   		bwameth.py index /Users/gyongyosilab/Documents/ccfDNA_Projekt/Reference_Sequences/hg38.fa 

Important: add --read group information, these annotations help downstram tools distinguish sampels or batches of reads

Execute bwameethalignment using *snakemake* workflow management system.


