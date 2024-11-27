import os
import snakemake.io
import glob

(SAMPLES,READS,) = glob_wildcards("reads/{sample}_{read}.fq.gz")
READS=["1","2"]
REF="/Users/gyongyosilab/Documents/ccfDNA_Projekt/Reference_Sequences/hg38.fa"

rule all:
    input: expand("reads/{sample}.bam",sample=SAMPLES)

rule bwameth:
    input:
        ref=REF,
        # determine `r1` based on the {sample} wildcard defined in `output`
        # and the fixed value `1` to indicate the read direction
        r1="reads/{sample}_1.fq.gz",
        # determine `r2` based on the {sample} wildcard similarly
        r2="reads/{sample}_2.fq.gz"

    output: 
        "reads/{sample}.bam"

    # better to pass in the threads than to hardcode them in the shell command
    threads: 8
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    conda: 
        "methylenv"
    shell: 
        "bwameth.py --reference {input.ref}  {input.r1} {input.r2} -t {threads} --read-group '{params.rg}' | samtools view -b -h - > {output}"

