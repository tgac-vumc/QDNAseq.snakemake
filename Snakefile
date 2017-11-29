import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
REF = "/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/version0.6.0/genome.fa"

SAMPLES, = glob_wildcards("fastq/{sample}.fastq.gz")
BINSIZES = "1 5 10 15 30 100 1000".split()  

rule all:
     input:
        expand('binAnnotations/{binSize}k.rds', binSize=BINSIZES)

rule bwa_mem:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "bam/{sample}.bam"
    threads: 12
    shell:
        "bwa mem -t 12 {REF} {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        "bam/{sample}.bam"
    output:
        "bam/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"

rule mark_duplicates:
    input:
        "bam/{sample}.sorted.bam"
    output:
        "bam/{sample}.sorted.md.bam"
    shell:
        "picard MarkDuplicates I={input} O={output} M=stats/{wildcards.sample}.sorted.md.metrics "
        "AS=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"

rule generate_stats:
    input:
        "bam/{sample}.sorted.md.bam"
    output:
        "stats/{sample}.reads.all"
    shell:
        "scripts/stats.sh {input} {wildcards.sample} stats {output}"

rule QDNAseq:
    input:
        bams=expand("bam/{sample}.sorted.md.bam", sample=SAMPLES)
    output:
        "rds/{sample}kbp-raw.rds"
    script:
        "scripts/QDNAseq.R"

rule QDNAseq_dewaved:
    input:
        #bams=expand("bam/{sample}.sorted.md.bam", sample=SAMPLES)
	"rds/{sample}kbp-raw.rds"
    output:
        "rds/{sample}kbp-dewaved.rds"
    script:
        "scripts/QDNAseq.R"

rule QDNAseq_CreateBinAnnotations:
    input:
        bams=expand("bam/{sample}.sorted.md.bam", sample=SAMPLES),
        mappability="rds/mappabilityBins-{binSize}k.rds"
    output:
        "binAnnotations/{binSize}k.rds"
    script:
        "scripts/CreateBinAnnotations.R"

rule QDNAseq_CreateMappability:
    input:
        bw="data/mappability/hg38.50mer.bw"
    output:
        "rds/mappabilityBins-{binSize}k.rds"
    params:
        binSize="{binSize}"
    log:
        "logs/QDNAseq_CreateMappability_{binSize}.log"
    script:
        "scripts/CreateMappability.R"

