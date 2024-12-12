import re
import os
import glob

configfile: "config.yaml"


DIR_FASTQ = os.path.join(config["path"]["dir_fastq"],"")
DIR_OUT = os.path.join(config["path"]["dir_out"],"")
DIR_QC = os.path.join(config["path"]["dir_qc"],"")
DIR_BAM = os.path.join(config["path"]["dir_bam"],"")
DIR_STATS = os.path.join(config["path"]["dir_stats"],"")
DIR_LOG = os.path.join(config["path"]["dir_log"],"")

(wholenames,) = glob_wildcards(DIR_FASTQ+"{wholename}_R1.fastq.gz")
profiletypes = config["summary"]["profiletypes"]
BINSIZES=config["QDNAseq"]["BINSIZES"]
imagetype=config["ACE"]["imagetype"]
ACEBINSIZES=config["ACE"]["ACEBINSIZES"]
setting = config["pipeline"]["setting"]

def getnames():
    SAMPLES=dict()
    for wholename in wholenames:
        sample = wholename
        #fastqfile=DIR_FASTQ+wholename+"*.fastq.gz"
        fastqfile=glob.glob(DIR_FASTQ+wholename+"*.fastq.gz")
        SAMPLES[sample]=fastqfile
    return(SAMPLES)

SAMPLES=getnames()
print(SAMPLES)

if setting == "service": #rule service
    rule service:
        input:
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
elif setting == "research": #rule research
    rule research:    
        input:
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
            #expand(DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed",binSize=BINSIZES), # CGH branch
            expand(DIR_OUT + "{binSize}kbp/BED/{sample}_annotate_focalCNA.bed", binSize=BINSIZES ,sample=SAMPLES.keys()), # annotate_focalCNA
            expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()), # postanalysisloop_ACE
            expand(DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds", ACEbinSize=ACEBINSIZES), # CNA_call_cellularity_based
            expand(DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/index.html",ACEbinSize=ACEBINSIZES), # lightBox
            #expand(DIR_OUT + DIR_QC + "qc-fastq/{sample}_fastqc.html", sample=wholenames), # qcfastq
            expand(DIR_OUT + DIR_QC + "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()), # qcbam
else: #rule all
    rule all:
        input: 
            expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES), # summary
            expand(DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed",binSize=BINSIZES), # CGHregions branch
            expand(DIR_OUT + "{ACEbinSize}kbp/CGHtest/combined.png", ACEbinSize=ACEBINSIZES), # CGHtest branch
            expand(DIR_OUT + "{binSize}kbp/BED/{sample}_annotate_focalCNA.bed", binSize=BINSIZES ,sample=SAMPLES.keys()), # annotate_focalCNA
            expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()), # postanalysisloop_ACE
            expand(DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds", ACEbinSize=ACEBINSIZES), # CNA_call_cellularity_based
            expand(DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/index.html",ACEbinSize=ACEBINSIZES), # lightBox
            #expand(DIR_OUT + DIR_QC + "qc-fastq/{sample}_fastqc.html", sample=wholenames), # qcfastq
            expand(DIR_OUT + DIR_QC + "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()), # qcbam
            ##expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/{sample}/summary_{sample}.{imagetype}", imagetype=imagetype ,binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()),

rule bwa_aln:
    input:
        lambda wildcards: SAMPLES[wildcards.sample],
    output:
        bam=temp(DIR_OUT + DIR_BAM + "{sample}.align.bam")
    params:
        ref=config['bwa']['REF'],
        n=config['bwa']['max_edit_distance'],
        q=config['bwa']['read_trimming_param'],
    threads: config['pipeline']['THREADS']
    log: DIR_OUT + DIR_LOG + "bwa/{sample}.log"
    shell:
        "bwa mem -t {threads} {params.ref} {input} | samtools view -bS - > {output.bam} 2> {log}"

rule samtools_sort:
    input:
        bam=DIR_OUT + DIR_BAM + "{sample}.align.bam"
    output:
        temp(DIR_OUT + "tmp/{sample}.all.bam")
    params:
        output=DIR_OUT + "tmp/{sample}.all"
    log: DIR_OUT + DIR_LOG + "samtools/{sample}.log"
    shell:
        "samtools sort {input.bam} -o {output} 2>> {log}"

rule mark_duplicates:
    input:
        DIR_OUT + "tmp/{sample}.all.bam"
    output:
        bam=DIR_OUT + DIR_BAM + "{sample}.bam",
        bai=DIR_OUT + DIR_BAM + "{sample}.bai",
        bambai=DIR_OUT + DIR_BAM + "{sample}.bam.bai",
        metrics_file=DIR_OUT + DIR_STATS + "{sample}.md.metrics"
    log: DIR_OUT + DIR_LOG + "mark_duplicates/{sample}.log"
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics_file}"
        " ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true &> {log}; "
        "ln -s {output.bai} {output.bambai}"

rule generate_stats:
    input:
        bam=DIR_OUT + DIR_BAM + "{sample}.bam",
        script="scripts/Run_generate_stats.sh"
    output:
        DIR_OUT + DIR_STATS + "{sample}.reads.all"
    params:
        outdir=DIR_OUT + DIR_STATS
    shell:
        "{input.script} {input.bam} {wildcards.sample} {params.outdir} {output}"

rule QDNAseq_binReadCounts:
    input:
        bams=expand(DIR_OUT + DIR_BAM + "{sample}.bam", sample=SAMPLES.keys()),
    output:
        binReadCounts=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-raw.rds"
    params:
        genome=config["QDNAseq"]["genome"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/binReadCounts.log"
    script:
        "scripts/Run_QDNAseq_binReadCounts.R"

rule QDNAseq_normalize:
    input:
        binReadCounts=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-raw.rds",
    output:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/corrected/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/corrected/",
        chrom_filter=config["QDNAseq"]["chrom_filter"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/normalizeBins.log"
    script:
        "scripts/Run_QDNAseq_normalize.R"

rule deWave:
    input:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
    output:
        dewaved=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-dewaved.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/dewaved/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/dewaved/",
        dewave_data=config["QDNAseq"]["dewave_data"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/dewave.log"
    script:
        "scripts/Run_deWave.R"

rule QDNAseq_segment:
    input:
        dewaved=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-dewaved.rds" if config['QDNAseq']['dewave'] else DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
    output:
        segmented=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segmented.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/segmented/{samples}.png",samples=SAMPLES.keys()),
        copynumbers=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-copynumbers.igv",
        segments=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segments.igv",
        copynumbersbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-copynumbers.bed",samples=SAMPLES.keys()),
        segmentsbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-segments.bed",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/segmented/",
        failed=DIR_OUT + "{binSize}kbp/failed_samples.txt",
        minimal_used_reads=config["QDNAseq"]["minimal_used_reads"],
        copynumbersbed=DIR_OUT + "{binSize}kbp/BED/%s-copynumbers.bed",
        segmentsbed=DIR_OUT + "{binSize}kbp/BED/%s-segments.bed",
        bedfolder=DIR_OUT + "{binSize}kbp/BED/",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/segment.log"
    script:
        "scripts/Run_QDNAseq_segment.R"

rule CNA_call:
    input:
        segmented=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segmented.rds",
    output:
        called=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-called.rds",
        freqplot=DIR_OUT + "{binSize}kbp/profiles/freqPlot/freqPlot_{binSize}kbp.png",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/called/{samples}.png",samples=SAMPLES.keys()),
        calls=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-calls.igv"
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/called/",
        failed=DIR_OUT + "{binSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{binSize}kbp/call.log"
    script:
        "scripts/Run_CNA_call.R"

rule CNA_call_cellularity_based:
    input:
        called=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-called.rds",
        fitpicker=expand(DIR_OUT + "{{ACEbinSize}}kbp/ACE/{main_ploidy}N/fitpicker_{main_ploidy}N.tsv", main_ploidy=config["ACE"]["main_ploidy"]),
    output:
        recalled=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds",
        calls=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-calls_cellularity_based.igv",
        allprofiles=expand(DIR_OUT + "{{ACEbinSize}}kbp/profiles/call_cellularity_based/{samples}.png",samples=SAMPLES.keys()),
        freqplot=DIR_OUT + "{ACEbinSize}kbp/profiles/freqPlot/freqPlot_{ACEbinSize}kbp_cellularity.png"
    params:
        profiles=DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/",
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        minimum_cellularity=config["QDNAseq"]["minimum_cellularity"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{ACEbinSize}kbp/call_cellularity_based.log"
    script:
        "scripts/Run_CNA_call_cellularity_based.R"

rule CNA_bedfiles:
    input:
        #recalled=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-reCalled.rds",
        recalled=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds",
    output:
        bedfile=expand(DIR_OUT + "{{ACEbinSize}}kbp/BED/{sample}_allCNAsPerBin.bed", sample=SAMPLES.keys()),
        focalCNA=expand(DIR_OUT + "{{ACEbinSize}}kbp/BED/{sample}_focalCNAs.bed", sample=SAMPLES.keys()),
        CNAs=expand(DIR_OUT + "{{ACEbinSize}}kbp/BED/{sample}_CNAs.bed", sample=SAMPLES.keys()),
    params:
        beddir=DIR_OUT + "{ACEbinSize}kbp/BED/",
        cytobands=config["CGHregions"]["cytobands"],
        max_focal_size_bed=config["BED"]["max_focal_size_bed"],
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{ACEbinSize}kbp/bedfiles.log"
    script:
        "scripts/Run_makeCNAbedFile.R"

rule annotate_focalCNA:
 #TODO remove in annotateFocalCNAbed file hardcoded location: /net/nfs/PAT/home/stef/code/ENSEMBL_API/ensembl74/ensembl/modules/
 #TODO remove from addCosmicCencus census location: /net/nfs/PAT/home/matias/data/ref/cosmic/hg19_okt2015/CosmicMutantExport.tsv
    input:
        bedfile=DIR_OUT + "{binSize}kbp/BED/{sample}_focalCNAs.bed",
        script="scripts/annotateFocalCNAbed.sh"
    output:
        DIR_OUT + "{binSize}kbp/BED/{sample}_annotate_focalCNA.bed"
    params:
        outdir=DIR_OUT + "{binSize}kbp/BED/"
    log: DIR_OUT + DIR_LOG + "CNA/{binSize}kbp/{sample}_annotate_focalCNA.log"
    shell:
        "{input.script} {input.bedfile} {wildcards.sample} {params.outdir} {output} {log} 2> {log}"

rule CGHregions:
    input:
        recalled=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds",
    output:
        RegionsCGH=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-RegionsCGH.rds",
        profiles=DIR_OUT + "{ACEbinSize}kbp/profiles/freqPlot/freqPlotREGIONS_{ACEbinSize}kbp.png"
    params:
        averr=config["CGHregions"]["averror"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "{ACEbinSize}kbp/CGHregions.log"
    script:
        "scripts/Run_CGHregions.R"

rule makeCGHregionsTable:
    input:
        RegionsCGH=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-RegionsCGH.rds",
    output:
        allRegions=DIR_OUT + "{binSize}kbp/profiles/freqPlot/allRegions.txt",
        allFocalRegions=DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.bed"
    params:
        min_freq_focal=config["CGHregions"]["min_freq_focal"],
        max_focal_size_mb=config["CGHregions"]["max_focal_size_mb"],
        cytobands=config["CGHregions"]["cytobands"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "{binSize}kbp/CGHregionstable.log"
    script:
        "scripts/Run_makeCGHregionstable.R"

rule annotate_RegionsFocalCNAbed:
    input:
        bedfile=DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.bed",
        script="scripts/annotateFocalCNAbed.sh"
    output:
        output=DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed"
    params:
        outdir=DIR_OUT + "{binSize}kbp/profiles/freqPlot/",
        filename="allFocalRegions"
    log: DIR_OUT + DIR_LOG + "CNA/{binSize}kbp/annotate_RegionsFocalCNAbed.log"
    shell:
        "{input.script} {input.bedfile} {params.filename} {params.outdir} {output} 2> {log}"

rule lightBox:
    input:
        sample=expand(DIR_OUT + "{{binSize}}kbp/profiles/{{profiletype}}/{sample}.png", sample=SAMPLES.keys()),
        script="scripts/createLightBox.sh",
    output:
        index=DIR_OUT + "{binSize}kbp/profiles/{profiletype}/index.html",
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/{profiletype}/",
        lb2dir="lb2/"
    shell:
        "{input.script} {params.profiles} {params.lb2dir} > {output.index}"

rule summary:
#TODO script lane-summary does contain relative links to files - maybe better to change
    input:
        stats=expand(DIR_OUT + DIR_STATS + "{sample}.reads.all", sample=SAMPLES.keys()),
        index=expand(DIR_OUT + "{{binSize}}kbp/profiles/{profiletype}/index.html", profiletype=profiletypes),
        script="scripts/Run_lane-summary.sh",
        qcfastq=expand(DIR_OUT + DIR_QC +  "qc-bam/{wholename}_fastqc.html", wholename=wholenames),
        bamqc=expand(DIR_OUT + DIR_QC +  "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
    output:
        DIR_OUT + "{binSize}kbp/summary.html"
    params:
        bamfolder=DIR_OUT + DIR_BAM,
        statsfolder = DIR_OUT + DIR_STATS,
        qcFastqfolder= DIR_OUT + DIR_QC + "qc-fastq/",
        qcBamfolder= DIR_OUT + DIR_QC + "qc-bam/"
    shell:
        "{input.script} {wildcards.binSize}kpb {params.bamfolder} {params.statsfolder} {params.qcFastqfolder} {params.qcBamfolder}> {output}"

rule qcfastq:
    input:
        fastq=DIR_FASTQ + "{sample}.fastq.gz"
        #fastq=expand(DIR_FASTQ + "fastq/{{sample}}.fastq.gz", sample=SAMPLES.keys())
    output:
        qcfastq=DIR_OUT + DIR_QC +  "qc-fastq/{sample}_fastqc.html",
        qczip=temp(DIR_OUT + DIR_QC +  "qc-fastq/{sample}_fastqc.zip")
        #qcfastq=expand(DIR_OUT + DIR_QC +  "qc-fastq/{{sample}}.fastqc.html", sample=SAMPLES.keys()),
        #qczip=temp(expand(DIR_OUT + DIR_QC +  "qc-fastq/{{sample}}_fastqc.zip", sample=SAMPLES.keys()))
    threads: config['pipeline']['THREADS']
    params:
        qcfastqdir=DIR_OUT + DIR_QC + "qc-fastq/"
    log: DIR_OUT + DIR_LOG + DIR_QC + "qc-fastq/{sample}.log"
    shell:
        "fastqc {input.fastq} --outdir {params.qcfastqdir} -t {threads} 2>> {log}"

rule qcbam:
    input:
        bam=DIR_OUT + DIR_BAM + "{sample}.bam",
    output:
        bamqc=DIR_OUT + DIR_QC +  "qc-bam/{sample}_fastqc.html",
        bamqczip=temp(DIR_OUT + DIR_QC +  "qc-bam/{sample}_fastqc.zip"),
    threads: config['pipeline']['THREADS']
    params:
        qcbamdir=DIR_OUT + DIR_QC + "qc-bam/"
    log: DIR_OUT + DIR_LOG + DIR_QC + "qc-bam/{sample}.log"
    shell:
        "fastqc {input.bam} --format bam_mapped --outdir {params.qcbamdir} -t {threads} 2>> {log}"

rule ACE:
    input:
        segmented=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-segmented.rds",
    output:
        #ACE=expand(DIR_OUT + "{{ACEbinSize}}kbp/ACE/{{ploidy}}N/summary_files/summary_{sample}.{{imagetype}}", sample=SAMPLES.keys()),
        fitpicker=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/fitpicker_{ploidy}N.tsv",
    params:
        outputdir=DIR_OUT + "{ACEbinSize}kbp/ACE/",
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "ACE/{ACEbinSize}kbp/{ploidy}N/ACE_log.tsv"
    script:
        "scripts/Run_ACE.R"

rule postanalysisloop_ACE:
    input:
        segmented=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-segmented.rds",
        fitpicker=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/fitpicker_{ploidy}N.tsv",
    output:
        ACE_post=expand(DIR_OUT + "{{ACEbinSize}}kbp/ACE/{{ploidy}}N/segmentfiles/{sample}_segments.tsv", sample=SAMPLES.keys()),
    params:
        outputdir=DIR_OUT + "{ACEbinSize}kbp/ACE/{ploidy}N/",
        failed=DIR_OUT + "{ACEbinSize}kbp/failed_samples.txt",
        suppressMessages=config["pipeline"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "ACE/{ACEbinSize}kbp/{ploidy}N/ACE_post_log.tsv"
    script:
        "scripts/Run_postanalysisloop_ACE.R"

rule CGHtest:
    input:
        RegionsCGH=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-RegionsCGH.rds",
    output:
        freqPlotCompare=DIR_OUT + "{ACEbinSize}kbp/CGHtest/freqPlotCompare.png", 
        CGHtest=DIR_OUT + "{ACEbinSize}kbp/CGHtest/CGHtest.Rdata",
        plotPFDR=DIR_OUT + "{ACEbinSize}kbp/CGHtest/plotPFDR.png",
        combined=DIR_OUT + "{ACEbinSize}kbp/CGHtest/combined.png"
    params:
        outputdir=DIR_OUT + "{ACEbinSize}kbp/CGHtest/",
        clinicaldataPath=config["CGHtest"]["clinicaldataPath"],
        columnSampleNames=config["CGHtest"]["columnSampleNames"],
        ClassSamples=config["CGHtest"]["ClassSamples"],
        columnClassSamples=config["CGHtest"]["columnClassSamples"],
        suppressMessages=config["pipeline"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "{ACEbinSize}kbp/CGHtest_log.tsv"
    script:
        "scripts/Run_CGHtest.R"

