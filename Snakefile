import re
import os

configfile: "config.yaml"


DIR_FASTQ = os.path.join(config["path"]["dir_fastq"],"")
DIR_OUT = os.path.join(config["path"]["dir_out"],"")
DIR_QC = os.path.join(config["path"]["dir_qc"],"")
DIR_BAM = os.path.join(config["path"]["dir_bam"],"")
DIR_STATS = os.path.join(config["path"]["dir_stats"],"")
DIR_LOG = os.path.join(config["path"]["dir_log"],"")

#(wholenames,) = glob_wildcards("../fastq/{wholename}.fastq.gz")
(wholenames,) = glob_wildcards(DIR_FASTQ+"{wholename}.fastq.gz")
profiletypes = config["all"]["profiletypes"]
BINSIZES=config["QDNAseq"]["BINSIZES"]
imagetype=config["ACE"]["imagetype"]
ACEBINSIZES=config["ACE"]["ACEBINSIZES"]

def getnames():
    SAMPLES=dict()
    for wholename in wholenames:
        sample = re.match('[a-zA-Z0-9\-]*', wholename).group(0)
        #fastqfile="../fastq/"+wholename+".fastq.gz"
        fastqfile=DIR_FASTQ+wholename+".fastq.gz"
        SAMPLES[sample]=fastqfile
    return(SAMPLES)
SAMPLES=getnames()

setting = config["all"]["setting"]
SettingService = [
        expand(DIR_OUT + "{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed",binSize=BINSIZES),
        expand(DIR_OUT + "{binSize}kbp/summary.html", binSize=BINSIZES),
        expand(DIR_OUT + "{binSize}kbp/BED/{sample}_annotate_focalCNA.bed", binSize=BINSIZES ,sample=SAMPLES.keys()),
        #expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/{sample}/summary_{sample}.{imagetype}", imagetype=imagetype ,binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()),
        expand(DIR_OUT + "{binSize}kbp/ACE/{ploidy}N/segmentfiles/{sample}_segments.tsv", binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys()),
        expand(DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds", ACEbinSize=ACEBINSIZES),
        expand(DIR_OUT + "{ACEbinSize}kbp/profiles/call_cellularity_based/index.html",ACEbinSize=ACEBINSIZES),
        expand(DIR_OUT + DIR_QC + "qc-fastq/{sample}_fastqc.html", sample=wholenames),
        expand(DIR_OUT + DIR_QC + "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
        ]
SettingResearch = [
        SettingService
        ]

rule all:
    input:
        SettingService if setting=="service" else SettingResearch 

rule bwa_aln:
    input:
        lambda wildcards: SAMPLES[wildcards.sample],
    output:
        sai=temp(DIR_OUT + DIR_BAM +  "{sample}.sai"),
        samse=temp(DIR_OUT + DIR_BAM + "{sample}.samse.sam")
    params:
        ref= config['all']['REF'],
        n=config['bwa']['max_edit_distance'],
        q=config['bwa']['read_trimming_param'],
    threads: config['all']['THREADS']
    log: DIR_OUT + DIR_LOG + "bwa/{sample}.log"
    shell:
        "bwa aln -n {params.n} -t {threads} -q {params.q} {params.ref} {input} > {output.sai} 2> {log}; "
        "bwa samse -f {output.samse} -r '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}'" 
        " {params.ref} {output.sai} {input} 2>> {log}"

rule samtools_sort:
    input:
        samse=DIR_OUT + DIR_BAM + "{sample}.samse.sam",
	    sai=DIR_OUT + DIR_BAM + "{sample}.sai"
    output:
        temp(DIR_OUT + "tmp/{sample}.all.bam")
    params:
        output=DIR_OUT + "tmp/{sample}.all"
    log: DIR_OUT + DIR_LOG + "samtools/{sample}.log"
    shell:
        "samtools view -uS {input.samse} 2> {log}| samtools sort - -o {output} 2>> {log}"

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
        script="scripts/stats.sh"
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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/binReadCounts.log"
    script:
        "scripts/binReadCounts.R"

rule QDNAseq_normalize:
    input:
        binReadCounts=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-raw.rds",
    output:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/corrected/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/corrected/",
        chrom_filter=config["QDNAseq"]["chrom_filter"],
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/normalizeBins.log"
    script:
        "scripts/QDNAseq_normalize.R"

rule deWave:
    input:
        corrected=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-corrected.rds",
    output:
        dewaved=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-dewaved.rds",
        allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/dewaved/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles=DIR_OUT + "{binSize}kbp/profiles/dewaved/",
        dewave_data=config["QDNAseq"]["dewave_data"],
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/dewave.log"
    script:
        "scripts/deWave.R"

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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "QDNAseq/{binSize}kbp/segment.log"
    script:
        "scripts/QDNAseq_segment.R"

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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{binSize}kbp/call.log"
    script:
        "scripts/CNA_call.R"

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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{ACEbinSize}kbp/call_cellularity_based.log"
    script:
        "scripts/CNA_call_cellularity_based.R"

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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "CNA/{ACEbinSize}kbp/bedfiles.log"
    script:
        "scripts/makeCNAbedFile.R"

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
        "{input.script} {input.bedfile} {wildcards.sample} {params.outdir} {output} 2> {log} "

rule CGHregions:
#TODO change reCalled to normal calls or cellularity based
    input:
        #recalled=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-reCalled.rds",
        recalled=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-call_cellularity_based.rds",
    output:
        RegionsCGH=DIR_OUT + "{ACEbinSize}kbp/data/{ACEbinSize}kbp-RegionsCGH.rds",
        profiles=DIR_OUT + "{ACEbinSize}kbp/profiles/freqPlot/freqPlotREGIONS_{ACEbinSize}kbp.png"
    params:
        averr=config["CGHregions"]["averror"],
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "{ACEbinSize}kbp/CGHregions.log"
    script:
        "scripts/CGHregions.R"

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
        suppressMessages=config["all"]["suppressMessages"]
    log: DIR_OUT + DIR_LOG + "{binSize}kbp/CGHregionstable.log"
    script:
        "scripts/makeCGHregionstable.R"

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
        script="scripts/lane-summary.sh",
        qcfastq=expand(DIR_OUT + DIR_QC +  "qc-fastq/{wholename}_fastqc.html", wholename=wholenames),
        bamqc=expand(DIR_OUT + DIR_QC +  "qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
    output:
        DIR_OUT + "{binSize}kbp/summary.html"
    params:
        bamfolder=DIR_OUT + DIR_BAM + "",
    shell:
        "{input.script} {wildcards.binSize}kpb {params.bamfolder} > {output}"

rule qcfastq:
    input:
        fastq=DIR_FASTQ + "{sample}.fastq.gz"
        #fastq=expand(DIR_FASTQ + "fastq/{{sample}}.fastq.gz", sample=SAMPLES.keys())
    output:
        qcfastq=DIR_OUT + DIR_QC +  "qc-fastq/{sample}_fastqc.html",
        qczip=temp(DIR_OUT + DIR_QC +  "qc-fastq/{sample}_fastqc.zip")
        #qcfastq=expand(DIR_OUT + DIR_QC +  "qc-fastq/{{sample}}.fastqc.html", sample=SAMPLES.keys()),
        #qczip=temp(expand(DIR_OUT + DIR_QC +  "qc-fastq/{{sample}}_fastqc.zip", sample=SAMPLES.keys()))
    threads: config['all']['THREADS']
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
    threads: config['all']['THREADS']
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
        suppressMessages=config["all"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "ACE/{ACEbinSize}kbp/{ploidy}N/log.tsv"
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
        suppressMessages=config["all"]["suppressMessages"]
    log:DIR_OUT + DIR_LOG + "{ACEbinSize}kbp/ACE/{ploidy}N_log.tsv"
    script:
        "scripts/Run_postanalysisloop_ACE.R"

# #TODO change reCalled to normal calls or cellularity based
# rule CNA_recall:
#     input:
#         called=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-called.rds",
#     output:
#         recalled=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-reCalled.rds",
#         #copynumbers=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-copynumbers.igv",
#         #segments=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-segments.igv",
#         calls=DIR_OUT + "{binSize}kbp/data/{binSize}kbp-calls.igv",
#         allprofiles=expand(DIR_OUT + "{{binSize}}kbp/profiles/reCalled/{samples}.png",samples=SAMPLES.keys()),
#         #copynumbersbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-copynumbers.bed",samples=SAMPLES.keys()),
#         #segmentsbed=expand(DIR_OUT + "{{binSize}}kbp/BED/{samples}-segments.bed",samples=SAMPLES.keys()),
#     params:
#         profiles=DIR_OUT + "{binSize}kbp/profiles/reCalled/",
#         #copynumbersbed=DIR_OUT + "{binSize}kbp/BED/%s-copynumbers.bed",
#         #segmentsbed=DIR_OUT + "{binSize}kbp/BED/%s-segments.bed",
#         #bedfolder=DIR_OUT + "{binSize}kbp/BED/",
#         failed=DIR_OUT + "{binSize}kbp/logs/failed_samples.txt",
#     log: DIR_OUT + "{binSize}kbp/logs/recall.log"
#     script:
#         "scripts/CNA_recall.R"
