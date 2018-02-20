import re

configfile: "config.yaml"

(wholenames,) = glob_wildcards("../fastq/{wholename}.fastq.gz")
profiletypes = config["all"]["profiletypes"]
BINSIZES=config["QDNAseq"]["BINSIZES"]
imagetype=config["ACE"]["imagetype"]
ACEBINSIZES=config["ACE"]["ACEBINSIZES"]

def getnames():
    SAMPLES=dict()
    for wholename in wholenames:
        sample = re.match('[a-zA-Z0-9\-]*', wholename).group(0)
        fastqfile="../fastq/"+wholename+".fastq.gz"
        SAMPLES[sample]=fastqfile
    return(SAMPLES)

SAMPLES=getnames()

rule all:
    input:
        expand("../{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed",binSize=BINSIZES),
        expand("../{binSize}kbp/summary.html", binSize=BINSIZES),
        expand("../{binSize}kbp/BED/{sample}_Cosmic.bed", binSize=BINSIZES ,sample=SAMPLES.keys()),
        expand("../{binSize}kbp/ACE/{ploidy}N/{sample}/summary_{sample}.{imagetype}", imagetype=imagetype ,binSize=ACEBINSIZES, ploidy=config["ACE"]["ploidies"], sample=SAMPLES.keys())


rule bwa_aln:
    input:
        lambda wildcards: SAMPLES[wildcards.sample],
    output:
        sai=temp("../bam/{sample}.sai"),
        samse=temp("../bam/{sample}.samse.sam")
    params:
        ref= config['all']['REF'],
        n=config['bwa']['max_edit_distance'],
        q=config['bwa']['read_trimming_param'],
    threads: config['all']['THREADS']
    log: "../logs/bwa/{sample}.log"
    shell:
        "bwa aln -n {params.n} -t {threads} -q {params.q} {params.ref} {input} > {output.sai} 2> {log};"
        "bwa samse -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -f {output.samse}"
        " {params.ref} {output.sai} {input} 2>> {log}"

rule samtools_sort:
    input:
        samse="../bam/{sample}.samse.sam",
	    sai="../bam/{sample}.sai"
    output:
        temp("../tmp/{sample}.all.bam")
    params:
        output="../tmp/{sample}.all"
    log: "../logs/samtools/{sample}.log"
    shell:
        "samtools view -uS {input.samse} 2> {log}| samtools sort - -o {output} 2>> {log}"

rule mark_duplicates:
    input:
        "../tmp/{sample}.all.bam"
    output:
        bam="../bam/{sample}.bam",
        bai="../bam/{sample}.bai",
        bambai="../bam/{sample}.bam.bai",
        metrics_file="../stats/{sample}.md.metrics"
    log: "../logs/mark_duplicates/{sample}.log"
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics_file} "
        "ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT &> {log} ;"
        "ln -s {output.bai} {output.bambai}"

rule generate_stats:
    input:
        bam="../bam/{sample}.bam",
        script="scripts/stats.sh"
    output:
        "../stats/{sample}.reads.all"
    params:
        outdir="../stats"
    shell:
        "{input.script} {input.bam} {wildcards.sample} {params.outdir} {output}"

rule QDNAseq_binReadCounts:
    input:
        bams=expand("../bam/{sample}.bam", sample=SAMPLES.keys()),
        script="scripts/binReadCounts.R"
    output:
        binReadCounts="../{binSize}kbp/data/{binSize}kbp-raw.rds"
    params:
        genome=config["QDNAseq"]["genome"],
    log: "../{binSize}kbp/logs/binReadCounts.log"
    script:
        "{input.script}"

rule QDNAseq_normalize:
    input:
        binReadCounts="../{binSize}kbp/data/{binSize}kbp-raw.rds",
        script="scripts/QDNAseq_normalize.R",
    output:
        corrected="../{binSize}kbp/data/{binSize}kbp-corrected.rds",
        allprofiles=expand("../{{binSize}}kbp/profiles/corrected/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles="../{binSize}kbp/profiles/corrected/",
    log: "../{binSize}kbp/logs/normalizeBins.log"
    script:
        "{input.script}"

rule deWave:
    input:
        corrected="../{binSize}kbp/data/{binSize}kbp-corrected.rds",
        script="scripts/deWave.R",
    output:
        dewaved="../{binSize}kbp/data/{binSize}kbp-dewaved.rds",
        allprofiles=expand("../{{binSize}}kbp/profiles/dewaved/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles="../{binSize}kbp/profiles/dewaved/",
        dewave_data=config["QDNAseq"]["dewave_data"],
    log: "../{binSize}kbp/logs/dewave.log"
    script:
        "{input.script}"

rule QDNAseq_segment:
    input:
        dewaved="../{binSize}kbp/data/{binSize}kbp-dewaved.rds",
        script="scripts/QDNAseq_segment.R",
    output:
        segmented="../{binSize}kbp/data/{binSize}kbp-segmented.rds",
        allprofiles=expand("../{{binSize}}kbp/profiles/segmented/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles="../{binSize}kbp/profiles/segmented/",
        failed="../{binSize}kbp/logs/failed_samples.txt",
        minimal_used_reads=config["QDNAseq"]["minimal_used_reads"]
    log: "../{binSize}kbp/logs/segment.log"
    script:
        "{input.script}"

rule CNA_call:
    input:
        segmented="../{binSize}kbp/data/{binSize}kbp-segmented.rds",
        script="scripts/CNA_call.R",
    output:
        called="../{binSize}kbp/data/{binSize}kbp-called.rds",
        freqplot="../{binSize}kbp/profiles/freqPlot/freqPlot_{binSize}kbp.png",
        allprofiles=expand("../{{binSize}}kbp/profiles/called/{samples}.png",samples=SAMPLES.keys()),
    params:
        profiles="../{binSize}kbp/profiles/called/",
        failed="../{binSize}kbp/logs/failed_samples.txt",
    log: "../{binSize}kbp/logs/call.log"
    script:
        "{input.script}"

rule CNA_recall:
    input:
        called="../{binSize}kbp/data/{binSize}kbp-called.rds",
        script="scripts/CNA_recall.R",
    output:
        recalled="../{binSize}kbp/data/{binSize}kbp-reCalled.rds",
        copynumbers="../{binSize}kbp/data/{binSize}kbp-copynumbers.igv",
        segments="../{binSize}kbp/data/{binSize}kbp-segments.igv",
        calls="../{binSize}kbp/data/{binSize}kbp-calls.igv",
        allprofiles=expand("../{{binSize}}kbp/profiles/reCalled/{samples}.png",samples=SAMPLES.keys()),
        copynumbersbed=expand("../{{binSize}}kbp/BED/{samples}-copynumbers.bed",samples=SAMPLES.keys()),
        segmentsbed=expand("../{{binSize}}kbp/BED/{samples}-segments.bed",samples=SAMPLES.keys()),
    params:
        profiles="../{binSize}kbp/profiles/reCalled/",
        copynumbersbed="../{binSize}kbp/BED/%s-copynumbers.bed",
        segmentsbed="../{binSize}kbp/BED/%s-segments.bed",
        bedfolder="../{binSize}kbp/BED/",
        failed="../{binSize}kbp/logs/failed_samples.txt",
    log: "../{binSize}kbp/logs/recall.log"
    script:
        "{input.script}"

rule CNA_bedfiles:
    input:
        recalled="../{binSize}kbp/data/{binSize}kbp-reCalled.rds",
        script="scripts/makeCNAbedFile.R",
    output:
        bedfile=expand('../{{binSize}}kbp/BED/{sample}_allCNAsPerBin.bed', sample=SAMPLES.keys()),
        focalCNA=expand('../{{binSize}}kbp/BED/{sample}_focalCNAs.bed', sample=SAMPLES.keys()),
        CNAs=expand('../{{binSize}}kbp/BED/{sample}_CNAs.bed',sample=SAMPLES.keys()),
    params:
        beddir='../{binSize}kbp/BED/',
        cytobands=config["CGHregions"]["cytobands"],
        max_focal_size_bed=config["BED"]["max_focal_size_bed"],
        failed="../{binSize}kbp/logs/failed_samples.txt",
    log: "../{binSize}kbp/logs/bedfiles.log"
    script:
        "{input.script}"

#TODO remove in annotateFocalCNAbed file hardcoded location: /net/nfs/PAT/home/stef/code/ENSEMBL_API/ensembl74/ensembl/modules/
#TODO remove from addCosmicCencus census location: /net/nfs/PAT/home/matias/data/ref/cosmic/hg19_okt2015/CosmicMutantExport.tsv
rule annotate_focalCNA:
    input:
        bedfile='../{binSize}kbp/BED/{sample}_focalCNAs.bed',
        script="scripts/annotateFocalCNAbed.sh"
    output:
        "../{binSize}kbp/BED/{sample}_Cosmic.bed"
    params:
        outdir="../{binSize}kbp/BED/"
    log: "../{binSize}kbp/logs/annotate_focalCNA.log"
    shell:
        "{input.script} {input.bedfile} {wildcards.sample} {params.outdir} {output} 2> {log} "

rule CGHregions:
    input:
        recalled="../{binSize}kbp/data/{binSize}kbp-reCalled.rds",
        script="scripts/CGHregions.R"
    output:
        RegionsCGH="../{binSize}kbp/data/{binSize}kbp-RegionsCGH.rds",
        profiles='../{binSize}kbp/profiles/freqPlot/freqPlotREGIONS_{binSize}kbp.png'
    params:
        averr=config["CGHregions"]["averror"],
    log: "../{binSize}kbp/logs/CGHregions.log"
    script:
        "{input.script}"

rule makeCGHregionsTable:
    input:
        RegionsCGH="../{binSize}kbp/data/{binSize}kbp-RegionsCGH.rds",
        script="scripts/makeCGHregionstable.R",
    output:
        allRegions="../{binSize}kbp/profiles/freqPlot/allRegions.txt",
        allFocalRegions="../{binSize}kbp/profiles/freqPlot/allFocalRegions.bed"
    params:
        min_freq_focal=config["CGHregions"]["min_freq_focal"],
        max_focal_size_mb=config["CGHregions"]["max_focal_size_mb"],
        cytobands=config["CGHregions"]["cytobands"],
    log: "../{binSize}kbp/logs/CGHregionstable.log"
    script:
        "{input.script}"

rule annotate_RegionsFocalCNAbed:
    input:
        bedfile="../{binSize}kbp/profiles/freqPlot/allFocalRegions.bed",
        script="scripts/annotateFocalCNAbed.sh"
    output:
        output="../{binSize}kbp/profiles/freqPlot/allFocalRegions.Cosmic.bed"
    params:
        outdir="../{binSize}kbp/profiles/freqPlot/",
        filename="allFocalRegions"
    log: "../{binSize}kbp/logs/annotatateregions.log"
    shell:
        "{input.script} {input.bedfile} {params.filename} {params.outdir} {output} 2> {log}"

rule lightBox:
    input:
        sample=expand("../{{binSize}}kbp/profiles/{{profiletype}}/{sample}.png", sample=SAMPLES.keys()),
        script="scripts/createLightBox.sh",
    output:
        index="../{binSize}kbp/profiles/{profiletype}/index.html",
    params:
        profiles="../{binSize}kbp/profiles/{profiletype}/",
        lb2dir="lb2/"
    shell:
        "{input.script} {params.profiles} {params.lb2dir} > {output.index}"

#TODO script lane-summary does contain relative links to files - maybe better to change
rule summary:
    input:
        stats=expand("../stats/{sample}.reads.all", sample=SAMPLES.keys()),
        index=expand("../{{binSize}}kbp/profiles/{profiletype}/index.html", profiletype=profiletypes),
        script='scripts/lane-summary.sh',
        qcfastq=expand("../qc-fastq/{wholename}_fastqc.html", wholename=wholenames),
        bamqc=expand("../qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
    output:
        "../{binSize}kbp/summary.html"
    params:
        bamfolder="../bam/",
    shell:
        "{input.script} {wildcards.binSize}kpb {params.bamfolder} > {output}"

rule qcfastq:
    input:
        expand("../fastq/{sample}.fastq.gz", sample=wholenames)
    output:
        qcfastq=expand("../qc-fastq/{sample}_fastqc.html", sample=wholenames),
        qczip=temp(expand("../qc-fastq/{sample}_fastqc.zip", sample=wholenames))
    threads: config['all']['THREADS']
    params:
        qcfastqdir="../qc-fastq/"
    log: "../logs/fastqc.log"
    shell:
        "fastqc {input} --outdir {params.qcfastqdir} -t {threads} 2> {log}"

rule qcbam:
    input:
        bam=expand("../bam/{sample}.bam", sample=SAMPLES.keys()),
    output:
        bamqc=expand("../qc-bam/{sample}_fastqc.html", sample=SAMPLES.keys()),
        bamqczip=temp(expand("../qc-bam/{sample}_fastqc.zip", sample=SAMPLES.keys())),
    threads: config['all']['THREADS']
    params:
        qcbamdir="../qc-bam/"
    log: "../logs/fastqc-bam.log"
    shell:
        "fastqc {input.bam} --format bam_mapped --outdir {params.qcbamdir} -t {threads} 2> {log}"

rule ACE:
    input:
        segmented="../{ACEbinSize}kbp/data/{ACEbinSize}kbp-segmented.rds",
        script="scripts/Run_ACE.R"
    output:
        ACE=expand("../{{ACEbinSize}}kbp/ACE/{{ploidy}}N/{sample}/summary_{sample}.{{imagetype}}", sample=SAMPLES.keys())
    params:
        outputdir="../{ACEbinSize}kbp/ACE/",
        failed="../{ACEbinSize}kbp/logs/failed_samples.txt",
    log:"../{ACEbinSize}kbp/ACE/{ploidy}N/log.tsv"
    script:
        "{input.script}"
