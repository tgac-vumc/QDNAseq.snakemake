#!/bin/sh

lib="/net/nfs/PAT/lib/"
ref=${ref:-"/net/nfs/PAT/data/ref/hg19/ensembl/Homo_sapiens.GRCh37.69.dna_sm.primary_assembly.fa.gz"}
fastqc="$lib/FastQC/fastqc --threads=2"
bwa="$lib/bwa/bwa-0.5.9/bwa"
samtools="$lib/samtools/samtools-0.1.18/samtools"
picard="$lib/picard-tools/picard-tools-1.61"
threads=${bwaThreads:-"4"}

for xin in "$@"
do
  xout=`basename $xin`

  # ($fastqc ${xin}.fastq* -o ../qc-fastq; rm ../qc-fastq/*.zip) &

  $bwa aln -n 2 -t $threads -q 40 "$ref" ${xin}.f* > "${xout}.sai"
  $bwa samse -r "@RG\tID:${xout}\tSM:${xout}" -f "${xout}.samse.sam" "$ref" "${xout}.sai" ${xin}.f*
  rm "${xout}.sai"

# sort
  $samtools view -uS "${xout}.samse.sam" | $samtools sort - "${xout}.all"
  rm "${xout}.samse.sam"

#  $samtools view -c "${xout}.all.bam" > "../stats/${xout}.reads.all"

# mark duplicates
  java -Xmx2g \
    -jar $picard/MarkDuplicates.jar \
    I="${xout}.all.bam" O="${xout}.md.bam" M="../stats/${xout}.md.metrics" AS=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
  rm "${xout}.all.bam"

  $samtools view -F 0x0404 "${xout}.md.bam" | awk "BEGIN {FS=\"\t\";uni=0;q1=0;q37=0} \
	{if (\$2 ~ /[^d]/ && \$5 >= 1) q1 += 1; \
	if (\$2 ~ /[^d]/ && \$5 >= 37) q37 += 1 } \
	END {print NR > \"../stats/${xout}.reads.unique\";\
	print q1 > \"../stats/${xout}.reads.q1\";\
	print q37 > \"../stats/${xout}.reads.q37\"}"

  $samtools idxstats "${xout}.md.bam" | awk "BEGIN {FS=\"\t\";aligned=0;unaligned=0} {aligned += \$3; unaligned += \$4} END {print aligned > \"../stats/${xout}.reads.aligned\"; print aligned + unaligned > \"../stats/${xout}.reads.all\" } "

  mv "${xout}.md.bam" "${xout}.bam"
  mv "${xout}.md.bai" "${xout}.bai"
  ln -s "${xout}.bai" "${xout}.bam.bai"
  # ($fastqc "${xout}.bam" -o ../qc-bam; rm ../qc-bam/*.zip) &
done

# EOF
