#!/bin/bash

input=$1
output=$2
outdir=$3

samtools view -F 0x0404 "${input}" | awk "BEGIN {FS=\"\t\";uni=0;q1=0;q37=0} \
        {if (\$2 ~ /[^d]/ && \$5 >= 1) q1 += 1; \
        if (\$2 ~ /[^d]/ && \$5 >= 37) q37 += 1 } \
        END {print NR > \"${outdir}/${output}.reads.unique\";\
        print q1 > \"${outdir}/${output}.reads.q1\";\
        print q37 > \"${outdir}/${output}.reads.q37\"}"

samtools idxstats "${input}" | awk "BEGIN {FS=\"\t\";aligned=0;unaligned=0} {aligned += \$3; unaligned += \$4} END {print aligned > \"${outdir}/${output}.reads.aligned\"; print aligned + unaligned > \"${outdir}/${output}.reads.all\" } "

