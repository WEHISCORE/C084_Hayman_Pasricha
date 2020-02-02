#!/usr/bin/env bash
# Run FastQC (and MultiQC on the output).
# Peter Hickey
# 2019-11-12

module load fastqc

OUTDIR="../output/FastQC"
mkdir -p ${OUTDIR}

fastqc -o ${OUTDIR} \
  --threads 10 \
  ../extdata/NN179/*fastq.gz

multiqc --title NN179 \
  --outdir ${OUTDIR} \
  --no-megaqc-upload \
  ${OUTDIR}
