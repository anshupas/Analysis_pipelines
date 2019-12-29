#!/bin/bash

## sets PIPELINE_HOME
source /package/sequencer/ngspipeline/profile

if [ $# -eq 3 ]; then
    FASTQ_FILE="$1"
    SAMPLE="$2"
    REF="$3"
else
  echo "Usage: $0 <FASTQ_FILE> <SAMPLE> <REFERENCE_NAME>"
  exit
fi

set -e
#set -x

:<<'RRBS_WORKFLOW'

  DESCRIBE

RRBS_WORKFLOW

## init RRBS pipeline :-)
DATE_FLAG=$(date +%Y-%m-%d)

## config
WORKFLOW="RRBS"
WORK_DIR=$(pwd)/${WORKFLOW}_${DATE_FLAG}_${SAMPLE}
GENOME="/project/genomes/${REF}/Sequence/WholeGenomeFasta/genome.fa"
FASTQ2_FILE="${FASTQ_FILE/_R1/_R2}"

FASTQC_DIR="00_FastQC"
FASTQ1_PATH=$(realpath $FASTQ_FILE)
FASTQ2_PATH=$(realpath $FASTQ2_FILE)

CURRENT_DIR=$(pwd)


## software
#FASTQC="/package/sequencer/fastqc/current/fastqc"
FASTQC="/project/bioinf_meissner/src/fastQC/fastqc_v0.11.5/fastqc"
JAVA_BIN="/package/sequencer/java/8/bin/java"
PICARD_TOOLS="/package/sequencer/picard-tools/current/picard.jar"
BSMAP="/package/sequencer/bin/bsmap"
SAMTOOLS="/package/sequencer/bin/samtools"
MCALL="/package/sequencer/ngspipeline/bin/mcall"

PIGZ="/package/sequencer/bin/pigz"

## computing
PROCESS_THREADS=8
MAPPING_THREADS=20

if [ ! -r $GENOME ]; then
    echo "$REF not found ($GENOME)"
    exit
fi

## here we go
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -s $FASTQ1_PATH "${SAMPLE}_R1_00_input.fq.gz"
ln -s $FASTQ2_PATH "${SAMPLE}_R2_00_input.fq.gz"

## (1) Fastqc
IN_FILE1="${SAMPLE}_R1_00_input.fq.gz"
IN_FILE2="${SAMPLE}_R2_00_input.fq.gz"
mkdir -p $FASTQC_DIR
$FASTQC --outdir $FASTQC_DIR --threads $PROCESS_THREADS $IN_FILE1 $IN_FILE2


## (2) Trimming

##  (a) Illumina Adapter Trimming
## long : AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
## short: AGATCGGAAGAGC
ADAPT1="AGATCGGAAGAGC"
ADAPT2="AAATCAAAAAAAC"
ADAPT3="CGGAAGAGCACAC"
BIN=/package/sequencer/anaconda3/envs/ngs/bin/cutadapt

OUT_FILE1="${SAMPLE}_R1_01_adp-trim.fq.gz"
OUT_FILE2="${SAMPLE}_R2_01_adp-trim.fq.gz"
$BIN \
  --quality-cutoff 20 \
  --overlap 3 \
  --minimum-length 25 \
  -u 3 \
  -U 3 \
  -a $ADAPT1 \
  -a $ADAPT2 \
  -a $ADAPT3 \
  -A $ADAPT1 \
  -A $ADAPT2 \
  -A $ADAPT3 \
  --cores $PROCESS_THREADS \
  -o $OUT_FILE1 \
  -p $OUT_FILE2 \
  $IN_FILE1 \
  $IN_FILE2 \
  &> "${OUT_FILE1%_R1_01_adp-trim.fq.gz}_01_adp-trim.log"


## (x) New Fastqc after Trimming
IN_FILE1=$OUT_FILE1
IN_FILE2=$OUT_FILE2
$FASTQC --outdir $FASTQC_DIR --threads $PROCESS_THREADS $IN_FILE1 $IN_FILE2


## (3) Mapping
IN_FILE1=$OUT_FILE1
IN_FILE2=$OUT_FILE2
OUT_FILE="${SAMPLE}_03_bsmap_rrbs.bam"
LOG="${OUT_FILE}.log"

# -s = seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
# In pipeline this was set to 1 (WGBS mode); here I revert it back to default, 12
#
# -u = report unmapped reads, default=off
# -R = print corresponding reference sequences in SAM output, default=off
# -m = minimal insert size allowed, default=28 (here set to 0 to allow overlapping mates)
# removed -D C-CGG \
($BSMAP \
  -v 0.1 -s 12 -q 20 -w 100 -S 1 -R -u -m 0 \
  -p $MAPPING_THREADS \
  -d $GENOME \
  -a $IN_FILE1 \
  -b $IN_FILE2 \
  | $SAMTOOLS sort -m 4G -@ $PROCESS_THREADS -O BAM -o $OUT_FILE - ) 2>> $LOG

$SAMTOOLS index $OUT_FILE


## (4) Statistics
IN_FILE=$OUT_FILE # BAM
OUT_TEXT="${IN_FILE}.metrics.txt"
$JAVA_BIN \
  -Xmx14g \
  -jar $PICARD_TOOLS CollectAlignmentSummaryMetrics \
  INPUT=$IN_FILE \
  OUTPUT=$OUT_TEXT \
  TMP_DIR=$(pwd)/tmp.java \
  REFERENCE_SEQUENCE=$GENOME \
  IS_BISULFITE_SEQUENCED=True

rm -fr "tmp.java"


## (5) Methylation Calling
BAM_FILE=$OUT_FILE

PREFIX="${SAMPLE}_05_mcall_rrbs" # AM-RRBS-085_03_bsmap_rrbs.bam.G.bed
$MCALL \
  --threads $MAPPING_THREADS \
  --reference $GENOME \
  --sampleName $PREFIX \
  --mappedFiles $BAM_FILE \
  --outputDir $WORK_DIR \
  --webOutputDir $WORK_DIR

## AM-RRBS-085_05_mcall_rrbs.G.bed -> AM-RRBS-085_03_bsmap_rrbs.bam.G.bed
rm -f "${PREFIX}.G.bed"
rm -f "${BAM_FILE}.HG.bed"
mv "${BAM_FILE}.G.bed" "${PREFIX}.CpG.bed"
rm -f run.config.*

## BigBed Conversion
## FIX -- SUBSTITUTE !!
CONVERT="/package/sequencer/ngspipeline/open_epp/convertmCallBigBed.py"
python $CONVERT \
  "$WORK_DIR/${PREFIX}.CpG.bed" \
  "${SAMPLE}_05_mcall_rrbs.CpG.bb" \
  "/project/genomes/${REF}/Sequence/WholeGenomeFasta/genome.chrom.sizes"


banner DONE
