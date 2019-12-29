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

## init RRBS pipeline
DATE_FLAG=$(date +%Y-%m-%d)

## config
WORKFLOW="RRBS"
WORK_DIR=$(pwd)/${WORKFLOW}_${DATE_FLAG}_${SAMPLE}
GENOME="/project/genomes/${REF}/Sequence/WholeGenomeFasta/genome.fa"
FASTQC_DIR="00_FastQC"
FASTQ_PATH=$(realpath $FASTQ_FILE)
CURRENT_DIR=$(pwd)

FASTQ_INDEX_FILE="${FASTQ_FILE/_R1_/_R2_}"
FASTQ_INDEX_PATH=$(realpath $FASTQ_INDEX_FILE)


## software
FASTQC="/package/sequencer/fastqc/current/fastqc"
JAVA_BIN="/package/sequencer/java/8/bin/java"
PICARD_TOOLS="/package/sequencer/picard-tools/current/picard.jar"
NUGEN_TRIM="/package/sequencer/ngspipeline/deps/nugen/trimRRBSdiversityAdaptCustomers.py"
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
ln -s $FASTQ_PATH "${SAMPLE}_00_input.fq.gz"

## (1) Fastqc
IN_FILE="${SAMPLE}_00_input.fq.gz"
mkdir -p $FASTQC_DIR
$FASTQC --outdir $FASTQC_DIR --threads $PROCESS_THREADS $IN_FILE


## (2) Trimming

##  (a) Illumina Adapter Trimming
## long : AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
## short: AGATCGGAAGAGC
ADAPT1="AGATCGGAAGAGC"
ADAPT2="AAATCAAAAAAAC"
BIN=/package/sequencer/anaconda3/envs/ngs/bin/cutadapt
OUT_FILE="${SAMPLE}_01_adp-trim.fq.gz"
$BIN \
  --quality-cutoff 20 \
  --overlap 3 \
  --minimum-length 25 \
  --adapter $ADAPT1 \
  --adapter $ADAPT2 \
  --cores $PROCESS_THREADS \
  --output $OUT_FILE \
  $IN_FILE \
  &> "${OUT_FILE%.fq.gz}.log"


## (b) Nugene Diversity Adapter Trimming
IN_FILE=$OUT_FILE
OUT_FILE="${SAMPLE}_02_div-trim.fq.gz"
LOG_FILE="${SAMPLE}_02_div-trim.log"

python \
  $NUGEN_TRIM \
  -1 $IN_FILE \
  -o $LOG_FILE

mv "${IN_FILE%.gz}_trimmed.fq.gz" $OUT_FILE



## (x) New Fastqc after Trimming
IN_FILE=$OUT_FILE
$FASTQC --outdir $FASTQC_DIR --threads $PROCESS_THREADS $IN_FILE


## (4) Mapping
IN_FILE=$OUT_FILE
OUT_FILE="${SAMPLE}_03_bsmap_rrbs.bam"

# -s = seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
# In pipeline this was set to 1 (WGBS mode); here I revert it back to default, 12
#
# -u = report unmapped reads, default=off
# -R = print corresponding reference sequences in SAM output, default=off
$BSMAP \
  -v 0.1 -s 12 -q 20 -w 100 -S 1 -u -R \
  -p $MAPPING_THREADS \
  -d $GENOME \
  -D C-CGG \
  -a $IN_FILE \
  | $SAMTOOLS sort -m 4G -@ $PROCESS_THREADS -O BAM -o $OUT_FILE -

$SAMTOOLS index $OUT_FILE


## (5) Statistics
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

## (xx) Deduplication Nugen
## samtools view -f 0x900 AM-RRBS-085_03_bsmap_rrbs.bam
## must return nothing. This is 0x100 (secondary) and 0x800 (supplementary)
mkdir -p $(pwd)/tmp.nugene
BIN=/package/sequencer/ngspipeline/deps/nugen/nugentechnologies-nudup-468c62e/nudup.py
IN_FILE=$OUT_FILE # BAM
PREFIX="${SAMPLE}_04_dedup"
$BIN \
  -f $FASTQ_INDEX_PATH \
  --out $PREFIX \
  --start 6 \
  --length 6 \
  -T $(pwd)/tmp.nugene \
  $IN_FILE

rm -fr $(pwd)/tmp.nugene

## result = AM-RRBS-085_04_dedup.sorted.dedup.bam

## (6o7) Methylation Calling
BAM_FILE=${PREFIX}.sorted.dedup.bam
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
