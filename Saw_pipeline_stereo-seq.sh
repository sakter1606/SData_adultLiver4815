
#!/usr/bin/env bash

#!/usr/bin/env bash
1.User‑tunable variables
MODE="${MODE:-Q40}"                   # Q40 or Q4
SN="${SN:-SS200000000TL_X1}"          # chip serial number
SPECIES="${SPECIES:-cattle}"           # species name for CIDCount
GENOME_SIZE="${GENOME_SIZE:-2736692181}" # e.g., mouse ~2.7e9; set appropriately
THREADS="${THREADS:-24}"

DATA_DIR="/data"    REF_DIR="/ref"    OUT_DIR="/work/result"  SIF="/env/SAW_v7.0.sif"

MASK_H5="$DATA_DIR/mask/${SN}.barcodeToPos.h5"
FASTA="$REF_DIR/genome/genome.fa"
GTF="$REF_DIR/genes/genes.gtf"
STAR_DIR="$REF_DIR/STAR_SJ100"

FQ1_LIST="/data/reads/L01_read_1.fq.gz,/data/reads/L02_read_1.fq.gz"
FQ2_LIST="/data/reads/L01_read_2.fq.gz,/data/reads/L02_read_2.fq.gz"
export SINGULARITY_BIND="$DATA_DIR,$REF_DIR,$OUT_DIR"

log() { echo "[$(date +'%F %T')] $*"; }
##Build the STAR index
build_star_index() {
  mkdir -p "$STAR_DIR"
  if [ ... index exists ... ]; then skip; else
    singularity exec "$SIF" mapping \
      --runMode genomeGenerate \
      --genomeDir "$STAR_DIR" \
      --genomeFastaFiles "$FASTA" \
      --sjdbGTFfile "$GTF" \
      --sjdbOverhang 99 \
      --runThreadN "$THREADS"
  fi
}
 ##--- STEP 1: Q4-only splitMask ---
split_mask_q4() {
  local SPLIT_DIR="$OUT_DIR/00.splitmask"
  mkdir -p "$SPLIT_DIR"
  log "Running splitMask for Q4 → $SPLIT_DIR"
  singularity exec "$SIF" splitMask \
      "$MASK_H5" \
      "$SPLIT_DIR"

#CIDCount (estimate CID & memory) ---
cidcount() {
  INPUT_BIN_OR_H5=... # Q4 uses a split .bin, Q40 uses .h5
  singularity mapping (STAR alignment inside SAW)
  exec "$SIF" CIDCount -i "$INPUT_BIN_OR_H5" -s "$SPECIES" -g "$GENOME_SIZE" > "$OUT_DIR/00.cidcount.stat"

##7) mapping (STAR alignment inside SAW)
#Q40 mapping

mapping_q40() {
  # Make a bcPara for each lane (Q40)
  in=$MASK_H5
  in1=<lane_read_1.fq.gz>
  in2=<lane_read_2.fq.gz> (optional)
  barcodeReadsCount=<lane>.barcodeReadsCount.txt
  barcodeStart=0
  barcodeLen=25
  umiStart=25
  umiLen=10
  mismatch=1
  polyAnum=15
  mismatchInPolyA=2

  # Run 'mapping' with STAR flags (same for each lane)
  singularity exec "$SIF" mapping \
    --outSAMattributes spatial \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir "$STAR_DIR" \
    --runThreadN "$THREADS" \
    --outFileNamePrefix "$MAP_DIR/${lane}." \
    --sysShell /bin/bash \
    --stParaFile "$MAP_DIR/${lane}.bcPara" \
    --readNameSeparator " " \
    --limitBAMsortRAM 38582880124 \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize=280000000 \
    --outBAMsortingBinsN 50 \
    --outSAMmultNmax -1
}
##Q4 mapping

mapping_q4() {
  # Per-index bcPara (e.g., 01..16)
  in=$OUT_DIR/00.splitmask/${idx}.${SN}.barcodeToPos.bin
  in1=$DATA_DIR/reads/L01_${idx}.fq.gz
  barcodeReadsCount=$MAP_DIR/${idx}.barcodeReadsCount.txt
  barcodeLen=24  # Q4 uses 24
  umiLen=10
  ...
  # mapping command identical to Q40, different bcPara inputs
}
#merge (only if multiple mapping outputs)

merge_cid_lists() {
  singularity exec "$SIF" merge \
    "$MASK_H5" \
    "$MAP_DIR/"*.barcodeReadsCount.txt \
    "$OUT_DIR/02.merge/${SN}.merge.barcodeReadsCount.txt"
}
# count (dedup, annotate, create GEF)

count_all() {
  singularity exec "$SIF" count \
    -i $BAM_GLOB \
    -o "$COUNT_DIR/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam" \
    -a "$GTF" \
    -s "$COUNT_DIR/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat" \
    --umi_len 10 \
    --sn "$SN" \
    --umi_on \
    --save_lq \
    -c "$THREADS" \
    --multi_map
##

tissue_cut() {
  singularity exec "$SIF" tissueCut \
    -i "$OUT_DIR/03.count/${SN}.raw.gef" \
    -o "$OUT_DIR/05.tissuecut" \
    -O Transcriptomics \
    -d \
    --sn "$SN"
}

spatial_cluster() {
  singularity exec "$SIF" spatialCluster \
    -i "$OUT_DIR/05.tissuecut/${SN}.tissue.gef" \
    -r 1.0 \
    -s 200
}
