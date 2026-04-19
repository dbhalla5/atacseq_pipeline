#!/bin/bash

# =========================
# TOBIAS BINDetect: merge replicates + run BINDetect + combine TF activity
# =========================

# Activate environment
# conda activate deeptools_env

ulimit -n 8192


# ------------------------
# Paths
# ------------------------
GENOME="/path/to/genomes/hg38/hg38.fa"
MOTIFS="/path/to/motifDBs/JASPAR_2024.meme"

BINDIR="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/BINDetect"
DAR_DIR="/path/to/ATAC/downstream_analysis/DiffBind_DARs"

# FootprintScores files per replicate

DMSO_REP2_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/DMSO_REP2/FootprintScores/DMSO_REP2_footprints.bw"
DMSO_REP3_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/DMSO_REP3/FootprintScores/DMSO_REP3_footprints.bw"

Treatment_1_REP2_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_1_REP2/FootprintScores/Treatment_1_REP2_footprints.bw"
Treatment_1_REP3_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_1_REP3/FootprintScores/Treatment_1_REP3_footprints.bw"

Treatment_1_2_REP2_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_1_2_REP2/FootprintScores/Treatment_1_2_REP2_footprints.bw"
Treatment_1_2_REP3_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_1_2_REP3/FootprintScores/Treatment_1_2_REP3_footprints.bw"

Treatment_2_REP2_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_2_REP2/FootprintScores/Treatment_2_REP2_footprints.bw"
Treatment_2_REP3_BW="/path/to/ATAC/downstream_analysis/tobias_res/day_2_3_only/Treatment_2_REP3/FootprintScores/Treatment_2_REP3_footprints.bw"

# DAR BED files

#DAR_Treatment_1_Treatment_1_2="$DAR_DIR/Day_2_3_only_summits_100_Diffbind_DARs_Treatment_1_and_Treatment_1_2_merged.bed"
#DAR_Treatment_2_Treatment_1_2="$DAR_DIR/Day_2_3_only_summits_100_Diffbind_DARs_Treatment_2_and_Treatment_1_2_merged.bed"

DAR_DMSO_Treatment_1="$DAR_DIR/sorted_DiffBind_DARs_Day_2_3_only_summits_100_DMSO_vs_Treatment_1.bed"
DAR_DMSO_Treatment_2 ="$DAR_DIR/sorted_DiffBind_DARs_Day_2_3_only_summits_100_DMSO_vs_Treatment_2.bed"
DAR_DMSO_Treatment_1_2 ="$DAR_DIR/sorted_DiffBind_DARs_Day_2_3_only_summits_100_DMSO_vs_Treatment_1_2.bed"



# Number of cores
CORES=4

# Create BINDetect output directory
# mkdir -p "$BINDIR"

# ------------------------
# Step 1: Merge replicates using deepTools bigwigCompare
# ------------------------
#  echo "Merging Treatment_1 replicates..."
#  bigwigCompare -b1 "$Treatment_1_REP2_BW" -b2 "$Treatment_1_REP3_BW" --operation mean -o "$BINDIR/Treatment_1_merged.bw"

#  echo "Merging Treatment_1_2 replicates..."
#  bigwigCompare -b1 "$Treatment_1_2_REP2_BW" -b2 "$Treatment_1_2_REP3_BW" --operation mean -o "$BINDIR/Treatment_1_2_merged.bw"

#  echo "Merging Treatment_2 replicates..."
#  bigwigCompare -b1 "$Treatment_2_REP2_BW" -b2 "$Treatment_2_REP3_BW" --operation mean -o "$BINDIR/Treatment_2_merged.bw"

#  echo "Merging DMSO replicates..."
#  bigwigCompare -b1 "$DMSO_REP2_BW" -b2 "$DMSO_REP3_BW" --operation mean -o "$BINDIR/DMSO_merged.bw"

# ------------------------
# Step 2: Run BINDetect
# ------------------------
run_bindetect() {
    local NAME=$1
    local PEAKS=$2
    local SIG1=$3
    local SIG2=$4
    local COND1=$5
    local COND2=$6

    OUTDIR="$BINDIR/$NAME"
    mkdir -p "$OUTDIR"

    echo "Running BINDetect for $NAME ..."
    TOBIAS BINDetect \
        --motifs "$MOTIFS" \
        --signals "$SIG1" "$SIG2" \
        --peaks "$PEAKS" \
        --genome "$GENOME" \
        --cond_names "$COND1" "$COND2" \
        --outdir "$OUTDIR" \
        --cores "$CORES" \
        --skip_excel True \
        --time_series False

    echo "Finished BINDetect for $NAME"
}

# 1. DMSO vs Treatment_1 DARs
run_bindetect "DMSO_vs_Treatment_1_DARs" "$DAR_DMSO_Treatment_1" \
    "$BINDIR/DMSO_merged.bw" "$BINDIR/Treatment_1_merged.bw" \
    DMSO Treatment_1

# 2. DMSO vs Treatment_2 DARs
run_bindetect "DMSO_vs_Treatment_2_DARs" "$DAR_DMSO_Treatment_2" \
    "$BINDIR/DMSO_merged.bw" "$BINDIR/Treatment_2_merged.bw" \
    DMSO Treatment_2
    
# 3. DMSO vs Treatment_1_2 DARs

run_bindetect "DMSO_vs_Treatment_1_2_DARs" "$DAR_DMSO_Treatment_1_2" \
    "$BINDIR/DMSO_merged.bw" "$BINDIR/Treatment_1_2_merged.bw" \
    DMSO Treatment_1_2

