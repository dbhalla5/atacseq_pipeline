#!/usr/bin/env Rscript

# Script: prepare_diffbind_rds.R  (for cluster)
# Purpose: Run DiffBind normalization + DE analysis for ATAC-seq 
# Data:  ATAC-seq (SLX-12345)


# to run this script
# 1. start a tmux session
# 2. activate conda (conda activate R-433)
# 3. Rscript prepare_diffbind_rds_on_cluster.R > log_diffbind.txt 2>&1

suppressPackageStartupMessages({
  library(DiffBind)
  library(tidyverse)
  library(biomaRt)
  library(BiocParallel)
})

# --- Setup ---
setwd("/path/to/ATAC/downstream_analysis")
sample_sheet <- "diffbind_sample_sheet_ATAC.csv"
samples <- read.csv(sample_sheet, stringsAsFactors = FALSE)
dir.create("rds", showWarnings = FALSE)

# --- Parallel params ---
param <- MulticoreParam(workers = 12)

# --- Load full dataset ---
dba_all <- dba(sampleSheet = samples)
dba_all <- dba.count(dba_all, summits = 250, bUseSummarizeOverlaps = TRUE, bParallel = TRUE)
# summits can be set to 75, if needed.  will have 151 bp window (summit - 75 to summit +75)
# this is based on the fragment size in the MultiQC file.

# --- Set config for DESeq2 and normalization ---
dba_all$config$Normalization <- DBA_NORM_NATIVE
dba_all$config$AnalysisMethod <- DBA_DESEQ2

# --- Explicit global contrasts: DMSO vs all ---
# Set the design formula
dba_all <- dba.contrast(dba_all, design = ~Condition)

# Add each contrast with the factor name included in contrast vector (no need for name1, name2)
dba_all <- dba.contrast(dba_all, contrast = c("Condition", "DMSO", "Treatment_1"))
dba_all <- dba.contrast(dba_all, contrast = c("Condition", "DMSO", "Treatment_2"))
dba_all <- dba.contrast(dba_all, contrast = c("Condition", "DMSO", "Treatment_1_2"))
# --- Normalize and analyze ---
dba_all <- dba.normalize(dba_all, method = DBA_DESEQ2, normalize = DBA_NORM_NATIVE, background = TRUE)
dba_all <- dba.analyze(dba_all, bParallel = TRUE)

# --- Save full analyzed object ---
saveRDS(dba_all, file = "rds/dba_all_summits_250_analyzed_all_days.rds")


# === Subset to Day2 + Day3 only ===
samples_d2d3 <- subset(samples, Day %in% c("Day2", "Day3"))
dba_d2d3 <- dba(sampleSheet = samples_d2d3)
dba_d2d3 <- dba.count(dba_d2d3, summits = 250, bUseSummarizeOverlaps = TRUE, bParallel = TRUE) 
# summits can be set to 75, if needed.  will have 151 bp window (summit - 75 to summit +75)
# this is based on the fragment size in the MultiQC file.

# Config
dba_d2d3$config$Normalization <- DBA_NORM_NATIVE
dba_d2d3$config$AnalysisMethod <- DBA_DESEQ2

# Explicit global contrasts again


# Define model matrix for Day2 + Day3 subset
dba_d2d3 <- dba.contrast(dba_d2d3, design = ~Condition)

# Add DMSO vs treatment contrasts with the factor name included
dba_d2d3 <- dba.contrast(dba_d2d3, contrast = c("Condition", "DMSO", "Treatment_1"))
dba_d2d3 <- dba.contrast(dba_d2d3, contrast = c("Condition", "DMSO", "Treatment_2"))
dba_d2d3 <- dba.contrast(dba_d2d3, contrast = c("Condition", "DMSO", "Treatment_1_2"))


# Normalize + analyze
dba_d2d3 <- dba.normalize(dba_d2d3, method = DBA_DESEQ2, normalize = DBA_NORM_NATIVE, background = TRUE)
dba_d2d3 <- dba.analyze(dba_d2d3, bParallel = TRUE)

# Save Day2+3 object
saveRDS(dba_d2d3, file = "rds/dba_d2d3_summits_250_analyzed_day23_only.rds")

message("All .rds files with DMSO-vs-all contrasts saved.")
