#### ATAC-seq data downstream analysis script
### 
### original sequencing file ID: ATAC-seq_SLX-12345


# conda install -c bioconda bioconductor-biomaRt
# conda install bioconda::bioconductor-diffbind

#library(remoter)
#remoter::client("address.of.cluster")

#remoter::client("localhost", port = 55555, verbose = TRUE)

# Local script for visualizing and analyzing DiffBind results

# ATAC-seq downstream analysis script using DiffBind-analyzed object
# Load libraries
# === Load Libraries ===

# ===== ATAC-seq downstream analysis using DiffBind =====

# basic plots:
#   
#  after reading in

# 1. PCA plot
# 
# > dba.plotPCA(dba_obj,
#               +             attributes=DBA_CONDITION,
#               +             label=DBA_CONDITION,
#               +             vColors=c("chartreuse4", "cornflowerblue", "darkorange1","black"),
#               +             labelSize=0.8,
#               +             dotSize=1.5)


# 2. heat map

# dba.plotHeatmap(dba_obj, attributes=c(DBA_CONDITION, DBA_REPLICATE))


# 3. Volcano plot
# 
# > dba.plotVolcano(dba_obj,contrast =1) # DMSO vs Treatment_1 
# > dba.plotVolcano(dba_obj,contrast =2) # DMSO vs Treatment_2 
# > dba.plotVolcano(dba_obj,contrast =3) # DMSO vs Treatment_1_2 


# 4. MA plot

# dba.plotMA(dba_obj, contrast=1)  # DMSO vs Treatment_1 
# dba.plotMA(dba_obj, contrast=2)  # DMSO vs Treatment_2 
# dba.plotMA(dba_obj, contrast=3)  # DMSO vs Treatment_1_2 

# 5. box plots (optional)

# dba.plotBox(dba_obj, contrast =1, method=dba_obj$config$AnalysisMethod, th=0.05, bUsePval = T)

# ============================ #
# DiffBind Downstream Analysis
# ============================ #
# === Libraries ===


library(DiffBind)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(tidyverse)
library(dplyr)

library(clusterProfiler)
library(org.Hs.eg.db)  # For human gene annotations
library(DOSE)          # For enrichment visualization


# === Setup ===
setwd("/path/to/ATAC_analysis")
dir.create("plots", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

dba_obj <- readRDS("path/to/rds/dba_d2d3_summits_250_analyzed_day23_only.rds")
prefix <- "Day_2_3"

#dba_obj <- readRDS("result_files/rds/dba_all_summits_250_analyzed_all_days.rds")
#prefix <- "all_days"


# Clear broken peak file paths
dba_obj$samples$peaks <- rep(NA_character_, nrow(dba_obj$samples))

padj_cutoff <- 0.05  # padj < 0.05
lfc_cutoff <- 1      # |log2FC| >1

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

report <- dba.report(dba_obj)

annotated_peaks_t <- annotatePeak(report, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
peak_anno_df <- as.data.frame(annotated_peaks_t)
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")  # For human


####### Differential accessibility analysis and Peak Annotation



for (i in seq_along(dba_obj$contrasts)) {
  contrast_name <- paste(dba_obj$contrasts[[i]]$name2, "vs", dba_obj$contrasts[[i]]$name1, sep = "_")
  message("Processing: ", contrast_name)
  
  sig_peaks <- dba.report(dba_obj, contrast = i, th = padj_cutoff, fold = lfc_cutoff)
  if (length(sig_peaks) == 0) {
    message("No significant peaks in contrast ", contrast_name)
    next
  }
  
  # === Add Peak IDs ===
  peak_ids <- names(sig_peaks)
  if (is.null(peak_ids)) peak_ids <- paste0("peak_", seq_along(sig_peaks))
  names(sig_peaks) <- peak_ids
  
  # Prepare DE statistics for merging
  stats_df <- as.data.frame(mcols(sig_peaks))[, c("Fold", "FDR", "p.value")]
  stats_df$peak_id <- peak_ids
  
  # === Annotation with ChIPseeker ===
  peakAnno <- annotatePeak(sig_peaks, TxDb = txdb, annoDb = "org.Hs.eg.db")
  anno_df <- as.data.frame(peakAnno)
  anno_df$peak_id <- rownames(anno_df)
  
  # === Save ChIPseeker Annotation Plots into Single PDF ===
  pdf(file = paste0("plots/Annotation_plots_", prefix, "_", contrast_name, ".pdf"),
      width = 10, height = 8)
  
  # 1. Barplot of annotation categories
  bar_plot <- plotAnnoBar(peakAnno)
  print(bar_plot)
  
  # 2. UpSet plot with Venn pie
  upset_plot <- upsetplot(peakAnno, vennpie = TRUE)
  print(upset_plot)
  
  # 3. Distance to TSS plot (base R, no print() needed)
  TSSplot <- plotDistToTSS(peakAnno,
                title = "Distribution of transcription factor-binding loci\nrelative to TSS")
  print(TSSplot)
  
  dev.off()
  

  
  
  # === Transcript ID Cleanup ===
  if ("transcriptId" %in% colnames(anno_df)) {
    anno_df$transcript_clean <- sub("\\.\\d+$", "", trimws(anno_df$transcriptId))
  }
  
  # === Merge Annotation + Stats ===
  anno_merged <- left_join(anno_df, stats_df, by = "peak_id")
  
  # Rename for clarity
#  colnames(anno_merged)[colnames(anno_merged) == "Fold"] <- "log2FoldChange"
#  colnames(anno_merged)[colnames(anno_merged) == "p.value"] <- "pvalue"
  colnames(anno_merged)[colnames(anno_merged) == "Entrez_Gene_ID"] <- "geneId"
  
  
  # === Output CSV ===
  write.csv(anno_merged,
            file = paste0("tables/Annotated_sigPeaks_", prefix, "_", contrast_name, ".csv"),
            row.names = FALSE)
  
}


#####################################################################################################

## GO and KEGG enrichment

annot_file <- "tables/Annotated_sigPeaks_Day_2_3_Treatment_1_2_vs_DMSO_summits_250.csv"

# Explicitly define contrast_name to avoid stale values
contrast_name <- "Treatment_1_2_vs_DMSO"

anno_df <- read.csv(annot_file)
entrez_ids <- unique(anno_df$geneId)
entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]

# Continue with enrichment analysis and plotting...



### GO
# Check and filter gene list
gene_list <- anno_df$geneId[!is.na(anno_df$geneId)]
gene_list <- unique(gene_list)

if (length(gene_list) < 10) {
  message("Not enough genes for enrichment in ", contrast_name)
  next
}

# GO BP Enrichment
ego <- tryCatch({
  enrichGO(gene = gene_list,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           readable = TRUE)
}, error = function(e) {
  message("GO enrichment failed for ", contrast_name, ": ", e$message)
  return(NULL)
})

# KEGG Enrichment
ekegg <- tryCatch({
  enrichKEGG(gene = gene_list,
             organism = "hsa",
             keyType = "kegg",
             pvalueCutoff = 0.05)
}, error = function(e) {
  message("KEGG enrichment failed for ", contrast_name, ": ", e$message)
  return(NULL)
})

# === Plotting ===


  barplot(ego, showCategory = 20, title = paste("GO Enrichment -", contrast_name))
  dotplot(ego, showCategory = 20, title = paste("GO Enrichment -", contrast_name))
  
 ##### fallback plan. if the above doesn't work. 
  
  #### bar plot
  library(ggplot2)
  
  # Simplify top GO terms
  top_terms <- head(ego@result[order(ego@result$pvalue), ], 10)
  
  ggplot(top_terms, aes(x = reorder(Description, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Top GO Enrichment -", contrast_name),
         x = "GO Term", y = "Gene Count") +
    theme_minimal()
  
  
  
  ######### dot plot
  
  library(ggplot2)
  
  # Take top 10 GO terms by p-value
  top_terms <- head(ego@result[order(ego@result$pvalue), ], 10)
  
  # Convert GeneRatio from string to numeric (e.g., "5/277" → 5/277)
  top_terms$GeneRatioNum <- sapply(top_terms$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  # Plot
  ggplot(top_terms, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
    scale_size(range = c(3, 8)) +
    labs(
      title = paste("GO Enrichment Dot Plot -", contrast_name),
      x = "Gene Ratio", y = "GO Term"
    ) +
    theme_minimal()
  
  #########################
  
  
  
  

  barplot(ekegg, showCategory = 20, title = paste("KEGG Enrichment -", contrast_name))
  dotplot(ekegg, showCategory = 20, title = paste("KEGG Enrichment -", contrast_name))
  
  
  # Simplify top KEGG terms
  top_terms <- head(ekegg@result[order(ekegg@result$pvalue), ], 10)
  
  ggplot(top_terms, aes(x = reorder(Description, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Top KEGG Enrichment -", contrast_name),
         x = "KEGG Term", y = "Gene Count") +
    theme_minimal()
  
  
  # Plot
  # Take top 10 KEGG terms by p-value
  top_terms <- head(ekegg@result[order(ekegg@result$pvalue), ], 10)
  
  # Convert GeneRatio from string to numeric (e.g., "5/277" → 5/277)
  top_terms$GeneRatioNum <- sapply(top_terms$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  
  ggplot(top_terms, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
    scale_size(range = c(3, 8)) +
    labs(
      title = paste("KEGG Enrichment Dot Plot -", contrast_name),
      x = "Gene Ratio", y = "KEGG Term"
    ) +
    theme_minimal()
  
  
  
