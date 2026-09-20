# ============================================================
# Core analysis: DESeq2 differential expression + global VST
# ------------------------------------------------------------
# Produces:
#   ../output/DESeq2_output.xlsx  (6 pairwise comparisons vs. control;
#                                  matches Supplementary Table 1)
#   ../output/vst_matrix.xlsx     (intermediate file; input to 02_wgcna.R)
#
# Run from a "code" folder that sits next to "input" and "output":
#   input/counts.xlsx, input/metadata.xlsx
#   code/01_deseq2.R   (this script)
#   output/            (created automatically)
# ============================================================

# ---- paths ----
IN_DIR  <- "../input"
OUT_DIR <- "../output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COUNTS_FILE <- file.path(IN_DIR, "counts.xlsx")
META_FILE   <- file.path(IN_DIR, "metadata.xlsx")

# ---- packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(DESeq2)
})
options(stringsAsFactors = FALSE)
set.seed(42)

# ---- group definitions ----
# 13 groups / 39 samples used for the global VST matrix
GROUP_ORDER_VST13 <- c("ctrl",
                       "cru_om", "pur_om",
                       "cru_pf", "pur_pf",
                       "cru_px", "pur_px",
                       "aav", "lps", "ec", "pa",
                       "iav", "iav_hk")

# 7 groups / 21 samples used to fit the DESeq2 model
GROUP_ORDER_DESEQ7 <- c("ctrl",
                        "cru_om", "pur_om",
                        "cru_pf", "pur_pf",
                        "cru_px", "pur_px")

# ---- read data ----
raw       <- read_excel(COUNTS_FILE)
meta_full <- read_excel(META_FILE)

stopifnot("gene_id" %in% colnames(raw))
count_full <- raw %>% column_to_rownames("gene_id") %>% dplyr::select(-gene_name)
gene_anno  <- raw %>% dplyr::select(gene_id, gene_name) %>% distinct()

# ---- sample sets ----
meta_vst <- meta_full %>%
  filter(group %in% GROUP_ORDER_VST13) %>%
  mutate(group = factor(group, levels = GROUP_ORDER_VST13)) %>%
  arrange(group) %>%
  as.data.frame()
rownames(meta_vst) <- meta_vst$sample
count_vst <- count_full[, meta_vst$sample, drop = FALSE]
message("VST sample set: ", ncol(count_vst), " samples (expect 39)")

meta_deseq <- meta_full %>%
  filter(group %in% GROUP_ORDER_DESEQ7) %>%
  mutate(group = factor(group, levels = GROUP_ORDER_DESEQ7)) %>%
  arrange(group) %>%
  as.data.frame()
rownames(meta_deseq) <- meta_deseq$sample
count_deseq <- count_full[, meta_deseq$sample, drop = FALSE]
message("DESeq2 sample set: ", ncol(count_deseq), " samples (expect 21)")

# ---- shared pre-filter, computed on the 39-sample set ----
keep_genes      <- rowSums(count_vst >= 5) >= 3
message("Genes retained after pre-filter: ", sum(keep_genes))
count_vst_flt   <- count_vst[keep_genes, ]
count_deseq_flt <- count_deseq[keep_genes, ]

# ────────────────────────────────────────────────────────────
# DESeq2: 7-group model (21 samples), 6 comparisons vs. control
# ────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = count_deseq_flt,
                              colData   = meta_deseq,
                              design    = ~ group)
dds$group <- relevel(dds$group, ref = "ctrl")
dds <- DESeq(dds)
message("DESeq2 model fit complete")

comparisons <- list(
  list(name = "Crd.OMKO_vs_Ctrl", num = "cru_om", ref = "ctrl"),
  list(name = "Crd.Pf_vs_Ctrl",   num = "cru_pf", ref = "ctrl"),
  list(name = "Crd.PhiX_vs_Ctrl", num = "cru_px", ref = "ctrl"),
  list(name = "Pur.OMKO_vs_Ctrl", num = "pur_om", ref = "ctrl"),
  list(name = "Pur.Pf_vs_Ctrl",   num = "pur_pf", ref = "ctrl"),
  list(name = "Pur.PhiX_vs_Ctrl", num = "pur_px", ref = "ctrl")
)
LFC_THRESH <- 0.585   # padj < 0.05 AND |log2FC| > 0.585

results_list <- list()
for (comp in comparisons) {
  res <- results(dds,
                 contrast = c("group", comp$num, comp$ref),
                 alpha    = 0.05,
                 pAdjustMethod = "BH")
  df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_anno, by = "gene_id") %>%
    mutate(significance = case_when(
      log2FoldChange >  LFC_THRESH ~ "Up",
      log2FoldChange < -LFC_THRESH ~ "Down",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > LFC_THRESH) %>%
    dplyr::select(gene_id, gene_name, significance,
                  baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(padj)
  results_list[[comp$name]] <- df
  message(sprintf("  %s: %d DEGs (padj<0.05, |log2FC|>0.585)",
                  comp$name, nrow(df)))
}

# ---- write DESeq2_output.xlsx (Arial 10, matches Supplementary Table 1) ----
write_deg_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  n <- nrow(df)
  addStyle(wb, sheet_name, createStyle(fontName = "Arial", fontSize = 10),
           rows = 1:(n + 1), cols = 1:ncol(df), gridExpand = TRUE, stack = FALSE)
  num_cols <- which(colnames(df) %in% c("baseMean", "log2FoldChange", "lfcSE", "stat"))
  sci_cols <- which(colnames(df) %in% c("pvalue", "padj"))
  if (n > 0) {
    addStyle(wb, sheet_name, createStyle(fontName = "Arial", fontSize = 10, numFmt = "0.0000"),
             rows = 2:(n + 1), cols = num_cols, gridExpand = TRUE, stack = TRUE)
    addStyle(wb, sheet_name, createStyle(fontName = "Arial", fontSize = 10, numFmt = "0.00E+00"),
             rows = 2:(n + 1), cols = sci_cols, gridExpand = TRUE, stack = TRUE)
  }
  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}

wb <- createWorkbook()
for (comp in comparisons) write_deg_sheet(wb, comp$name, results_list[[comp$name]])
saveWorkbook(wb, file.path(OUT_DIR, "DESeq2_output.xlsx"), overwrite = TRUE)
message("Saved: DESeq2_output.xlsx")

# ────────────────────────────────────────────────────────────
# Global VST (13 groups, 39 samples) — input for 02_wgcna.R
# ────────────────────────────────────────────────────────────
dds_vst <- DESeqDataSetFromMatrix(countData = count_vst_flt,
                                  colData   = meta_vst,
                                  design    = ~ group)
dds_vst    <- estimateSizeFactors(dds_vst)
dds_vst    <- estimateDispersions(dds_vst)
vst_obj    <- vst(dds_vst, blind = FALSE)
vst_matrix <- assay(vst_obj)
message("Global VST complete: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples")

vst_df_out <- as.data.frame(vst_matrix) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_anno, by = "gene_id") %>%
  dplyr::select(gene_id, gene_name, everything())
write.xlsx(vst_df_out, file.path(OUT_DIR, "vst_matrix.xlsx"), overwrite = TRUE)
message("Saved: vst_matrix.xlsx")
