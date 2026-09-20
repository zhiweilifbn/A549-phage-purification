# ============================================================
# Core analysis: WGCNA co-expression network + hub gene degree
# ------------------------------------------------------------
# Produces:
#   ../output/WGCNA_output.xlsx  (one sheet per module: rank, gene_id,
#                                 gene_name, weighted_degree, is_hub;
#                                 matches Supplementary Table 6)
#
# Requires ../output/vst_matrix.xlsx, produced by 01_deseq2.R — run
# that script first.
# ============================================================

# ---- paths ----
IN_DIR  <- "../input"
OUT_DIR <- "../output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VST_FILE  <- file.path(OUT_DIR, "vst_matrix.xlsx")
META_FILE <- file.path(IN_DIR, "metadata.xlsx")

# ---- packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(WGCNA)
})
options(stringsAsFactors = FALSE)
set.seed(42)

# ---- group definition: Ctrl + 6 phage preparations + AAV (8 groups, 24 samples) ----
GROUP_ORDER_8 <- c("ctrl",
                   "cru_om", "pur_om",
                   "cru_pf", "pur_pf",
                   "cru_px", "pur_px",
                   "aav")

# ────────────────────────────────────────────────────────────
# Read VST matrix + metadata, select the 8-group / 24-sample subset
# ────────────────────────────────────────────────────────────
vst_raw <- read_excel(VST_FILE)
stopifnot(all(c("gene_id", "gene_name") %in% colnames(vst_raw)))
gene_anno <- vst_raw %>% dplyr::select(gene_id, gene_name) %>% distinct()

vst_matrix <- vst_raw %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

meta_full <- read_excel(META_FILE)
meta8 <- meta_full %>%
  filter(group %in% GROUP_ORDER_8) %>%
  mutate(group = factor(group, levels = GROUP_ORDER_8)) %>%
  arrange(group) %>%
  as.data.frame()
rownames(meta8) <- meta8$sample

vst_matrix <- vst_matrix[, meta8$sample, drop = FALSE]
message("VST matrix: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples (expect 24)")
stopifnot(ncol(vst_matrix) == 24)

# ────────────────────────────────────────────────────────────
# Gene filter for WGCNA input (computed on the 24-sample subset)
# ────────────────────────────────────────────────────────────
vst_medians <- apply(vst_matrix, 1, median)
vst_expr    <- vst_matrix[vst_medians >= 1, ]
gene_mad    <- apply(vst_expr, 1, mad)
n_keep      <- ceiling(length(gene_mad) * 0.50)
top_genes   <- names(sort(gene_mad, decreasing = TRUE))[1:n_keep]
vst_wgcna   <- vst_expr[top_genes, ]
expr_for_wgcna <- t(vst_wgcna)
message("WGCNA input: ", nrow(expr_for_wgcna), " samples x ", ncol(expr_for_wgcna), " genes")

gsg <- goodSamplesGenes(expr_for_wgcna, verbose = 0)
if (!gsg$allOK) {
  expr_for_wgcna <- expr_for_wgcna[gsg$goodSamples, gsg$goodGenes]
  message("After goodSamplesGenes: ", nrow(expr_for_wgcna), " x ", ncol(expr_for_wgcna))
}

# ────────────────────────────────────────────────────────────
# Soft-thresholding power selection (justifies power = 12)
# ────────────────────────────────────────────────────────────
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(expr_for_wgcna,
                         powerVector = powers,
                         networkType = "signed",
                         RsquaredCut = 0.85,
                         verbose = 0)
sft_df <- as.data.frame(sft$fitIndices)
colnames(sft_df) <- c("Power", "SFT_Rsq", "slope", "truncated_Rsq",
                      "mean_k", "median_k", "max_k")
sft_df$signed_Rsq <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]

pick_beta_smart <- function(df, rsq_cut = 0.85, mk_lo = 50, mk_hi = 200) {
  both <- df$signed_Rsq >= rsq_cut & df$mean_k >= mk_lo & df$mean_k <= mk_hi
  if (any(both)) return(min(df$Power[both]))
  both80 <- df$signed_Rsq >= 0.80 & df$mean_k >= mk_lo & df$mean_k <= mk_hi
  if (any(both80)) return(min(df$Power[both80]))
  if (any(df$mean_k >= mk_lo & df$mean_k <= mk_hi)) return(min(df$Power[df$mean_k >= mk_lo & df$mean_k <= mk_hi]))
  NA
}
SELECTED_BETA <- pick_beta_smart(sft_df)
if (is.na(SELECTED_BETA)) SELECTED_BETA <- 12
message("Selected soft-thresholding power (beta) = ", SELECTED_BETA)

# ────────────────────────────────────────────────────────────
# Network construction
# ────────────────────────────────────────────────────────────
TOM_BASE <- file.path(OUT_DIR, "TOM_tmp")
net <- blockwiseModules(
  expr_for_wgcna,
  power = SELECTED_BETA,
  networkType = "signed", TOMType = "signed", corType = "pearson",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = TOM_BASE,
  verbose = 1,
  maxBlockSize = 30000
)
message("WGCNA network construction complete")

mod_colors <- labels2colors(net$colors)
names(mod_colors) <- colnames(expr_for_wgcna)
size_tab <- sort(table(mod_colors), decreasing = TRUE)
size_tab <- size_tab[names(size_tab) != "grey"]
color_to_name <- setNames(paste0("M", sprintf("%02d", seq_along(size_tab))),
                          names(size_tab))
color_to_name["grey"] <- "grey"

MOD_NAME <- color_to_name[mod_colors]
names(MOD_NAME) <- names(mod_colors)
ALL_MODULES <- sort(setdiff(unique(MOD_NAME), "grey"))
message("Modules found (excluding grey): ", length(ALL_MODULES))

membership_df <- data.frame(
  gene_id = names(MOD_NAME),
  mod_name = unname(MOD_NAME),
  stringsAsFactors = FALSE
)

# ────────────────────────────────────────────────────────────
# Load TOM, compute weighted degree per gene, define hub genes
# (top 10% weighted degree within each module)
# ────────────────────────────────────────────────────────────
load(paste0(TOM_BASE, "-block.1.RData"))
TOM_MAT <- as.matrix(TOM)
rownames(TOM_MAT) <- colnames(expr_for_wgcna)
colnames(TOM_MAT) <- colnames(expr_for_wgcna)
rm(TOM)

all_per_mod <- list()
for (mn in ALL_MODULES) {
  g_ids <- membership_df$gene_id[membership_df$mod_name == mn]
  if (length(g_ids) < 10) {
    message("  skipping ", mn, ": fewer than 10 genes")
    next
  }
  tom_sub <- TOM_MAT[g_ids, g_ids]
  diag(tom_sub) <- 0
  wd <- rowSums(tom_sub)
  n_hub <- max(1, ceiling(length(g_ids) * 0.10))
  hub_ids <- names(sort(wd, decreasing = TRUE))[1:n_hub]
  gn <- gene_anno$gene_name[match(g_ids, gene_anno$gene_id)]
  full_df <- data.frame(
    gene_id = g_ids, gene_name = gn,
    weighted_degree = round(wd, 4),
    is_hub = g_ids %in% hub_ids,
    stringsAsFactors = FALSE
  )
  full_df <- full_df[order(full_df$weighted_degree, decreasing = TRUE), ]
  full_df$rank <- seq_len(nrow(full_df))
  full_df <- full_df[, c("rank", "gene_id", "gene_name", "weighted_degree", "is_hub")]
  rownames(full_df) <- NULL
  all_per_mod[[mn]] <- full_df
}

# ---- write WGCNA_output.xlsx (Arial 10, matches Supplementary Table 6) ----
wb <- createWorkbook()
for (mn in names(all_per_mod)) {
  df <- all_per_mod[[mn]]
  addWorksheet(wb, mn)
  writeData(wb, mn, df)
  addStyle(wb, mn, createStyle(fontName = "Arial", fontSize = 10),
           rows = 1:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE, stack = FALSE)
  addStyle(wb, mn, createStyle(fontName = "Arial", fontSize = 10, numFmt = "0.0000"),
           rows = 2:(nrow(df) + 1), cols = which(colnames(df) == "weighted_degree"),
           gridExpand = TRUE, stack = TRUE)
  setColWidths(wb, mn, cols = 1:ncol(df), widths = "auto")
  freezePane(wb, mn, firstRow = TRUE)
}
saveWorkbook(wb, file.path(OUT_DIR, "WGCNA_output.xlsx"), overwrite = TRUE)
message("Saved: WGCNA_output.xlsx")

# ---- clean up intermediate TOM file ----
tom_file <- paste0(TOM_BASE, "-block.1.RData")
if (file.exists(tom_file)) file.remove(tom_file)
