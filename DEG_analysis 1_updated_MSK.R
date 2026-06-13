# =============================================================================
# SETUP
# =============================================================================

if (!("BiocManager" %in% installed.packages())) install.packages("BiocManager")
pkgs <- c("DESeq2", "ggplot2", "RColorBrewer", "pheatmap", "tidyverse",
          "GEOquery", "data.table", "AnnotationDbi", "org.Hs.eg.db",
          "apeglm", "EnhancedVolcano", "ggrepel", "ggVennDiagram",
          "clusterProfiler", "enrichplot", "RCy3", "msigdbr",
          "readr", "limma", "patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

library(DESeq2); library(ggplot2); library(RColorBrewer); library(pheatmap)
library(tidyverse); library(GEOquery); library(data.table)
library(AnnotationDbi); library(org.Hs.eg.db); library(apeglm)
library(EnhancedVolcano); library(ggrepel); library(ggVennDiagram)
library(clusterProfiler); library(enrichplot); library(RCy3)
library(msigdbr); library(readr); library(limma); library(patchwork)

# =============================================================================
# PARAMETERS
# =============================================================================

LOG2FC_CUTOFF  <- 1
PVAL_CUTOFF    <- 0.05
TOP_N_HEATMAP  <- 100
TOP_N_VOLCANO  <- 20
TOP_N_PATHWAYS <- 20

# =============================================================================
# DATA LOADING
# =============================================================================

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE151243",
              "file=GSE151243_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "&")
tbl <- as.matrix(data.table::fread(path, header = TRUE, colClasses = "integer"),
                 rownames = 1)

gse      <- getGEO("GSE151243", GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(phenoData(gse[[1]]))

metadata.filt <- metadata[, c("geo_accession", "race:ch1",
                              "patient_id:ch1", "characteristics_ch1.1")]
colnames(metadata.filt) <- c("geo_accession", "race", "patient_id", "tissue_type")
metadata.filt$tissue_type <- gsub("tissue_type: ", "", metadata.filt$tissue_type)
metadata.filt$tissue_type <- factor(metadata.filt$tissue_type,
                                    levels = c("Lesion", "Perilesion"))
metadata.filt$race       <- factor(metadata.filt$race,
                                   levels = c("W", "B", "W_A", "O"))
metadata.filt$patient_id <- factor(metadata.filt$patient_id,
                                   levels = c("HS_10","HS_11","HS_12","HS_13",
                                              "HS_14","HS_15","HS_16","HS_17",
                                              "HS_18","HS_19","HS_01","HS_21",
                                              "HS_22","HS_02","HS_04","HS_05",
                                              "HS_06","HS_07","HS_08","HS_09"))

# Annotation
data_all <- as.data.frame(tbl)
anno <- AnnotationDbi::select(org.Hs.eg.db,
                              keys    = rownames(data_all),
                              columns = c("ENSEMBL", "SYMBOL"),
                              keytype = "ENTREZID")
anno           <- anno[!duplicated(anno$ENTREZID), ]
rownames(anno) <- anno$ENTREZID

data_all <- data_all[rownames(data_all) %in% anno$ENTREZID, ]
anno      <- anno[rownames(data_all), ]

# Per-race count matrices and metadata
metadata.b <- metadata.filt[metadata.filt$race == "B", ]
metadata.w <- metadata.filt[metadata.filt$race == "W", ]

data.b <- data_all[, colnames(data_all) %in% metadata.b$geo_accession, drop = FALSE]
data.w <- data_all[, colnames(data_all) %in% metadata.w$geo_accession, drop = FALSE]

stopifnot(all(colnames(data.b) == metadata.b$geo_accession),
          all(colnames(data.w) == metadata.w$geo_accession))

# Low-count filter
data.b.filt <- data.b[rowSums(data.b) > 10, ]
data.w.filt <- data.w[rowSums(data.w) > 10, ]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Run DESeq2, annotate results, add sig labels and display label
run_deseq <- function(counts, meta, anno) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = meta,
                                design    = ~ patient_id + tissue_type)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("tissue_type", "Lesion", "Perilesion"))
  res <- res[order(res$padj), ]
  summary(res)
  
  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::mutate(sig = dplyr::case_when(
      padj < PVAL_CUTOFF & log2FoldChange > LOG2FC_CUTOFF ~ "Up in Lesion",
      padj < PVAL_CUTOFF & log2FoldChange < -LOG2FC_CUTOFF ~ "Up in Perilesion",
      TRUE ~ "Not Significant"
    ))
  
  res_df_id <- dplyr::left_join(res_df, anno, by = c("gene" = "ENTREZID")) |>
    dplyr::filter(!is.na(SYMBOL)  & SYMBOL  != "",
                  !is.na(ENSEMBL) & ENSEMBL != "") |>
    dplyr::mutate(label = dplyr::if_else(!is.na(SYMBOL) & SYMBOL != "",
                                         SYMBOL, gene))
  
  vsd     <- vst(dds, blind = TRUE)   # for QC plots
  vsd_res <- vst(dds, blind = FALSE)  # for results plots
  
  list(dds = dds, res = res, res_df_id = res_df_id,
       vsd = vsd, vsd_res = vsd_res)
}

#' Volcano plot with EnhancedVolcano
plot_volcano <- function(res_df_id, title) {
  fc_max <- max(abs(res_df_id$log2FoldChange), na.rm = TRUE)
  p_max  <- max(-log10(res_df_id$padj), na.rm = TRUE)
  n_degs <- sum(res_df_id$sig != "Not Significant")
  EnhancedVolcano(res_df_id,
                  lab      = res_df_id$label,
                  x        = "log2FoldChange",
                  y        = "padj",
                  pCutoff  = PVAL_CUTOFF,
                  FCcutoff = LOG2FC_CUTOFF,
                  labSize  = 3,
                  xlim     = c(-fc_max, fc_max),
                  ylim     = c(0, p_max),
                  title    = paste0(title, " (", n_degs, " DEGs)"),
                  subtitle = "Lesion vs Perilesion",
                  col      = c("grey70", "steelblue", "firebrick", "darkred"),
                  colAlpha = 0.6)
}

#' Heatmap of top DE genes
plot_heatmap <- function(res_df_id, vsd_res, title, n = TOP_N_HEATMAP) {
  top_de <- res_df_id |>
    dplyr::filter(padj < PVAL_CUTOFF) |>
    dplyr::slice_min(padj, n = n)
  
  mat <- assay(vsd_res)[top_de$gene, ]
  rownames(mat) <- top_de$SYMBOL
  mat <- mat - rowMeans(mat)
  
  anno_col <- data.frame(tissue_type = vsd_res$tissue_type,
                         row.names   = colnames(vsd_res))
  my_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  pheatmap(mat,
           color           = my_colors,
           annotation_col  = anno_col,
           show_rownames   = TRUE,
           show_colnames   = FALSE,
           scale           = "none",
           fontsize_row    = 6,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           main = paste0("Top ", n, " DE genes (VST, row-centred) - ", title))
}

#' PCA plot
plot_pca <- function(vsd, intgroup, title, color_var, shape_var = NULL,
                     color_vals = NULL, shape_vals = NULL) {
  pca  <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  pct  <- round(100 * attr(pca, "percentVar"))
  p <- ggplot(pca, aes(x = PC1, y = PC2,
                       color  = .data[[color_var]],
                       shape  = if (!is.null(shape_var)) .data[[shape_var]] else NULL,
                       label  = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    xlab(paste0("PC1: ", pct[1], "% variance")) +
    ylab(paste0("PC2: ", pct[2], "% variance")) +
    ggtitle(title) +
    theme_bw()
  if (!is.null(color_vals)) p <- p + scale_color_manual(values = color_vals)
  if (!is.null(shape_vals)) p <- p + scale_shape_manual(values = shape_vals)
  p
}

#' Run WikiPathways enrichment and produce plots
run_pathway_enrichment <- function(degs, genesets, title) {
  res_wp <- clusterProfiler::enricher(degs$gene, TERM2GENE = genesets,
                                      pAdjustMethod = "fdr",
                                      pvalueCutoff  = PVAL_CUTOFF,
                                      minGSSize = 5, maxGSSize = 400)
  res_wp_df  <- as.data.frame(res_wp)
  res_wp_sim <- enrichplot::pairwise_termsim(res_wp)
  
  print(treeplot(res_wp_sim, showCategory = 30, label_format = 30) +
          ggtitle(paste("Treeplot -", title)))
  print(dotplot(res_wp, showCategory = TOP_N_PATHWAYS, label_format = 40) +
          ggtitle(paste("WikiPathways Enrichment -", title)) +
          theme(axis.text.y = element_text(size = 8)))
  print(barplot(res_wp, showCategory = TOP_N_PATHWAYS) +
          ggtitle(paste("WikiPathways Enrichment -", title)) +
          theme(axis.text.y = element_text(size = 8)))
  
  list(res_wp = res_wp, res_wp_df = res_wp_df, res_wp_sim = res_wp_sim)
}

#' Build cyemapplot input data frame
make_mapped_data <- function(res_df_id) {
  res_df_id |>
    dplyr::select(gene, log2FoldChange, pvalue, padj) |>
    dplyr::rename(ID = gene, log2FC = log2FoldChange) |>
    dplyr::filter(!is.na(log2FC), !is.na(pvalue), !is.na(padj)) |>
    dplyr::distinct(ID, .keep_all = TRUE) |>
    dplyr::mutate(ranking = log2FC * -log10(pvalue))
}

# =============================================================================
# STEP 1: DESEQ2
# =============================================================================

cat("Running DESeq2 for Black patients...\n")
res_b <- run_deseq(data.b.filt, metadata.b, anno)

cat("Running DESeq2 for White patients...\n")
res_w <- run_deseq(data.w.filt, metadata.w, anno)

# Convenient aliases
res_df_id.b <- res_b$res_df_id
res_df_id.w <- res_w$res_df_id

# =============================================================================
# PLOTS: QC
# =============================================================================

# PCA per cohort (blind VST)
plot_pca(res_b$vsd, c("tissue_type", "patient_id"),
         "PCA - Black patients (blind VST)", "tissue_type")
plot_pca(res_w$vsd, c("tissue_type", "patient_id"),
         "PCA - White patients (blind VST)", "tissue_type")

# Dispersion estimates
plotDispEsts(res_b$dds, main = "Dispersion - Black patients")
plotDispEsts(res_w$dds, main = "Dispersion - White patients")

# MA plots with shrinkage
for (cohort in list(list(dds = res_b$dds, label = "Black"),
                    list(dds = res_w$dds, label = "White"))) {
  res_shrunk <- lfcShrink(cohort$dds,
                          coef = "tissue_type_Perilesion_vs_Lesion",
                          type = "apeglm")
  DESeq2::plotMA(res_shrunk,
                 main = paste("MA plot (apeglm) -", cohort$label,
                              "- Perilesion vs Lesion"),
                 ylim = c(-5, 5))
}

# =============================================================================
# PLOTS: RESULTS
# =============================================================================

# Volcano plots
plot_volcano(res_df_id.b, "Volcanoplot - Black patients")
plot_volcano(res_df_id.w, "Volcanoplot - White patients")

# Heatmaps
plot_heatmap(res_df_id.b, res_b$vsd_res, "Black patients")
plot_heatmap(res_df_id.w, res_w$vsd_res, "White patients")

# PCA on results VST
plot_pca(res_b$vsd_res, c("tissue_type", "patient_id"),
         "PCA after VST - Black patients", "tissue_type")
plot_pca(res_w$vsd_res, c("tissue_type", "patient_id"),
         "PCA after VST - White patients", "tissue_type")

# =============================================================================
# COMBINED PCA (both cohorts)
# =============================================================================

common_genes     <- intersect(rownames(data.b.filt), rownames(data.w.filt))
data_combined    <- cbind(data.b.filt[common_genes, ], data.w.filt[common_genes, ])
metadata_combined <- rbind(metadata.b, metadata.w)

dds_combined <- DESeqDataSetFromMatrix(countData = data_combined,
                                       colData   = metadata_combined,
                                       design    = ~ race + tissue_type)
dds_combined <- estimateSizeFactors(dds_combined)
vsd_combined <- vst(dds_combined, blind = FALSE)

plot_pca(vsd_combined, c("race", "tissue_type"),
         "PCA - Black & White combined (uncorrected)", "race",
         shape_var   = "tissue_type",
         color_vals  = c("B" = "firebrick", "W" = "steelblue"),
         shape_vals  = c("Lesion" = 16, "Perilesion" = 17))

# Patient-effect correction with limma
assay(vsd_combined) <- removeBatchEffect(assay(vsd_combined),
                                         batch = vsd_combined$patient_id)
plot_pca(vsd_combined, c("race", "tissue_type"),
         "PCA - Black & White combined (patient effect removed)", "race",
         shape_var   = "tissue_type",
         color_vals  = c("B" = "firebrick", "W" = "steelblue"),
         shape_vals  = c("Lesion" = 16, "Perilesion" = 17))

# Scree plot
pca_full    <- prcomp(t(assay(vsd_combined)))
pca_var_pct <- round(pca_full$sdev^2 / sum(pca_full$sdev^2) * 100, 1)
pca_df      <- data.frame(PC       = factor(paste0("PC", 1:20), levels = paste0("PC", 1:20)),
                          Variance = pca_var_pct[1:20])
ggplot(pca_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1), color = "firebrick", linewidth = 0.8) +
  geom_point(color = "firebrick", size = 2) +
  xlab("Principal Component") + ylab("% Variance Explained") +
  ggtitle("Screeplot - Combined PCA") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5))

# =============================================================================
# VENN DIAGRAM + UNIQUE DEG DIRECTIONS
# =============================================================================

up_black   <- res_df_id.b[res_df_id.b$log2FoldChange > 1 & res_df_id.b$padj < 0.05,"SYMBOL"]
down_black   <- res_df_id.b[res_df_id.b$log2FoldChange < -1 & res_df_id.b$padj < 0.05,"SYMBOL"]

up_white  <- res_df_id.w[res_df_id.w$log2FoldChange > 1 & res_df_id.w$padj < 0.05,"SYMBOL"]
down_white   <- res_df_id.w[res_df_id.w$log2FoldChange < -1 & res_df_id.w$padj < 0.05,"SYMBOL"]

venn_sets <- list(
  "Up Black"   = up_black,
  "Down Black" = down_black,
  "Up White"   = up_white,
  "Down White" = down_white
)

venn_plot <- ggVennDiagram(venn_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("DEG overlap: Up/Down in Black vs White patients") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(venn_plot)

# =============================================================================
# STEP 2: PATHWAY ENRICHMENT ANALYSIS (PEA)
# =============================================================================

genesets.wp <- msigdbr(species = "Homo sapiens",
                       category    = "C2",
                       subcategory = "CP:WIKIPATHWAYS") |>
  dplyr::select(gs_name, entrez_gene)

degs.b <- res_df_id.b[abs(res_df_id.b$log2FoldChange) > LOG2FC_CUTOFF &
                        res_df_id.b$padj < PVAL_CUTOFF, ]
degs.w <- res_df_id.w[abs(res_df_id.w$log2FoldChange) > LOG2FC_CUTOFF &
                        res_df_id.w$padj < PVAL_CUTOFF, ]

pea_b <- run_pathway_enrichment(degs.b, genesets.wp, "Black patients")
pea_w <- run_pathway_enrichment(degs.w, genesets.wp, "White patients")

#GO ENRICHMENT
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")

library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
enrichFrame.b <- enrichGO(gene = res_df_id.b$ENSEMBL,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
enrichResult.b <- as.data.frame(enrichFrame.b)
barplot(enrichFrame.b,
        x = "GeneRatio",
        color = "p.adjust",
        title = "Top 20 of GO Enrichment HS_Black",
        showCategory = 20,
        label_format = 80
)

enrichFrame.w <- enrichGO(gene = res_df_id.w$ENSEMBL,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENSEMBL",
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
enrichResult.w <- as.data.frame(enrichFrame.w)
barplot(enrichFrame.w,
        x = "GeneRatio",
        color = "p.adjust",
        title = "Top 20 of GO Enrichment HS_White",
        showCategory = 20,
        label_format = 80
)

# =============================================================================
# CYEMAPPLOT
# =============================================================================

source("cyemapplot.R")
outdir <- "Output-plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cytoscapePing()

cyemapplot(pea_b$res_wp_sim,
           show_category    = nrow(pea_b$res_wp_df),
           min_edge         = 0.2, visualization = "deg",
           degs_data        = make_mapped_data(res_df_id.b),
           ig_layout        = igraph::layout_with_kk,
           layout_scale     = 1000, min_cluster_size = 1,
           plot_clusters    = TRUE, top_clusters     = 3,
           analysis_name    = "HS Black patients - WikiPathways")

cyemapplot(pea_w$res_wp_sim,
           show_category    = nrow(pea_w$res_wp_df),
           min_edge         = 0.2, visualization = "deg",
           degs_data        = make_mapped_data(res_df_id.w),
           ig_layout        = igraph::layout_with_kk,
           layout_scale     = 1000, min_cluster_size = 1,
           plot_clusters    = TRUE, top_clusters     = 3,
           analysis_name    = "HS White patients - WikiPathways")
