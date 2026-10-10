# Bulk RNA-seq DEG analysis workflow
## step 3. DESeq2: unnormalized gene estimated counts.
## step 4. EnhancedVolcano: log2FC and nominal P values.


library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
count_file <- normalizePath(args[1], mustWork = TRUE)
sample_file <- normalizePath(args[2], mustWork = TRUE)
out <- args[3]
label_file <- args[4]
dir.create(file.path(out, "counts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "differential"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, file) write.table(
  x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write_matrix <- function(x, file) write_tsv(
  data.frame(gene_symbol = rownames(x), x, check.names = FALSE), file
)

#--------------------- DESeq2 prep counts ---------------------#
# Read unnormalized gene estimated counts exported by 01
count <- as.matrix(read.delim(count_file, row.names = 1, check.names = FALSE))
sample_info <- read.delim(sample_file, stringsAsFactors = FALSE)
stopifnot(nrow(sample_info) == 12, !anyDuplicated(sample_info$sample),
          sum(sample_info$condition == "control") == 6,
          sum(sample_info$condition == "test") == 6,
          !anyDuplicated(colnames(count)), setequal(colnames(count), sample_info$sample))
# Set condition and align count columns to the sample table
sample_info$condition <- factor(sample_info$condition, levels = c("control", "test"))
rownames(sample_info) <- sample_info$sample
count <- count[, rownames(sample_info), drop = FALSE]
stopifnot(identical(colnames(count), rownames(sample_info)))
# Feng count preparation: round to 2 decimals, then integers; keep mean > 10
stopifnot(is.numeric(count), all(is.finite(count)), all(count >= 0))
count <- round(round(count, 2))
storage.mode(count) <- "integer"
count <- count[rowMeans(count) > 10, , drop = FALSE]
write_matrix(count, file.path(out, "counts/gene_DESeq2_counts_mean_gt_10.tsv"))

#--------------------- DEG between groups ---------------------#
dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = sample_info,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "test", "control"))

#--------------------- Export DEG results ---------------------#
deg <- data.frame(gene_symbol = rownames(res), as.data.frame(res), check.names = FALSE)
deg <- deg[order(deg$pvalue, na.last = TRUE), ]
write_tsv(deg, file.path(out, "differential/DESeq2_unshrunk_all_genes_test_vs_control.tsv"))
write_tsv(
  deg[!is.na(deg$pvalue) & deg$pvalue < 0.05, ],
  file.path(out, "differential/DEG_unshrunk_pvalue_lt_0.05.tsv")
)
write_tsv(
  deg[!is.na(deg$padj) & deg$padj < 0.05, ],
  file.path(out, "differential/DEG_unshrunk_padj_lt_0.05.tsv")
)
write_tsv(
  deg[!is.na(deg$pvalue) & deg$pvalue < 0.05 & abs(deg$log2FoldChange) >= 0.5, ],
  file.path(out, "differential/DEG_unshrunk_pvalue_lt_0.05_abs_log2FC_ge_0.5.tsv")
)
write_tsv(
  deg[!is.na(deg$padj) & deg$padj < 0.05 & abs(deg$log2FoldChange) >= 0.5, ],
  file.path(out, "differential/DEG_unshrunk_padj_lt_0.05_abs_log2FC_ge_0.5.tsv")
)
write_matrix(counts(dds, normalized = TRUE), file.path(out, "counts/DESeq2_normalized_counts.tsv"))

# Plot dispersion and MA diagnostics
pdf(file.path(out, "figures/DESeq2_dispersion.pdf"), width = 6, height = 6)
plotDispEsts(dds)
dev.off()
pdf(file.path(out, "figures/DESeq2_MA.pdf"), width = 6, height = 6)
plotMA(res, alpha = 0.05, ylim = c(-6, 6))
dev.off()

#--------------------- EnhancedVolcano ---------------------#
label_font_size <- 3.6
plot_data <- deg[is.finite(deg$log2FoldChange) & !is.na(deg$pvalue), ]
plot_data$pvalue_plot <- ifelse(plot_data$pvalue == 0, 1e-50, plot_data$pvalue)
rownames(plot_data) <- make.unique(plot_data$gene_symbol)

p <- EnhancedVolcano(
  plot_data,
  lab = rep("", nrow(plot_data)), selectLab = character(0),
  x = "log2FoldChange", y = "pvalue_plot",
  pCutoff = 0.05, FCcutoff = 0.5,
  xlim = c(-2, 2), ylim = c(0, 15),
  col = c("#333333", "#55A868", "#4C72B0", "#E64B35"),
  colAlpha = 0.78, pointSize = 1.6,
  cutoffLineType = "dashed", cutoffLineCol = "#666666",
  title = "CRISPRi of the deletion", subtitle = NULL, caption = NULL,
  xlab = expression(log[2] * "FC"), ylab = expression(-log[10] * italic(P)),
  legendPosition = "none", drawConnectors = FALSE,
  gridlines.major = FALSE, gridlines.minor = FALSE, border = "partial"
) + theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

#--------------------- Volcano plot ---------------------#
label_genes <- unique(scan(label_file, what = "character", quiet = TRUE))
label_data <- plot_data[match(label_genes, plot_data$gene_symbol), , drop = FALSE]
label_data <- label_data[!is.na(label_data$gene_symbol), , drop = FALSE]
label_data$y <- -log10(label_data$pvalue_plot)

for (left_side in c(TRUE, FALSE)) {
  z <- label_data[(label_data$log2FoldChange < 0) == left_side, , drop = FALSE]
  if (nrow(z) == 0) next
  target_x <- if (left_side) -0.6 else 0.6
  p <- p + geom_text_repel(
    data = z, aes(log2FoldChange, y, label = gene_symbol),
    nudge_x = target_x - z$log2FoldChange,
    direction = "y", hjust = if (left_side) 1 else 0,
    size = label_font_size, color = "black", segment.color = "black",
    min.segment.length = 0, box.padding = 0.36, point.padding = 0.18,
    segment.curvature = 0, arrow = NULL, max.overlaps = Inf, seed = 2026
  )
}

prefix <- file.path(out, "figures/volcano_unshrunk_P0.05_LFC0.5_EnhancedVolcano")
pdf(paste0(prefix, "_AI_editable.pdf"), width = 6, height = 5.8, useDingbats = FALSE)
print(p)
dev.off()
ggsave(paste0(prefix, ".png"), p, width = 6, height = 5.8, dpi = 600)
write_tsv(
  label_data[, c("gene_symbol", "log2FoldChange", "pvalue", "padj")],
  paste0(prefix, "_labeled_genes.tsv")
)
