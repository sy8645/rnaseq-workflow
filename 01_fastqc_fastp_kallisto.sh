# Bulk RNA-seq analysis workflow
## step 0. FastQC, fastp filtering, and clean-read FastQC.
## step 1. Kallisto transcript-level quantification.
## step 2. Tximport transcript-to-gene summarization in R.

SAMPLE_TABLE=$1
REFMRNA=$2
REFGENE=$3
OUT=$4
THREADS=${5:-8}

mkdir -p "$OUT"/{fastqc_raw,clean,fastqc_clean,kallisto,reference}

#--------------------- Read quality control ---------------------#
# 1) FastQC for raw reads
# 2) Trim adapters and low-quality reads with fastp defaults
# 3) FastQC for clean reads: verify quality after filtering
tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample condition r1 r2; do
  mkdir -p "$OUT/fastqc_raw/$sample" "$OUT/clean/$sample" "$OUT/fastqc_clean/$sample"
  fastqc -t 4 --extract -o "$OUT/fastqc_raw/$sample" "$r1" "$r2"

  c1="$OUT/clean/$sample/${sample}_R1.clean.fq.gz"
  c2="$OUT/clean/$sample/${sample}_R2.clean.fq.gz"
  fastp -i "$r1" -I "$r2" -o "$c1" -O "$c2" -w "$THREADS" \
    -j "$OUT/clean/$sample/${sample}.fastp.json" \
    -h "$OUT/clean/$sample/${sample}.fastp.html"

  fastqc -t 4 --extract -o "$OUT/fastqc_clean/$sample" "$c1" "$c2"
done

#--------------------- Kallisto ---------------------#
# 1) Build the hg38 refMrna transcriptome index
kallisto index -i "$OUT/reference/refMrna.kallisto.idx" "$REFMRNA"

# 2) Quantify clean paired-end reads
printf 'sample\tcondition\tabundance_tsv\n' > "$OUT/samples.kallisto.tsv"
tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample condition r1 r2; do
  c1="$OUT/clean/$sample/${sample}_R1.clean.fq.gz"
  c2="$OUT/clean/$sample/${sample}_R2.clean.fq.gz"
  kallisto quant -i "$OUT/reference/refMrna.kallisto.idx" \
    -o "$OUT/kallisto/$sample" -t "$THREADS" "$c1" "$c2"
  printf '%s\t%s\tkallisto/%s/abundance.tsv\n' "$sample" "$condition" "$sample" \
    >> "$OUT/samples.kallisto.tsv"
done

#--------------------- Tximport ---------------------#
# 1) Make the RefSeq transcript-to-gene table
gzip -dc "$REFGENE" | awk -F '\t' 'BEGIN{OFS="\t"} NF>=13 && $2!="" && $13!="" {print $2,$13}' \
  | LC_ALL=C sort -u > "$OUT/reference/tx2gene.tsv"

# 2) Import Kallisto and summarize gene-level counts and TPM (R)
# Rscript runs the R code between the two RSCRIPT markers.
Rscript - "$OUT" <<'RSCRIPT'
library(tximport)
out <- normalizePath(commandArgs(trailingOnly = TRUE)[1])
samples <- read.delim(file.path(out, "samples.kallisto.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(samples) == 12, !anyDuplicated(samples$sample),
          sum(samples$condition == "control") == 6, sum(samples$condition == "test") == 6)
files <- setNames(file.path(out, samples$abundance_tsv), samples$sample)
tx2gene <- read.delim(file.path(out, "reference/tx2gene.tsv"), header = FALSE)
stopifnot(!anyDuplicated(tx2gene[[1]]), all(file.exists(files)))
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                countsFromAbundance = "no", dropInfReps = TRUE)
stopifnot(identical(colnames(txi$counts), samples$sample))
dir.create(file.path(out, "counts"), showWarnings = FALSE)
write_matrix <- function(x, name) write.table(
  data.frame(gene_symbol = rownames(x), x, check.names = FALSE),
  file.path(out, "counts", name), sep = "\t", quote = FALSE, row.names = FALSE
)
write_matrix(txi$counts, "gene_estimated_counts_all.tsv")
write_matrix(txi$abundance, "gene_TPM_all.tsv")
write_matrix(txi$length, "gene_effective_length_all.tsv")
saveRDS(txi, file.path(out, "counts/tximport_gene_level.rds"))
RSCRIPT

