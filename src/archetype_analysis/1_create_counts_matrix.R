#!/usr/bin/env Rscript

# tcga_counts_normalize.R
# Build a counts matrix from STAR ReadsPerGene.out.tab files and normalize with DESeq2.
# Assumes STAR output format: first 4 header rows, columns: geneID, Unstranded, Stranded1, Stranded2.
# Defaults: use column 2 (Unstranded), skip first 4 rows.

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required. Please install it (e.g., BiocManager::install('DESeq2')).")
  }
  library(DESeq2)
})


args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}

counts_dir  <- normalizePath(get_arg("--counts_dir", "../data/counts"), mustWork = FALSE)
pattern     <- get_arg("--pattern", NULL)             # e.g. "ReadsPerGene.out.tab$"
out_dir     <- normalizePath(get_arg("--out_dir", "../data"), mustWork = FALSE)
out_prefix  <- get_arg("--out_prefix", "tcga")
skip_header <- as.integer(get_arg("--skip_header", "4"))
column_idx  <- as.integer(get_arg("--column", "2"))   # 2 = Unstranded
min_nonzero <- as.integer(get_arg("--min_nonzero", "1"))

if (!dir.exists(counts_dir)) stop("counts_dir does not exist: ", counts_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Parameters:")
message("  counts_dir  : ", counts_dir)
message("  pattern     : ", ifelse(is.null(pattern), "(none)", pattern))
message("  out_dir     : ", out_dir)
message("  out_prefix  : ", out_prefix)
message("  skip_header : ", skip_header)
message("  column      : ", column_idx)
message("  min_nonzero : ", min_nonzero)


all_files <- list.files(counts_dir, full.names = TRUE)
if (!is.null(pattern)) {
  all_files <- all_files[grepl(pattern, basename(all_files))]
}
if (length(all_files) == 0) stop("No count files found in: ", counts_dir)


read_star_counts <- function(path, skip_header = 4, col = 2) {
  df <- utils::read.table(path, header = FALSE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  if (nrow(df) <= skip_header) stop("File has fewer rows than 'skip_header': ", path)
  df <- df[(skip_header + 1):nrow(df), , drop = FALSE]
  if (ncol(df) < col) stop("Requested column ", col, " not found in file: ", path)
  
  gene_ids <- df[[1]]
  counts   <- df[[col]]
  
  
  counts <- suppressWarnings(as.numeric(counts))
  if (anyNA(counts)) stop("Non-numeric counts encountered in: ", path)
  
  names(counts) <- gene_ids
  counts
}


sample_name_from_file <- function(fname) sub("_.*", "", fname)

samples <- basename(all_files)
sample_names <- vapply(samples, sample_name_from_file, character(1))

message("Reading ", length(all_files), " files...")
counts_list <- vector("list", length(all_files))
genes_ref <- NULL

for (i in seq_along(all_files)) {
  vec <- read_star_counts(all_files[i], skip_header = skip_header, col = column_idx)
  
  if (is.null(genes_ref)) {
    genes_ref <- names(vec)
  } else {
    # Align genes if needed
    if (!identical(genes_ref, names(vec))) {
      # fallback to intersection with order of reference
      common <- intersect(genes_ref, names(vec))
      if (length(common) == 0) stop("No overlapping genes between files. Check inputs.")
      vec <- vec[common]
      if (length(common) < length(genes_ref)) {
        genes_ref <- common
      }
    }
  }
  counts_list[[i]] <- vec
}


common_genes <- Reduce(intersect, lapply(counts_list, names))
if (length(common_genes) == 0) stop("No common genes across files.")
counts_list <- lapply(counts_list, function(v) v[common_genes])

cts <- do.call(cbind, counts_list)
colnames(cts) <- sample_names
rownames(cts) <- common_genes


cts <- cts[order(rownames(cts)), , drop = FALSE]

unnorm_path <- file.path(out_dir, paste0(out_prefix, "_unnormalized_counts.csv"))
utils::write.csv(cts, unnorm_path, quote = TRUE)
message("Wrote unnormalized counts: ", unnorm_path)

# ---- Build DESeq2 object & normalize (size factors only) ---------------------
coldata <- DataFrame(dummy = factor(rep("Standard", ncol(cts))))
rownames(coldata) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = round(cts, 0), colData = coldata, design = ~ 1)

# Estimate size factors 
dds <- estimateSizeFactors(dds)

size_factors <- sizeFactors(dds)
sf_path <- file.path(dirname(out_dir), "supplementary", "sizeFactors.csv")
dir.create(dirname(sf_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(sample = names(size_factors), sizeFactor = size_factors),
                 sf_path, row.names = FALSE)
message("Wrote size factors: ", sf_path)

counts_norm <- counts(dds, normalized = TRUE)
counts_log  <- log2(1 + counts_norm)


if (!is.na(min_nonzero) && min_nonzero > 0) {
  keep <- rowSums(counts_norm > 0) >= min_nonzero
  counts_norm <- counts_norm[keep, , drop = FALSE]
  counts_log  <- counts_log[keep, , drop = FALSE]
  message("Kept ", sum(keep), " genes with at least ", min_nonzero, " nonzero counts.")
}

#  Write outputs 
norm_path <- file.path(out_dir, paste0(out_prefix, "_normalized_counts.csv"))
log_path  <- file.path(out_dir, paste0(out_prefix, "_log_normalized_counts.csv"))
samples_path <- file.path(out_dir, paste0(out_prefix, "_samples.csv"))

utils::write.csv(counts_norm, norm_path, quote = TRUE)
utils::write.csv(counts_log,  log_path,  quote = TRUE)
utils::write.csv(data.frame(file = basename(all_files), sample = sample_names),
                 samples_path, row.names = FALSE)

message("Wrote normalized counts: ", norm_path)
message("Wrote log2(1+norm) counts: ", log_path)
message("Wrote sample mapping: ", samples_path)


