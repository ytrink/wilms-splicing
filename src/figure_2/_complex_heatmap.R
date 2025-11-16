
# Figure 2: ComplexHeatmap Visualization
# ------------------------------------------------------------------------------
# Generates complex heatmaps for Wilms tumor samples using multiple clustering
# methods and distance metrics. Produces both PDF and SVG versions for figures.
# ==============================================================================

# for more information on the complex heatmap package - see: https://jokergoo.github.io/ComplexHeatmap-reference/book/



# ----------------------------- Libraries -------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in c("dplyr", "devtools", "ComplexHeatmap", "stringr",
              "matrixStats", "readxl", "svglite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(dplyr)
library(ComplexHeatmap)
library(stringr)
library(matrixStats)
library(readxl)
library(svglite)

# steps: 1) load data
# 2) fix sample names
# 3) process annotation col

# 4) generate heatmaps

# 1) load data
message("Loading expression matrix...")
counts <- read.csv("../data/tcga_counts_all_threeArcs.csv", row.names = 1)

message("Loading selected genes...")
selected_genes_table <- read.csv("../data/selected_genes_annotated.csv")
selected_genes <- selected_genes_table[, 1]
selected_genes_annotations <- data.frame(Gene.Annotation = selected_genes_table[, 4],
                                         row.names = selected_genes)

message("Preparing heatmap matrix...")
heatmap_matrix <- as.matrix(counts[selected_genes, ])
heatmap_matrix_stand <- t(scale(t(heatmap_matrix), center = TRUE, scale = TRUE))


message("Loading deconvolution output...")
deconv <- read.csv("../deconvolution_output/cpm_cell_type_predictions.csv", row.names = 1)
sample_names <- rownames(deconv)

# 2) Fix sample names
sample_names[1:136] <- str_sub(sample_names[1:136], 2)
rownames(deconv) <- sample_names

message("Loading metadata...")
annotationCol <- read.csv("../data/metadataAll.csv")
rownames(annotationCol) <- str_replace_all(annotationCol$`File Name`, "-", "_")

# 3) process annotationCol
annotationCol <- annotationCol %>%
  select(Sample.Type, Stage,
         Histologic.Classification.of.Primary.Tumor.Normals.Labelled,
         Cap.mesenchyme, Endothelium, Loop.of.Henle, Podocyte,
         Proliferating.cap.mesenchyme, Proximal.tubule, UB, Macrophage,
         Fibroblast, Myofibroblast, Renal.Vesicle, S.shaped.body,
         Manual.Classification, dist_from_Arc1_Blast, dist_from_Arc2_Stromal,
         dist_from_Arc3_Normal)

colnames(annotationCol) <- c("Sample Type", "Stage", "Histologic Classification",
                             "Cap mesenchyme", "Endothelium", "Loop of Henle",
                             "Podocyte", "Proliferating cap mesenchyme", "Proximal tubule",
                             "UB", "Macrophage", "Fibroblast", "Myofibroblast",
                             "Renal Vesicle", "S-shaped body", "Manual Classification",
                             "dist_from_Arc1_Blast", "dist_from_Arc2_Stromal", "dist_from_Arc3_Normal")

# Replace normals in stage with “Normal”
idx.normal <- which(annotationCol$`Sample Type` == "Solid Tissue Normal")
annotationCol$Stage[idx.normal] <- "Normal"

colors <- list(
  "Sample Type" = c("Primary Tumor" = "Blue", "Solid Tissue Normal" = "Green",
                    "Recurrent Tumor" = "Red", "Metastatic" = "Purple", "Arc" = "Yellow"),
  "Stage" = c("I" = "Blue", "II" = "Grey", "IIIB" = "Brown", "III" = "Red",
              "IV" = "White", "Normal" = "Green", "Arc" = "Yellow"),
  "Histologic Classification" = c("DAWT" = "Blue", "FHWT" = "Red",
                                  "Arc" = "Yellow", "Solid Tissue Normal" = "White"),
  "Manual Classification" = c("Epithelial_1" = "Blue", "Epithelial_2" = "Red",
                              "Unspecified Continuum" = "Green", "Normal" = "White", "Arc" = "Yellow")
)

colors.rows <- list(
  "Gene.Annotation" = c("Epithelial" = "Blue", "Nephrogenic Stroma" = "Red",
                        "Cap Mesenchyme" = "Green", "Immune" = "Black",
                        "Cell cycle" = "Yellow", "Housekeeping" = "Purple",
                        "Muscle" = "White", "Endothelial" = "Brown")
)

topAnnotation <- HeatmapAnnotation(df = annotationCol, col = colors, na_col = "black")
sideAnnotation <- rowAnnotation(df = selected_genes_annotations, col = colors.rows)

# helper 
save_heatmap <- function(filename, method = "complete", dist = "spearman", format = "pdf") {
  if (format == "pdf") pdf(file = filename, width = 20, height = 12)
  if (format == "svg") svglite(file = filename, width = 20, height = 12)
  
  draw(Heatmap(heatmap_matrix_stand,
               show_row_names = TRUE,
               clustering_distance_rows = dist,
               clustering_distance_columns = dist,
               clustering_method_rows = method,
               clustering_method_columns = method,
               top_annotation = topAnnotation,
               right_annotation = sideAnnotation,
               name = "Expression (z-score)",
               column_title = paste(method, "/", dist),
               na_col = "black"))
  
  dev.off()
  message("Saved: ", filename)
}


# 4) generate heatmaps 
dir.create("./heatmaps", showWarnings = FALSE)

methods <- c("complete", "average", "single")
dists <- c("spearman", "euclidean", "maximum", "pearson")

for (m in methods) {
  pdf_file <- paste0("./heatmaps/complex_heatmap_", m, ".pdf")
  pdf(pdf_file, width = 20, height = 12)
  for (d in dists) {
    draw(Heatmap(heatmap_matrix_stand,
                 show_row_names = TRUE,
                 clustering_distance_rows = d,
                 clustering_distance_columns = d,
                 clustering_method_rows = m,
                 clustering_method_columns = m,
                 top_annotation = topAnnotation,
                 right_annotation = sideAnnotation,
                 name = "Expression (z-score)",
                 column_title = paste(m, "/", d),
                 na_col = "black"))
  }
  dev.off()
  message("Saved PDF: ", pdf_file)
}

# ----------------------------- SVG for Figure --------------------------------
dir.create("./heatmaps_svg", showWarnings = FALSE)
save_heatmap("./heatmaps_svg/complex_heatmap_complete_pearson.svg",
             method = "complete", dist = "pearson", format = "svg")

message("All heatmaps completed successfully.")
