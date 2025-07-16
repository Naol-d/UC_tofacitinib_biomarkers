# Libraries Required
library(Seurat)
library(harmony)
library(kBET)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(data.table)
library(biomaRt)
library(GenomicFeatures)
library(edgeR)
library(cluster)
library(ggplot2)
library(ggthemes)

# Random Seed for reproducibility in clustering and visualisation downstream
options(Seurat.object.seed = 263452)
set.seed(263452)

# Setting up a base directory
dataset1_base_dir = "/scratch/prj/bmb_tofacitinib/data/outputs/scRNA_aligned_ouputs/dataset_1"
dataset2_base_dir = "/scratch/prj/bmb_tofacitinib/data/outputs/scRNA_aligned_ouputs/dataset_2"

# Load single cell matrices and Create Seurat Objects for samples in dataset1
GSM5525955 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525955", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM5525955")
GSM5525956 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525956", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM5525956")
GSM5525957 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525957", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM5525957")
GSM5525963 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525963", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM5525963")

# Load single cell matrices and Creating Seurat Objects for each sample (dataset2)

GSM6614354 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614354", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614354")
GSM6614355 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614355", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614355")
GSM6614356 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614356", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614356")
GSM6614357 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614357", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614357")
GSM6614358 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614358", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614358")
GSM6614359 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614359", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614359")


# Combine all of the matrices data into one
combined_data <- merge(GSM5525955, y = list(GSM5525956, GSM5525957, GSM5525963, GSM6614354, GSM6614355, GSM6614356,
                 GSM6614357, GSM6614358, GSM6614359),
                  add.cell.ids = c("GSM5525955","GSM5525956","GSM5525957","GSM5525963",
                  "GSM6614354","GSM6614355","GSM6614356","GSM6614357","GSM6614358","GSM6614359"),
                  project = "Combined")


# A check to see how many cells there are per batch
table(combined_data$orig.ident)


############
### General Quality Control (QC) Steps:

# Computing the mitochondrial content requires a gene_symbols
# A method that is safe but somewhat has a higher space complexitiy is...
# Loading the same Seurat object of the combined sample but use argument...
# gene.column = 2 which uses the gene symbols from the features.tsv.gz
# Since the structures would be generated in the exact same order they are effectively identical

#dataset1

MT_GSM5525955 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525955", "outs", "filtered_feature_bc_matrix"), gene.column =2), project = "GSM5525955")
MT_GSM5525956 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525956", "outs", "filtered_feature_bc_matrix"), gene.column =2), project = "GSM5525956")
MT_GSM5525957 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525957", "outs", "filtered_feature_bc_matrix"), gene.column =2), project = "GSM5525957")
MT_GSM5525963 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset1_base_dir, "GSM5525963", "outs", "filtered_feature_bc_matrix"), gene.column =2), project = "GSM5525963")

#dataset2

MT_GSM6614354 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614354", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614354")
MT_GSM6614355 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614355", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614355")
MT_GSM6614356 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614356", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614356")
MT_GSM6614357 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614357", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614357")
MT_GSM6614358 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614358", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614358")
MT_GSM6614359 <- CreateSeuratObject(counts = Read10X(data.dir = file.path(dataset2_base_dir, "GSM6614359", "outs", "filtered_feature_bc_matrix"), gene.column =1), project = "GSM6614359")


MT_combined_data <- merge(MT_GSM5525955, y = list(MT_GSM5525956, MT_GSM5525957, MT_GSM5525963, MT_GSM6614354, MT_GSM6614355, MT_GSM6614356, MT_GSM6614357, MT_GSM6614358, MT_GSM6614359),
                  add.cell.ids = c("GSM5525955","GSM5525956","GSM5525957","GSM5525963",
                  "GSM6614354","GSM6614355","GSM6614356","GSM6614357","GSM6614358","GSM6614359"),
                  project = "Combined")


# Check: Check to see if all of the cell names match between the two objects
if (identical(Cells(combined_data), Cells(MT_combined_data))) {
  cat("✅ Cell names in combined_data and MT_combined_data match exactly. Mitochondrial content can be safely transferred.\n")
} else {
  cat("⚠️ Cell names in combined_data and MT_combined_data do NOT match. DO NOT proceed with mitochondrial % calculation.\n")
}

cat("======= QC Filtering has begun =======")

# As objects are identical apart from the gene naming format, we can compute the mitochondrial content using the MT_combined object using...
# MT pattern recognition then pass the percentages into the original combined_data object containing ENSEMBL IDs instead
combined_data[["percent.mt"]] <- PercentageFeatureSet(MT_combined_data, pattern = "^MT-")

# Visualise QC metrics
VlnPlot(combined_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(combined_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(combined_data, feature1 = "percent.mt", feature2 = "nFeature_RNA")


# Based on visualisations the decided metrics are as follows within the subset code
filtered_combined <- subset(combined_data, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 &
           nCount_RNA > 1000 & nCount_RNA < 50000 &
           percent.mt < 20
)

# Visualisation of the filtering
VlnPlot(filtered_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cat("\n======= Doublet Detection has begun =======\n")

###################
### Doublet Detection
## Step 1:
# Due the incompatibility between Seurat v5 data structures and scDblFinder expected format...
# we are required to extract the counts data directly into a matrix and compute doublets analysis.

# Get all layer names (e.g., counts.GSM52...)
layer_names <- Layers(filtered_combined[["RNA"]])

# Extract matrices and label cells by sample
count_list <- list()
sample_vec <- c()


# Loop to iterate through the layers to ensure cell names are unique before doublet analysis occurs
for (ln in layer_names) {
  mtx <- GetAssayData(filtered_combined[["RNA"]], layer = ln) # Access each layer counts data
  count_list[[ln]] <- mtx
  sample_vec <- c(sample_vec, rep(ln, ncol(mtx))) # Stores the sample name and the number of cells for iteration of the loop
}

# Combine matrices
combined_counts <- do.call(cbind, count_list)

# Build SCE
sce <- SingleCellExperiment(assays = list(counts = combined_counts))
colData(sce)$sample <- sample_vec  # Important for sample-aware doublet detection

# Run scDblFinder
# NOTE: Can parallise this using MulticoreParam (linux) or SnowParam (Windows), using parameter 'workers='
sce <- scDblFinder(sce, samples = "sample", BPPARAM = SerialParam(RNGseed = 789456))

# We will then inspect doublet detection results
table(sce$scDblFinder.class)


## Step 2:

# To successfully subset the data we ensure the  column names in the seurat object are the same name as the single cell experiment object.
# As Seurat's merge function, add.cell.ids will automatically add ID, in this case sample name ensuring that each column is unique
# We can directly subset the only the singlet cells

# Get vector of singlet cells
singlets <- colnames(sce)[sce$scDblFinder.class == "singlet"]

# Subset Seurat object
combined_data_singlets <- subset(filtered_combined, cells = singlets)


# Check:
# To see if number of cells matches the amount of singlets
# Expected No. of cells = 27358, based on random seed used and dataset  

# The no. of singlets were found in the scDblFinder output
num_singlets_sce <- sum(sce$scDblFinder.class == "singlet")

# The no. of cells are in the Seurat object after filtering out doublets
num_singlets_seurat <- length(colnames(combined_data_singlets))

print(paste0("Singlets in Seurat object: ", num_singlets_seurat, 
             " / Singlets detected by scDblFinder: ", num_singlets_sce))

# Optional: Warn if there's a mismatch
if (num_singlets_seurat != num_singlets_sce) {
  warning("Mismatch detected: Check if cell names were correctly synchronized before filtering.")
}

# Check to see if the column names of the two data structures are the same to avoid any errors
print(paste0(
  "Column names match between Seurat and scDblFinder singlets: ",
  identical(colnames(combined_data_singlets), colnames(sce)[sce$scDblFinder.class == "singlet"])
))

cat("\n======= SCTransform is being applied to the Seurat Object =======\n")

##################
### Apply SCTransform to Normalise, Find Feature Variables and Scale
normalised_combined_data <- SCTransform(
  combined_data_singlets,
  assay = "RNA",
  new.assay.name = "SCT",
  verbose = FALSE
)

## Asigning Dataset Classification to see if batch effects exist between the two datasets
normalised_combined_data$dataset_batch <- ifelse(
  normalised_combined_data$orig.ident %in% c("GSM5525955", "GSM5525956", "GSM5525957", "GSM5525963"),
  "dataset1", "dataset2"
)

##################
### Dimensional Reductionality (PCA)
normalised_combined_data <- RunPCA(normalised_combined_data, assay = "SCT",verbose = FALSE)

##################
### Assessing Batch Effects
p_batch <- DimPlot(
  normalised_combined_data,
  reduction = "pca",
  group.by  = "dataset_batch"
) + ggtitle("PCA colored by datasets 1 or 2")

p_sample <- DimPlot(
  normalised_combined_data,
  reduction = "pca",
  group.by  = "orig.ident"
) + ggtitle("PCA colored by samples")


ggsave(
  filename = "../../../data/outputs/visualisations/scRNA_prep/BE_figures/pre-BE/pca_by_dataset_pre_BE.png",
  plot     = p_batch,
  width    = 6,
  height   = 4,
  dpi      = 600,
  units    = "in"
)

ggsave(
  filename = "../../../data/outputs/visualisations/scRNA_prep/pre-BE/pca_by_sample_pre_BE.png",
  plot     = p_sample,
  width    = 6,
  height   = 4,
  dpi      = 600,
  units    = "in"
)

cat("\n======= Computing ASW of pre-batch corrected single cell data =======\n")

# Make a copy of the normalised_combined data test differing batch correction resolutions
BT_normalised_combined_data <- normalised_combined_data

# Extracting PCA embeddings (cells × PCs)
pc_emb <- Embeddings(normalised_combined_data, "pca")[, 1:10]
# Creating a numeric batch label
batch_lab_ds <- as.numeric(factor(normalised_combined_data$dataset_batch))
batch_lab_sample <- as.numeric(factor(normalised_combined_data$orig.ident))
# Computing silhouette widths
pca_sil_ds <- silhouette(batch_lab_ds, dist(pc_emb))
pca_sil_sample <- silhouette(batch_lab_sample, dist(pc_emb))
# Average silhouette width across all cells
pca_mean_sil_ds <- mean(pca_sil_ds[, "sil_width"])
pca_mean_sil_samples <-  mean(pca_sil_sample[, "sil_width"])
cat("Average batch‐silhouette width on PCs 1–10 Pre-Batch Correction based on dataset:", round(pca_mean_sil_ds, 3), "\n")
cat("Average batch‐silhouette width on PCs 1–10 Pre-Batch Correction based on samples:", round(pca_mean_sil_samples, 3), "\n")

#kbet_pca_ds <- kBET(pc_emb, batch_lab_ds)
#kbet_pca_sample <- kBET(pc_emb, batch_lab_sample)

#pre_rej_ds <- kbet_pca_ds$stats$avg.observed
#pre_rej_sample <- kbet_pca_sample$stats$avg.observed

##################
### Batch Correction (harmony)
cat("/n Batch Correction starting (At Dataset Resolution)... \n")
normalised_combined_data <- RunHarmony(
  object = normalised_combined_data,
  group.by.vars = "dataset_batch",   # corrects across your batches
  reduction.use = "pca",              
  assay.use = "SCT",              
  dims.use = 1:30,                
  verbose = FALSE
)


# Checking batch mixing on Harmony embeddings (at dataset resolution)
ds_harmony_embeddings <- Embeddings(normalised_combined_data, "harmony")[,1:10]
batch_lab_ds <- as.numeric(factor(normalised_combined_data$dataset_batch))
batch_lab_samples <- as.numeric(factor(normalised_combined_data$orig.ident))
ds_sil_harmony <- silhouette(batch_lab_ds, dist(ds_harmony_embeddings))
ds_sil_harmony_samples <- silhouette(batch_lab_samples, dist(ds_harmony_embeddings))

cat("Mean batch‐silhouette on Harmony by dataset (PCs 1–10, after dataset-level correction):",
    round(mean(ds_sil_harmony[, "sil_width"]), 3), "\n")
cat("Mean batch‐silhouette on Harmony by sample (PCs 1–10, after dataset-level correction):",
    round(mean(ds_sil_harmony_samples[, "sil_width"]), 3), "\n")

#batch_ds <- normalised_combined_data$dataset_batch
#kbet_ds <- kBET(ds_harmony_embeddings, batch_ds)

#ds_mean_sil   <- mean(ds_sil_harmony[, "sil_width"])
#ds_rej_rate   <- kbet_ds$stats$avg.observed

cat("/n Batch Correction starting (At a Sample resolution)... \n")
BT_normalised_combined_data <- RunHarmony(
  object = BT_normalised_combined_data,
  group.by.vars = "orig.ident",   # corrects across batches (defined as individual samples)
  reduction.use = "pca",              
  assay.use = "SCT",              
  dims.use = 1:30,                
  verbose = FALSE
)

# Checking batch mixing on Harmony embeddings (at sample resolution)
sample_harmony_embeddings <- Embeddings(BT_normalised_combined_data, "harmony")[,1:10]
BT_batch_lab_samples <- as.numeric(factor(BT_normalised_combined_data$orig.ident))
BT_batch_lab_ds <- as.numeric(factor(normalised_combined_data$dataset_batch))
sample_sil_harmony <- silhouette(BT_batch_lab_samples, dist(sample_harmony_embeddings))
sample_sil_harmony_ds <- silhouette(BT_batch_lab_ds, dist(sample_harmony_embeddings))

cat("Mean batch‐silhouette on Harmony by sample (PCs 1–10, after sample-level correction):",
    round(mean(sample_sil_harmony[, "sil_width"]), 3), "\n")

cat("Mean batch‐silhouette on Harmony by dataset (PCs 1–10, after sample-level correction):",
    round(mean(sample_sil_harmony_ds[, "sil_width"]), 3), "\n")

#batch_samples <- BT_normalised_combined_data$orig.ident
#kbet_sample <- kBET(sample_harmony_embeddings, batch_samples)

#sample_mean_sil <- mean(sample_sil_harmony[, "sil_width"])
#sample_rej_rate <- kbet_sample$stats$avg.observed

### Visualisations of Batch Corrections
#plot_batch_mixing <- function(sil_width, kbet_rej, title = "Batch-mixing performance") {
#
#  df <- data.frame(
#    metric    = c("Silhouette width", "kBET rejection"),
#    value     = c(sil_width, kbet_rej)
#  )
  
  # flip kBET so that higher=better
#  df$value_adj <- with(df,
#                       ifelse(metric == "kBET rejection",
#                             1 - value,
#                              value))
#   labels
#  df$label <- with(df,
#                  ifelse(metric == "Silhouette width",
#                         sprintf("%.2f", value),
#                         sprintf("rej=%.2f\nmix=%.2f", value, value_adj))
# )
  
  
# kbet_plots <- ggplot(df, aes(x = metric, y = value_adj, fill = metric)) +
#    geom_col(width = 0.6, show.legend = FALSE) +
#    geom_text(aes(label = label),
#             vjust = -0.5, size = 3.5, lineheight = 0.9) +
#    scale_y_continuous(
#      name   = "Batch-mixing score\n(higher = better)",
#      limits = c(0, 1),
#      expand = expansion(mult = c(0, 0.1))
#    ) +
#    labs(
#      x     = NULL,
#      title = title
#    ) +
#    theme_few(base_size = 14) +
#    theme(
#      plot.title   = element_text(hjust = 0.5, face = "bold"),
#      axis.text.x  = element_text(face = "plain"),
#      axis.title.y = element_text(margin = margin(r = 10)),
#      panel.grid   = element_blank()
#    )
  
#  return(kbet_plots)
#}

#p_pre_sample <- plot_batch_mixing(pca_mean_sil_samples, pre_rej_sample, "Pre-Batch Correction (samples defined as batch)")
#p_pre_ds <- plot_batch_mixing(pca_mean_sil_ds, pre_rej_ds, "Pre-Batch Correction (datasets defined as batch)")
#p_post_ds <- plot_batch_mixing(ds_mean_sil, ds_rej_rate, "Post-Harmony BC (dataset batches)")
#p_post_sample <- plot_batch_mixing(sample_mean_sil, sample_rej_rate, "Post-Harmony BC (sample batches)")

#save_plots_png <- function(plots, dir, width = 6,height = 4, dpi = 600) {
#  if (!dir.exists(dir)) {
#    stop("Directory does not exist: ", dir)
#  }
#  for (nm in names(plots)) {
#    file <- file.path(dir, paste0(nm, ".png"))
#    ggsave(
#     filename = file,
#      plot     = plots[[nm]],
#      width    = width,
#      height   = height,
#      dpi      = dpi,
#      units    = "in"
#    )
#    message("Saved: ", file)
#  }
#}

#all_plots <- list(
#  pre_sample   = p_pre_sample,
#  pre_dataset  = p_pre_ds,
#  post_dataset = p_post_ds,
#  post_sample  = p_post_sample
#)

#save_plots_png(
#  plots = all_plots,
#  dir   = "../../../data/outputs/visualisations/scRNA_prep/BC_Assessment/",  
#  width = 7,
#  height = 4,
#  dpi = 600
#)


cat("\n============ Clustering has begun ============\n")

##################
### Clustering
normalised_combined_data <- FindNeighbors(normalised_combined_data,reduction = "harmony", dims = 1:30)
normalised_combined_data <- FindClusters(normalised_combined_data, resolution = 0.5)

# Visualisation using UMAP
normalised_combined_data <- RunUMAP(normalised_combined_data,reduction = "harmony", dims = 1:30, seed.use=263452,  verbose = FALSE)
cell_clustering_plot <- DimPlot(normalised_combined_data,reduction = "umap", label = TRUE, label.size = 4) + NoLegend()

ggsave(
  filename = "../../../data/outputs/visualisations/scRNA_prep/UMAP_cell_clustering_no_labels.png",
  plot     = cell_clustering_plot,
  width    = 6,
  height   = 4,
  dpi      = 600,
  units    = "in"
)

# Visualisation by Batch
batch_corrected_plot_bydataset <- DimPlot(
    normalised_combined_data,
    reduction = "umap",
    group.by = "dataset_batch"
  ) + ggtitle("UMAP (Harmony) by dataset 1 or 2")

ggsave(
  filename = "../../../data/outputs/visualisations/scRNA_prep/BE_figures/post-BE/UMAP_cell_cluster_visualisation_post-BE.png",
  plot     = batch_corrected_plot_bydataset,
  width    = 6,
  height   = 4,
  dpi      = 600,
  units    = "in"
)

cat("\n============ Cell Marker Identification has begun ============\n")

##################
### Cell Marker Identification
cat(" \n Cell Marker Identification Beginning... \n")
normalised_combined_data <- PrepSCTFindMarkers(
  normalised_combined_data,
  assay = "SCT",
)

cat ("\n Identifying all markers \n")
cluster_markers <- FindAllMarkers(
  normalised_combined_data,
  assay = "SCT",
  slot = "data",        # Pearson residuals
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  random.seed = 263452
)


# Subset markers to the top 15 markers per cluster
top100 <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100)

### Assess if any genes in bulkRNA (excluded due to batch influence in variance being over 8.3% or the top 5% of genes) are present as the top 30 markers
batch_related_genes <- readLines("../../../data/outputs/bulkRNA_processed/batch_top5pct_genes.txt")

top100_ranked <-
  top100 %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%    # strongest → weakest
  mutate(marker_rank = row_number()) %>%    # 1..100 per cluster
  ungroup() %>%
  mutate(overall_row = row_number())        # 1..nrow(top100)

overlap_details <-
  top100_ranked %>%
  filter(gene %in% batch_related_genes) %>%
  dplyr::select(overall_row, cluster, marker_rank, gene, avg_log2FC)
print(overlap_details)


## Ensembl Configuration to prevent a timeout
options(timeout = 120)               
httr::set_config(httr::timeout(120))  

#### Converting Ensembl IDs into gene symbols for literature review
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113)

# Get rid of Ensembl versioning
top100$gene <- sub("\\.\\d+$", "", top100$gene)
ensembl_ids <- top100$gene

# Ensembl Search and outputs a mapping of ensembl ID to gene symbols
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Checks:
print("Mapping details...")
nrow(mapping)
head(mapping)
setdiff(ensembl_ids, mapping$ensembl_gene_id)

# Merge the mapping with the top100 data structure
top100_mapped <- top100 %>%
  left_join(mapping, by = c("gene" = "ensembl_gene_id")) %>%
  mutate(gene_symbol = ifelse(hgnc_symbol == "", gene, hgnc_symbol))

# Remove the redudant rows and what is left is a gene_symbol column that will have
# both gene symbols and a fall back to the a ensembl ID if gene symbol doesn't exist
top100_mapped <- top100_mapped[, !(names(top100_mapped) %in% c("hgnc_symbol", "gene"))]

# Prints out the top 100 markers
cat("\n===  top 100 makers ===\n")
top100_mapped %>%
  group_by(cluster) %>%
  summarise(gene_symbol = paste(unique(gene_symbol), collapse = ", "), .groups="drop") %>%
  arrange(cluster) %>%
  mutate(line = paste0(cluster, ": ", gene_symbol)) %>%
  pull(line) %>%
  cat(sep = "\n")

print("#####################")

cat("\n============ Assignment of cell types has begun ============\n")

##################
### Assigning cell types to the cell names, post literature review

# Check current assigned names
levels(Idents(normalised_combined_data))

cluster_labels <- c(
  `0`  = "CD4_Treg_cell",
  `1`  = "CD8_T_cell",
  `2`  = "CD4_Central_Memory_T_cell",
  `3`  = "Fibroblast_cell",
  `4`  = "Inflammatory_Macrophage_cell",
  `5`  = "Plasma_cell",
  `6`  = "CD4_Tfh_cell",
  `7`  = "Plasma_cell",
  `8`  = "Plasmablast_cell",
  `9`  = "Naive_B_cell",
  `10` = "Memory_B_cell",
  `11` = "Naive_B_cell",
  `12` = "Plasma_cell",
  `13` = "Memory_B_cell",
  `14` = "Epithelial_cell",
  `15` = "Plasma_cell",
  `16` = "Plasma_cell",
  `17` = "Plasmablast_cell",
  `18` = "Endothelial_cell",
  `19` = "CD4_Naive_T_cell",
  `20` = "Plasma_cell",
  `21` = "Plasma_cell",
  `22` = "Memory_B_cell",
  `23` = "Plasma_cell",
  `24` = "Plasma_cell",
  `25` = "Cycling_S_G2M_cell",
  `26` = "Plasma_cell",
  `27` = "Plasmablast_cell",
  `28` = "Plasmablast_cell",
  `29` = "Plasma_cell",
  `30` = "Plasma_cell",
  `31` = "Mast_cell",
  `32` = "Plasma_cell",
  `33` = "Pericyte_cell",
  `34` = "Plasma_cell",
  `35` = "Plasma_cell",
  `36` = "Plasma_cell",
  `37` = "Plasma_cell",
  `38` = "Memory_B_cell",
  `39` = "Plasma_cell",
  `40` = "Plasma_cell",
  `41` = "Plasma_cell",
  `42` = "Plasmablast_cell",
  `43` = "Plasma_cell",
  `44` = "Plasma_cell"
)


print(cluster_labels)


# Check: To see if seurat_clusters has the same number of unique values as the clusters
# Expected Output = 45
length(unique(normalised_combined_data@meta.data$seurat_clusters))

# Step 1:

# Create a new column in the meta.data of the Seurat object called cell_type.
# Use a named vector, where the names correspond to Seurat cluster IDs (seurat_clusters) and the values are cell type labels, we map each cell’s cluster to its corresponding label.
# The label from the named vector is assigned to the cell_type column for each row (cell) where the seurat_clusters value matches the name in the vector.

normalised_combined_data@meta.data$cell_type <- cluster_labels[
  as.character(normalised_combined_data$seurat_clusters)
]


Idents(normalised_combined_data) <- "cell_type"

# Check: To see if a cell was assigned a cell_type correctly
cat("=== Cell Type Assignment Summary ===\n")
table(normalised_combined_data$cell_type)

# Check: To see if the total number of cells labelled equals the number of cells in the object
cat("\n=== Cell Type Labelling Consistency Check ===\n")
sum(table(normalised_combined_data$cell_type)) == ncol(normalised_combined_data)

cat("\n============ Creating reference files for Cibersortx Signature matrix generation has begun ============\n")

##################
### Creating Sample Reference file
expr_matrix <- GetAssayData(normalised_combined_data, assay = "RNA", layer = "counts")
cell_types <- normalised_combined_data@meta.data$cell_type
names(cell_types) <- colnames(normalised_combined_data)

# Check: Ensuring that the cell names match between counts and annotations 
cat("\n=== Cell names match between expr_matrix and cell_types ===\n")
all(colnames(expr_matrix) == names(cell_types))  # Should return TRUE


## Preparing the data into a delimited format for Cibersortx input

# Convert the sparse matrix into a dense one
dense_expr_matrix <- as.matrix(expr_matrix)

# Remove all of the zero values
dense_expr_matrix <- dense_expr_matrix[rowSums(dense_expr_matrix) > 0, ]

# Assign the cell_types names to the columns replacing the UMI barcode headers instead
colnames(dense_expr_matrix) <- as.character(cell_types)


############## Creating a Raw Reference File Version ##############
raw_dt_expr <- as.data.table(dense_expr_matrix)

# Bind the genes into a column and combine them with the data.table
raw_dt_expr <- cbind(Gene = rownames(dense_expr_matrix), raw_dt_expr)

# Get rid of the ensembl ID version numbers
raw_dt_expr$Gene <- sub("\\.\\d+$", "", raw_dt_expr$Gene)

fwrite(raw_dt_expr, file = "../../../data/outputs/scRNA_processed/integrated_raw_count_CIBERSORTx_reference.txt", sep = "\t",
       quote = FALSE, col.names = TRUE)


############## Creating a RPKM Reference Version ##############
##### Computing RPKM values Steps
txdb <- makeTxDbFromGFF("/scratch/prj/bmb_tofacitinib/data/reference_GRCh38_gencode_v47/gencode.v47.primary_assembly.annotation.gtf", format="gtf")
exonsByGene <- exonsBy(txdb, by="gene")
geneLengths <- sum(width(reduce(exonsByGene)))
gene_lengths_df <- data.frame(Geneid = names(geneLengths), Length = as.integer(geneLengths))

# Subset the gene length mappings to only include the geneID that exist after filtering for non-zero counts across all cells
subsetted_gene_lengths_df <- gene_lengths_df[gene_lengths_df$Geneid %in% rownames(dense_expr_matrix), ]

# Confirm the counts match:
#   If some versioned IDs were removed by the zero‐sum filter, subsetted_gene_lengths_mapping should lose exactly those rows.
if (nrow(subsetted_gene_lengths_df) != length(rownames(dense_expr_matrix))) {
  stop(
    "Mismatch: after filtering, some gene IDs in dense_expr_matrix are not found in gene_lengths_df$Geneid."
  )
}

# Re‐order subsetted_gene_lengths_mapping so that aligned_gene_length_mapping$Geneid[i] == subsetted_gene_lengths_mapping[i]:
gene_mapping <- subsetted_gene_lengths_df[match(rownames(dense_expr_matrix), subsetted_gene_lengths_df$Geneid), ]

# Sanity‐check:
if (!all(gene_mapping$Geneid == rownames(dense_expr_matrix))) {
  stop("❌ Mismatch detected: 'gene_mapping$Geneid' does not align with 'dense_expr_matrix gene IDs'. Please verify that gene identifiers are correctly ordered and matched.")
}

# NOTE: Since the rownames and aligned_gene_length_mapping are in the same order we can compute RPKMs efficiently using vectorisation

#############
### Step 2. Computing RPKM values using correct genes

# Storing gene lengths (bp) in a vector
gene_length_bp <- gene_mapping$Length
names(gene_length_bp) <- gene_mapping$Geneid

cat("==== Printing out the names of gene_length_bp and rownames of dense_expr_matrix as test====\n")
head(names(gene_length_bp))
head(rownames(dense_expr_matrix))

if (all(rownames(dense_expr_matrix) == names(gene_length_bp))) {
  cat("✅ Alignment Check Passed: Row names of dense_expr_matrix match names of gene_length_bp.\n")
} else {
  cat("❌ Alignment Check Failed: Row names of dense_expr_matrix do NOT match names of gene_length_bp.\n")
}

#  total number of reads per cell
library_sizes <- colSums(dense_expr_matrix)

# Check that lengths match matrix dimensions
if (length(gene_length_bp) != nrow(dense_expr_matrix)) {
  stop("Length of gene_length_bp does not match the number of rows (genes) in the expression matrix. CANNOT compute RKPMs")
}

if (length(library_sizes) != ncol(dense_expr_matrix)) {
  stop("Length of library_sizes does not match the number of columns (cells/samples) in the expression matrix. CANNOT compute RKPMs")
}

# Compute RPKM using standard formula: (1e9 * count) / (gene_length_bp * library_size)
rpkm_matrix <- (dense_expr_matrix * 1e9) / outer(gene_length_bp, library_sizes, "*")

#### A quick check to see if edgeR RKPMs computation matches rkpm matrix
dge <- DGEList(counts = dense_expr_matrix)
edgeR_rpkm <- rpkm(dge, gene.length = gene_length_bp)

max_abs_diff <- max(abs(rpkm_matrix - edgeR_rpkm))
if (isTRUE(all.equal(rpkm_matrix, edgeR_rpkm))) {
  cat("✅ Manual RPKM matches edgeR rpkm (within machine precision).\n")
} else {
  cat("❌ Discrepancy—max abs difference =", max_abs_diff, "\n")
}

# Convert the matrix into a data.table (data.table version of dataframe)
rpkm_dt_expr <- as.data.table(rpkm_matrix)

# Bind the genes into a column and combine them with the data.table
rpkm_dt_expr <- cbind(Gene = rownames(dense_expr_matrix), rpkm_dt_expr)

# Get rid of the ensembl ID version numbers 
rpkm_dt_expr$Gene <- sub("\\.\\d+$", "", rpkm_dt_expr$Gene)


fwrite(rpkm_dt_expr, file = "../../../data/outputs/scRNA_processed/integrated_RKPM_count_CIBERSORTx_reference.txt", sep = "\t",
       quote = FALSE, col.names = TRUE)


# Final check
# Ensuring that no duplicate headers exist
cat("\n=== Duplicate Gene Check ===\n")
if (any(duplicated(rownames(expr_matrix)))) {
  cat("❌ Duplicate gene names found! Please remove or resolve them.\n")
} else {
  cat(" ✅ No duplicate gene names — all rows are unique.\n")
}


# Ensure that number of columns matches the number of singlets
cat("\n=== Column Count Consistency Check ===\n")
num_columns_in_output <- ncol(dense_expr_matrix)
num_singlets_in_seurat <- sum(normalised_combined_data$cell_type %in% unique(cell_types))

cat("Columns in expression matrix:", num_columns_in_output, "\n")
cat("Number of singlet cells with cell type:", num_singlets_in_seurat, "\n")

if (num_columns_in_output == num_singlets_in_seurat) {
  cat("✅ Column count matches number of singlet cells — ready for CIBERSORTx.\n")
} else {
  cat("❌ Mismatch in column count — something may be wrong.\n")
}


################### Optional Sections ########################


### Optional 
# Save data in RDS file for further investigations

saveRDS(cluster_markers, file = "../../../data/outputs/scRNA_processed/RDS_files/integrated_all_cluster_markers.rds")

# Save top 30 markers per cluster
saveRDS(top100, file = "../../../data/outputs/scRNA_processed/RDS_files/integrated_top100_cluster_markers.rds")

# Save Seurat object
saveRDS(normalised_combined_data, file = "../../../data/outputs/scRNA_processed/RDS_files/integrated_normalised_combined_data.rds")

saveRDS(dense_expr_matrix, file = "../../../data/outputs/scRNA_processed/RDS_files/integrated_dense_expr_matrix.rds")



stratified_downsample <- function(dt_expr,
                                  keep_fraction = 0.35,
                                  seed = 263452) {
  ############### Downsampling Stratified ###############
  ## Used to generate a Sig Matrix from the

  # Total number of cells (columns − 1 for “Gene”)
  all_cell_cols <- colnames(dt_expr)[-1]
  N_total_cells <- length(all_cell_cols)
  cat("Total cells available in dt_expr:", N_total_cells, "\n\n")

  # ── 2) Decide on a target total number of cells to keep ───────────────────────────────
  set.seed(seed)
  target_total  <- round(N_total_cells * keep_fraction)
  pct <- round(keep_fraction * 100, 1)
  cat("Sampling", target_total, "cells out of", N_total_cells,
      "(", pct, "% )\n\n")

  # ── 3) Compute how many cells per cell type to sample ─────────────────────────────────
  counts_per_type <- table(all_cell_cols)
  cat("Original counts per cell type:\n")
  print(counts_per_type)
  cat("\n")

  # Compute proportional targets
  props           <- counts_per_type / sum(counts_per_type)
  target_per_type <- round(props * target_total)

  # Adjust for rounding so that sum(target_per_type) == target_total
  delta <- target_total - sum(target_per_type)
  if (delta != 0) {
    # Order types by abundance, descending
    ordered_types <- names(sort(counts_per_type, decreasing = TRUE))
    # If delta > 0, add 1 to the first `delta` types; if < 0, subtract 1
    if (delta > 0) {
      for (i in seq_len(abs(delta))) {
        target_per_type[ordered_types[i]] <- target_per_type[ordered_types[i]] + 1
      }
    } else {
      for (i in seq_len(abs(delta))) {
        tname <- ordered_types[i]
        target_per_type[tname] <- max(0, target_per_type[tname] - 1)
      }
    }
  }

  cat("Adjusted cells to keep per type:\n")
  print(target_per_type)
  cat("Total cells after adjustment:", sum(target_per_type), "\n\n")
  stopifnot(sum(target_per_type) == target_total)

  # ── 4) For each cell type, randomly sample that many column indices ───────────────────
  sampled_indices <- integer(0)

  for (ct in names(target_per_type)) {
    # Find all indices in dt_expr whose column name == ct
    idx_all <- which(colnames(dt_expr) == ct)
    n_keep  <- target_per_type[ct]
    if (length(idx_all) <= n_keep) {
      sampled_indices <- c(sampled_indices, idx_all)
    } else {
      sampled_indices <- c(
        sampled_indices,
        sample(idx_all, size = n_keep, replace = FALSE)
      )
    }
  }

  # Ensure we keep column 1 (“Gene”) plus the sorted sampled indices
  sampled_indices <- sort(sampled_indices)
  final_cols      <- c(1L, sampled_indices)

  down_dt_expr <- dt_expr[, ..final_cols]

  # Sanity checks
  n_cells_down <- ncol(down_dt_expr) - 1
  cat("Number of cells in downsampled table:", n_cells_down, "\n")
  stopifnot(n_cells_down == target_total)

  cat("Counts per type in downsampled table:\n")
  print(table(colnames(down_dt_expr)[-1]))
  cat("\n")

  #### Checks: build comparison table but do not write to disk

  # Counts per type in the original dt_expr
  orig_cols   <- colnames(dt_expr)[-1]   # drop “Gene” column
  orig_counts <- as.data.table(table(orig_cols), keep.rownames = TRUE)
  setnames(orig_counts, c("CellType", "Orig_Count"))

  # Counts per type in the downsampled table
  down_cols   <- colnames(down_dt_expr)[-1]
  down_counts <- as.data.table(table(down_cols), keep.rownames = TRUE)
  setnames(down_counts, c("CellType", "Down_Count"))

  # Merge them (fills with 0 if a type was completely dropped, though that shouldn't happen)
  comparison <- merge(
    orig_counts, down_counts,
    by    = "CellType",
    all.x = TRUE,
    all.y = TRUE
  )

  # Replace any NA with 0
  comparison[is.na(Orig_Count), Orig_Count := 0]
  comparison[is.na(Down_Count), Down_Count := 0]

  # 4) Add fraction columns
  total_orig <- sum(comparison$Orig_Count)
  total_down <- sum(comparison$Down_Count)

  comparison[, Orig_Frac := Orig_Count / total_orig ]
  comparison[, Down_Frac := Down_Count / total_down ]

  # Round values
  comparison[, `:=`(
    Orig_Frac  = round(Orig_Frac, 3),
    Down_Frac  = round(Down_Frac, 3)
  )]
  comparison <- comparison[order(-Orig_Count)]  # sort by original abundance
  print(comparison)

  return(list(
    down_dt_expr = down_dt_expr,
    comparison   = comparison
  ))
}

# Downsizing the raw matrix
res_raw <- stratified_downsample(
  dt_expr       = raw_dt_expr,
  keep_fraction = 0.35,
  seed          = 263452
)
down_raw   <- res_raw$down_dt_expr
summary_raw <- res_raw$comparison

# Downsampling the RPKM matrix 
res_rpkm <- stratified_downsample(
  dt_expr       = rpkm_dt_expr,
  keep_fraction = 0.35,
  seed          = 263452
)
down_rpkm   <- res_rpkm$down_dt_expr
summary_rpkm <- res_rpkm$comparison



fwrite(down_raw, file= ".../../../data/outputs/scRNA_processed/integrated_raw_stratified_downsampled_reference.txt", sep= "\t", quote= FALSE, col.names = TRUE)
fwrite(down_rpkm, file= "../../../data/outputs/scRNA_processed/integrated_rpkm_stratified_downsampled_reference.txt", sep= "\t", quote= FALSE, col.names = TRUE)














