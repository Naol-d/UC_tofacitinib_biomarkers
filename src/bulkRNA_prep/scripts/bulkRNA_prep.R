# Libraries
library(DESeq2)
library(sva)
library(ggplot2)
library(variancePartition)

# Load files 
raw_counts <-read.table("../../../data/outputs/bulkRNA_aligned_outputs/FC_combo.txt", header = TRUE, sep = "\t")
metadata <- read.table("../../../data/originals/bulkRNA/bulkRNA_sample_metadata/Meta_data_2.csv", header = TRUE, sep = ",")

# Clean the column names so that only the sample names appear (previously was just the full directory from the HPC)
clean_colnames <- function(x) {
  gsub(".*disease\\.(.*?)_Aligned.*", "\\1", x)
}
colnames(raw_counts)[7:ncol(raw_counts)] <- clean_colnames(colnames(raw_counts)[7:ncol(raw_counts)])
print(colnames(raw_counts))

# Output the cleaned FC_combo file for storage
write.table(raw_counts, "../../../data/outputs/bulkRNA_processed/FC_combo_latest.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Clean the raw_counts table so it becomes a raw count matrix
lengths <- raw_counts$Length
counts <- raw_counts[, 7:ncol(raw_counts)]  # Adjust column range as needed
rownames(counts) <- raw_counts$Geneid  # Set gene IDs as row names

# Generate a raw matrix txt for Deconvolution Purposes / Ensembl Gene ID depreciation check
write.table(counts, "../../../data/outputs/bulkRNA_processed/bulk_rnaseq_raw_count.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# Remove BR18A as there was a error with the sample metadata handling
counts <- counts[, !colnames(counts) %in% c("BR18A")] 
metadata <- metadata[!metadata$ID %in% c("BR18A"), ] 

# Creating a metadata_df which will be used for batch effects assessments
metadata_df <- data.frame(
  row.names = metadata$ID,  
  Batch = metadata$batch,   
  Condition = metadata$condition,
  Sex = metadata$Sex
)

metadata_df$Batch <- as.factor(metadata_df$Batch)
metadata_df$Condition <- as.factor(metadata_df$Condition)
metadata_df$Sex <- as.factor(metadata_df$Sex)
str(metadata_df)

# Check if they are now properly aligned
cat("Check if rownames of metadata_df match colnames of counts: ", all(rownames(metadata_df) == colnames(counts)), "\n")



#################
# Assessing Batch Effects on raw data
dds_raw <- DESeqDataSetFromMatrix(countData = counts, colData = metadata_df, design = ~ Batch + Condition + Sex)
dds_raw <- estimateSizeFactors(dds_raw)

# VST
vst_raw <- vst(dds_raw, blind = TRUE)
vst_raw_mat <- assay(vst_raw)

# Remove zero variance genes using a SD=0
zero_var_raw <- apply(vst_raw_mat, 1, function(x) sd(x) == 0)
vst_raw_filtered <- vst_raw_mat[!zero_var_raw, ]

pca_raw <- prcomp(t(vst_raw_filtered), scale. = TRUE)

pca_df_raw <- data.frame(
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  PC3 = pca_raw$x[, 3],  
  PC4 = pca_raw$x[, 4],
  PC5 = pca_raw$x[, 5],
  PC6 = pca_raw$x[, 6],
  Batch = metadata_df$Batch,
  Condition = metadata_df$Condition,
  Sex = metadata_df$Sex
)

# For visualisation labels are added to the duplicate samples to track their projections
label_samples <- c("BR10A", "BR13A", "2_BR10A", "2_BR13A")
pca_df_raw$Label <- ifelse(rownames(pca_df_raw) %in% label_samples, rownames(pca_df_raw), NA)

plot_pca_with_links <- function(pc_x, pc_y, pc_df, title, filename, width = 6, height = 4, dpi = 600, units = "in") {
  # Stores the coordinates for the duplicated samples to help connect downstream
  pair_lines <- data.frame(
    Sample1 = c("BR10A", "BR13A"),
    Sample2 = c("2_BR10A", "2_BR13A"),
    x1 = pc_df[c("BR10A", "BR13A"), pc_x],
    y1 = pc_df[c("BR10A", "BR13A"), pc_y],
    x2 = pc_df[c("2_BR10A", "2_BR13A"), pc_x],
    y2 = pc_df[c("2_BR10A", "2_BR13A"), pc_y]
  )
  
  graph_with_links <- ggplot(pc_df, aes_string(x = pc_x, y = pc_y, shape = "Condition", fill = "Batch", color = "Sex")) +
    geom_point(size = 3, alpha = 0.8, stroke = 1) +
    ggrepel::geom_text_repel(aes(label = Label), color = "black", size = 3, na.rm = TRUE) +
    geom_segment(data = pair_lines,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 inherit.aes = FALSE,
                 color = "black", linetype = "dashed") +
    scale_color_manual(name= "Sex", values = c("Male" = "blue", "Female" = "red")) +
    scale_fill_manual(name = "Batch", breaks = c("Batch_1", "Batch_2"), labels = c("Batch 1", "Batch 2"), values = c("Batch_1" = "green", "Batch_2" = "orange")) +
    scale_shape_manual(name = "Condition", values = c("Relapse" = 21, "Remission" = 24)) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 21,
          fill  = NA,
          stroke = 1
        )
      ),
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          color = "black",
          stroke = 0.5
        )
      ),
      shape = guide_legend(
        override.aes = list(
          fill = "grey",
          color = "black"
        )
      )
    ) +
    theme_minimal() +
    ggtitle(title)

  ggsave(filename = filename,
         plot     = graph_with_links,
         width    = width,
         height   = height,
         dpi      = dpi,
         units    = units)

}

plot_pca_with_links("PC1", "PC2", pca_df_raw, "PCA1 vs PCA2 BEFORE Batch Correction", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/pre-BE/PCA1vsPCA2_preBE.png")
plot_pca_with_links("PC3", "PC4", pca_df_raw, "PCA3 vs PCA4 BEFORE Batch Correction", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/pre-BE/PCA3vsPCA4_preBE.png")
plot_pca_with_links("PC5", "PC6", pca_df_raw, "PCA5 vs PCA6 BEFORE Batch Correction", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/pre-BE/PCA5vsPCA6_preBE.png")

cat("\n======== PCA Summary Stats Before BATCH CORRECTION===========\n\n")
summary(pca_raw)$importance
raw_scree_plot <- data.frame(PC = 1:length(pca_raw$sdev),
                         Variance = (pca_raw$sdev)^2 / sum((pca_raw$sdev)^2))

raw_scree_visual <- ggplot(raw_scree_plot, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity") +
    geom_line(aes(y = cumsum(Variance)), col = "red") +
    theme_minimal() +
    labs(title = "RNA-seq PCA: Variance Explained by Principal Components",
         y = "Proportion of Variance Explained",
         x = "Principal Component")

ggsave(
  filename = "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/pre-BE/raw_Scree-plot.png",
  plot = raw_scree_visual,
  width    = 6,
  height   = 4,
  dpi      = 600
)


# Removing the duplicated samples before batch correction to ensure that Combat-seq's batch correction isn't biased
counts <- counts[, !colnames(counts) %in% c("2_BR10A", "2_BR13A")]
metadata <- metadata[!metadata$ID %in% c("2_BR10A", "2_BR13A"),]

cat("Removed duplicate samples BR10 and BR13. Now Checking if rownames of metadata match colnames of counts.. OUTPUT: ", all(metadata$ID == colnames(counts)), "\n")


# Creating a second metadata_df which will be used for batch Correction and Normalisation updated with the duplicated sample removal
metadata_df2 <- data.frame(
  row.names = metadata$ID,  # Assign sample IDs as row names
  Batch = metadata$batch,   # Extract batch information
  Condition = metadata$condition,  # Extract condition column
  Sex = metadata$Sex
)

metadata_df2$Batch <- as.factor(metadata_df2$Batch)
metadata_df2$Condition <- as.factor(metadata_df2$Condition)
metadata_df2$Sex <- as.factor(metadata_df2$Sex)
str(metadata_df2)

# Ensure metadata is correctly ordered to match counts matrix column order
# metadata_df2 <- metadata_df2[colnames(counts), , drop = FALSE]

# Check if they are now properly aligned
cat("Check if rownames of metadata_df2 match colnames of counts: ", all(rownames(metadata_df2) == colnames(counts)), "\n")

###########  Batch Correction
# Extract batch and condition information from metadata
batch <- as.factor(metadata$batch)
condition <- as.factor(metadata$condition)

# Apply ComBat-Seq
adjusted_counts <- ComBat_seq(counts = as.matrix(counts), batch = batch, group = condition)

# Output batch Corrected Bulk RNA Matrix, Not Normalised for downstream analysis
write.table(adjusted_counts, "../../../data/outputs/bulkRNA_processed/batch_corrected_bulkRNA_matrix.txt", sep="\t", row.names = TRUE, quote = FALSE)

############ Normalisation 
# Normalising counts using DESeq2 + additional batch effect mitigation strategy added
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts, colData = metadata_df2, design = ~ Batch + Condition + Sex)

# Normalisation using size factors
dds <- estimateSizeFactors(dds) 
norm_counts <- counts(dds, normalized = TRUE)


# Stripping version number from Ensemble gene IDs
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df <- cbind(Geneid = rownames(norm_counts_df), norm_counts_df)

norm_counts_df$Geneid <- sub("\\.\\d+$", "", norm_counts_df$Geneid)
norm_counts_df <- cbind(norm_counts_df[, 1, drop = FALSE], description = NA, norm_counts_df[, -1])


# Check to make sure number of unique gene_ids stripped of their version numbers are equal to the total length of the rows
# Get unique Gene IDs
unique_gene_ids <- unique(norm_counts_df$Geneid)

# Count the number of unique Gene IDs
num_unique_gene_ids <- length(unique_gene_ids)

# Number of rows in raw counts e.g. number of genes
raw_count_row_length <- nrow(raw_counts)

# Output check message
cat("Data Integrity Check:\n")
cat("  - Unique Gene IDs in normalized counts: ", num_unique_gene_ids, "\n")
cat("  - Total rows in raw counts matrix:      ", raw_count_row_length, "\n")

# Optional: warning if there's a mismatch
if (num_unique_gene_ids != raw_count_row_length) {
  cat("WARNING: Number of unique Gene IDs does not match number of rows in raw counts!\n")
} else {
  cat("Check Passed: Unique Gene IDs match number of rows in raw counts.\n")
}


############# Assessing Batch Effect Strategy using PCA

# Variance Stablisation Transformation
vst_counts <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_counts)

# Removing rows with no variance
zero_var_samples <- apply(vst_matrix, 1, function(x) sd(x) == 0)
num_zero_var <- sum(zero_var_samples)

if (num_zero_var > 0) {
  cat("Number of samples with zero variance:", num_zero_var, "\n")
}

vst_mat_filtered <- vst_matrix[ !zero_var_samples, ] 

cat("\nThe dimensions of vst_mat_filtered:", paste(dim(vst_mat_filtered), collapse = " x "), "\n")

pca_result <- prcomp(t(vst_mat_filtered), scale. = TRUE)

pca_df <- data.frame(PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2],
                     PC3 = pca_result$x[, 3],  
                     PC4 = pca_result$x[, 4],
                     PC5 = pca_result$x[, 5],
                     PC6 = pca_result$x[, 6],
                     Batch = metadata_df2$Batch,
                     Condition = metadata_df2$Condition,
                     Sex = metadata_df2$Sex)

plot_normalized_pca <- function(pc_x, pc_y, pca_df, title, filename, width = 6, height = 4, dpi = 300, units = "in") {
  p <- ggplot(pca_df, aes_string(x = pc_x, y = pc_y, shape = "Condition", fill = "Batch", color = "Sex")) +
    geom_point(size = 3, alpha = 0.8, stroke = 1) +
    scale_color_manual(name= "Sex", values = c("Male" = "blue", "Female" = "red")) +
    scale_fill_manual(name = "Batch", breaks = c("Batch_1", "Batch_2"), labels = c("Batch 1", "Batch 2"), values = c("Batch_1" = "green", "Batch_2" = "orange")) +
    scale_shape_manual(name = "Condition", values = c("Relapse" = 21, "Remission" = 24)) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 21,
          fill  = NA,
          stroke = 1
        )
      ),
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          color = "black",
          stroke = 0.5
        )
      ),
      shape = guide_legend(
        override.aes = list(
          fill = "grey",
          color = "black"
        )
      )
    ) +
    theme_minimal() +
    ggtitle(title)

  ggsave(filename = filename,
         plot     = p,
         width    = width,
         height   = height,
         dpi      = dpi,
         units    = units)

}

plot_normalized_pca("PC1", "PC2", pca_df, "Batch Corrected: PCA1 vs PCA2", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/PCA1vsPCA2_postBE.png")
plot_normalized_pca("PC3", "PC4", pca_df, "Batch Corrected: PCA3 vs PCA4", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/PCA3vsPCA4_postBE.png")
plot_normalized_pca("PC5", "PC6", pca_df, "Batch Corrected: PCA5 vs PCA6", "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/PCA5vsPCA6_postBE.png")



# Assessing what sample seems to be the outlier
outlier_sample <- rownames(pca_result$x)[pca_result$x[,1] > 400]
print(outlier_sample)

cat("\n======== PCA Summary Stats AFTER COMBAT-SEQ BATCH CORRECTION===========\n\n")

summary(pca_result)
scree_plot <- data.frame(PC = 1:length(pca_result$sdev),
                         Variance = (pca_result$sdev)^2 / sum((pca_result$sdev)^2))

scree_visual <- ggplot(scree_plot, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity") +
    geom_line(aes(y = cumsum(Variance)), col = "red") +
    theme_minimal() +
    labs(title = "Batch Corrected RNA-seq PCA: Variance Explained by PCs",
         y = "Proportion of Variance Explained",
         x = "Principal Component")

ggsave(
  filename = "../outputs/visualisations/BE_figures/post-BE/BC_Scree-plot.png",
  plot = scree_visual,
  width    = 6,
  height   = 4, 
  dpi      = 600 
)

# Variance Parition on the pca to assess if the batch effects have been overestimated
pca_scores       <- pca_result$x    # samples × PCs
expr_for_varPart <- t(pca_scores)   # PCs × samples

# sanity check
stopifnot(all(colnames(expr_for_varPart) == rownames(metadata_df2)))

# 2) Computing per‐PC variance and filtering out near‐zero PCs
pc_vars        <- rowVars(expr_for_varPart)
keep           <- pc_vars > 1e-05
expr_filtered  <- expr_for_varPart[keep, ]

dropped <- rownames(expr_for_varPart)[!keep]
message("Dropped zero‐variance PCs: ", paste(dropped, collapse = ", "))

# 3) Fit the variance‐partition model on the **filtered** matrix
pca_var_explained_filtered <- fitExtractVarPartModel(
  formula = ~ Batch + Condition + Sex,
  exprObj = expr_filtered,
  data    = metadata_df2
)

# 4) Plot and save
pca_vp_plot <- plotVarPart(pca_var_explained_filtered)

ggsave("../outputs/visualisations/BE_figures/post-BE/PCA_Var_Partition1.png",
       plot   = pca_vp_plot,
       width  = 6,
       height = 4,
       dpi    = 600)


############# Assessing Batch Effect Mitigation Efforts using variance partition

# Ensure metadata is a data frame
metadata_df2 <- as.data.frame(metadata_df2)

# Check dimensions to match rows in metadata with columns in the expression matrix
if (!all(colnames(vst_matrix) == rownames(metadata_df2))) {
  stop("Sample names in vst_matrix and metadata must match!")
}

# Fit variance partition model
varPart <- fitExtractVarPartModel(
  formula = ~ Batch + Condition + Sex,  # Model formula
  exprObj = vst_mat_filtered,          # Expression matrix
  data = metadata_df2             # Metadata data frame
)

# Visualize variance explained by each factor
counts_varpar_plot <-plotVarPart(varPart)

ggsave("../outputs/visualisations/BE_figures/post-BE/Counts_Var_Partition.png",
       plot = counts_varpar_plot,
       width  = 6,
       height = 4,
       dpi    = 600)


# How much batch effect is influencing gene expression by identifying genes for which "Batch" explains more than 5% of variance in expression

# Identify outliers for a factor (e.g., Batch)
batch_outliers <- rownames(varPart)[varPart$Batch > 0.05]  # Outlier Genes past a 5% variance threshold

num_outliers <- length(batch_outliers)  # Number of outliers
total_genes <- nrow(vst_matrix)  # Total number of genes in the counts matrix
proportion_outliers <- (num_outliers / total_genes) * 100

print(paste("Number of proportional batch related outlier genes:",proportion_outliers))



# Identifying and quantifying genes where "Sex" explains more than 5% of the variance in expression
# Identify outliers for a factor (e.g., Sex)
sex_outliers <- rownames(varPart)[varPart$Sex > 0.05]  # Outlier Genes past a 5% variance threshold
# print(batch_outliers)


sex_num_outliers <- length(sex_outliers)  # Number of outliers
total_sex_genes <- nrow(vst_matrix)  # Total number of genes in the counts matrix
sex_proportion_outliers <- (sex_num_outliers / total_sex_genes) * 100

print(paste("Number of sex-related outlier genes (%):", sex_proportion_outliers))


# Genes that are both highly influenced by batch effects and also have a large amount of unexplained variance (residuals)
# Identify residual outliers
residual_outliers <- rownames(varPart)[varPart$Residuals > 0.8]

# Find the intersection
intersection_outliers <- intersect(batch_outliers, residual_outliers)

# Get the length of the intersection
length_of_intersection <- length(intersection_outliers)
print(paste("Number of genes in the intersection of batch and residual outliers:", length_of_intersection))


# Percentage of genes that are outliers for both Batch and Residual variance relative to the total number of genes
proportion_intersect <- (length_of_intersection / nrow(vst_matrix)) * 100

print(paste("Proportion of intersection outliers (%):", proportion_intersect))


# Performing PCA on the founded intersecting outlier batch genes and the unexplained residuals
intersect_vst_matrix <- vst_matrix[rownames(vst_matrix) %in% intersection_outliers, ]

dim(intersect_vst_matrix)

pca_intersect <- prcomp(t(intersect_vst_matrix), scale. = TRUE)

 
intersect_pca_df <- as.data.frame(pca_intersect$x) 
intersect_pca_df$Sample <- rownames(pca_intersect$x) 
 
intersect_pca_df$Batch     <- metadata_df2[intersect_pca_df$Sample, "Batch"]
intersect_pca_df$Condition <- metadata_df2[intersect_pca_df$Sample, "Condition"]
intersect_pca_df$Sex       <- metadata_df2[intersect_pca_df$Sample, "Sex"]

######## Visual if batch effects exist within intersect

intersect_pca_plot <- ggplot(intersect_pca_df,
            aes(x     = PC1,
                y     = PC2,
                shape = Condition,
                fill  = Batch,
                color = Sex)) +
  geom_point(size = 3, alpha = 0.8, stroke = 1) +
  scale_shape_manual(values = c(Relapse = 21, Remission = 24)) +
  scale_fill_manual( name= "Batch", breaks = c("Batch_1", "Batch_2"),  labels = c("Batch 1","Batch 2"), values = c("Batch_1" = "green", "Batch_2" = "orange")) +
  scale_color_manual(values = c(Male = "blue", Female = "red")) +
  guides(
    color = guide_legend(override.aes = list(shape = 21, fill = NA, stroke = 1)),
    fill  = guide_legend(override.aes = list(shape = 21, color = "black", stroke = 0.5)),
    shape = guide_legend(override.aes = list(fill = "grey", color = "black"))
  ) +
  labs(
    title = "PCA of Intersecting Genes",
    x     = "PC1",
    y     = "PC2"
  ) +
  theme_minimal()

print(intersect_pca_plot)


cat("\n======== PCA Summary Stats OF RESIDUE BATCH EFFECTS (VarPart Result PCA) ===========\n\n")
summary(pca_intersect)$importance

# Summary Stats of PCA

intersect_scree_plot <- data.frame(PC = 1:length(pca_intersect$sdev),
                         Variance = (pca_intersect$sdev)^2 / sum((pca_intersect$sdev)^2))

intersect_scree_visual <- ggplot(intersect_scree_plot, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity") +
    geom_line(aes(y = cumsum(Variance)), col = "red") +
    theme_minimal() +
    labs(title = "Batch Driven Genes PCA: Variance Explained by PCs",
         y = "Proportion of Variance Explained",
         x = "Principal Component")

ggsave(
  filename = "../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/intersect_Scree-plot.png",
  plot = intersect_scree_visual,
  width    = 6,
  height   = 4,
  dpi      = 600
)


ggsave("../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/PCA_intersect_PC1vsPC2.png",
       plot   = intersect_pca_plot,
       width  = 6,
       height = 4,
       units  = "in",
       dpi    = 300)


summary(pca_intersect)  # Check variance explained by PCs

####### Visual distribution of variance driven by batch
hist(varPart$Batch,
     breaks = 50,
     main   = "Batch explained variance across genes",
     xlab   = "Variance explained by Batch")
abline(v = c(0.05, 0.10), col = c("blue","red"), lwd = 2)

hist(varPart$Batch,
     breaks = 100,
     main   = "Batch-explained variance across genes (percentage)",
     xlab   = "Variance explained by Batch",
     xlim   = c(0, 0.2))

qs_perc <- 100 * quantile(varPart$Batch, probs = c(0.90, 0.95, 0.99))

cat(
  "Batch explained variance quantiles:\n",
  sprintf("  90th percentile: %5.1f%%  (only 10%% of genes exceed this)\n", qs_perc[1]),
  sprintf("  95th percentile: %5.1f%%  (only  5%% of genes exceed this)\n", qs_perc[2]),
  sprintf("  99th percentile: %5.1f%%  (only  1%% of genes exceed this)\n", qs_perc[3])
) 

qs <- quantile(varPart$Batch, probs = c(0.90, 0.95, 0.99))

hist(varPart$Batch,
     breaks = 50,
     main   = "Batch explained variance across genes",
     xlab   = "Variance explained by Batch (%)")
     abline(v = qs[1], col = "green", lwd = 2, lty = 2)  
     abline(v = qs[2], col = "blue",  lwd = 2, lty = 2)
     abline(v = qs[3], col = "red",   lwd = 2, lty = 2) 
     legend("topright",
       legend = c(
         paste0("90th pct = ", sprintf("%.1f%%", qs[1] * 100)),
         paste0("95th pct = ", sprintf("%.1f%%", qs[2] * 100)),
         paste0("99th pct = ", sprintf("%.1f%%", qs[3] * 100))
       ),
       col    = c("blue", "red", "darkgreen"),
       lty    = 2,
       lwd    = 2,
       bty    = "n")

# Compute the top 5% genes contributing to at least 8.3% of variance contributed to batch
cut95 <- quantile(varPart$Batch, probs = 0.95)
drop95thPercentile <- rownames(varPart)[ varPart$Batch >= cut95 ]
expr_orig    <- vst_mat_filtered

expr_95pct   <- expr_orig[ !rownames(expr_orig) %in% drop95thPercentile, ]

varPart_95 <- fitExtractVarPartModel(
  formula = ~ Batch + Condition + Sex,
  exprObj = expr_95pct,
  data    = metadata_df2
)

post95cut_plot <-plotVarPart(varPart_95)

ggsave("../../../data/outputs/visualisations/bulkRNA_prep/BE_figures/post-BE/Counts_Var_Partition_After_Advanced_BE.png",
       plot = post95cut_plot,
       width  = 6,
       height = 4,
       dpi    = 600)


# Genes will be removed in downstream analysis (deteremined as >8.3% of variance explained by batch)
writeLines(drop95thPercentile, 
           con = "../../../data/outputs/bulkRNA_processed/batch_top5pct_genes.txt")

write.csv(varPart, "../../../data/outputs/bulkRNA_processed/variancePartition_results.csv", quote = FALSE, row.names = TRUE)

# Remove all genes that had batch contributing to more than 8.3% of variance
norm_counts_df_filt <- norm_counts_df[ !rownames(norm_counts_df) %in% drop95thPercentile, , drop = FALSE]

# Outputted normalised (batch-corrected) counts as required for a GSEA input
write.table(norm_counts_df_filt, "../../../data/outputs/bulkRNA_processed/BCd+DESseq_normalized_counts.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

## Create a phenotype cls file for GSEA input
# Extract sample names and condition labels
samples <- metadata$ID
conditions <- metadata$condition

# Get unique condition names
unique_conditions <- unique(conditions)

# Create CLS file content
cls_header <- paste(length(samples), length(unique_conditions), 1, sep = " ")
cls_conditions <- paste("#", paste(unique_conditions, collapse = " "))
cls_labels <- paste(conditions, collapse = " ")

# Combine all parts
cls_content <- paste(cls_header, cls_conditions, cls_labels, sep = "\n")

# Save CLS file
writeLines(cls_content, "../../../data/outputs/bulkRNA_processed/phenotype.cls")













































