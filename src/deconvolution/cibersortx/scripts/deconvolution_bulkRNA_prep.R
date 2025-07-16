library(sva)

raw_counts <-read.table("../../../data/outputs/bulkRNA_aligned_outputs/FC_combo.txt", header = TRUE, sep = "\t")
metadata <- read.table("../../../data/originals/bulkRNA_sample_metadata/Meta_data_2.csv", header = TRUE, sep = ",")

# Clean the column names so that only the sample names appear (previously was just the full directory from the HPC)

clean_colnames <- function(x) {
  gsub(".*disease\\.(.*?)_Aligned.*", "\\1", x)
}

colnames(raw_counts)[7:ncol(raw_counts)] <- clean_colnames(colnames(raw_counts)[7:ncol(raw_counts)])

print(colnames(raw_counts))


counts <- raw_counts[, 7:ncol(raw_counts)]  # Adjust column range as needed
rownames(counts) <- raw_counts$Geneid  # Set gene IDs as row names

# Remove duplicate samples and erroneous before batch correction
counts <- counts[, !colnames(counts) %in% c("2_BR10A", "2_BR13A", "BR18A")]
metadata <- metadata[!metadata$ID %in% c("2_BR10A", "2_BR13A", "BR18A"), ]

sample_names <- colnames(counts)  # Get column names from count matrix

metadata_df <- data.frame(
  row.names = sample_names,  # Assign sample IDs as row names
  Condition = metadata$condition,  # Extract condition column
  Sex = metadata$Sex
)

# Extract batch and condition information from metadata
batch <- as.factor(metadata$batch)
condition <- as.factor(metadata$condition)

# Apply ComBat-Seq
adjusted_counts <- ComBat_seq(counts = as.matrix(counts), batch = batch, group = condition)

# Function creates a df suitable for cibersortx format, Makes a column called Genes and strips Ensembl ID version numbers
CibersortX_Structure_Formatter <- function(mat_or_df) {
  df <- as.data.frame(mat_or_df)
  df$Gene <- sub("\\..*", "", rownames(df))  # Strip version from Ensembl ID
  df <- df[, c("Gene", setdiff(names(df), "Gene"))]  # Move Gene column to front
  return(df)
}


# ======== Creating a TPM Matrix ======== #

# Subset only the gene names and the gene lengths
gene_lengths_df <- raw_counts[, c(1, 6)]

# Check to see if the gene ID names are identical (the same in order and content)
if (all(rownames(adjusted_counts) == gene_lengths_df$GeneID)) {
  cat("Gene IDs in counts and gene length table are identical and aligned.\n")
} else {
  cat("Gene IDs are not perfectly aligned, ensure you merge and realign.\n")
}


## TPM Computation
gene_lengths <- gene_lengths_df$Length
length_kb <- gene_lengths / 1000
rpk_matrix <- sweep(adjusted_counts, 1, length_kb, FUN = "/")
tpm_matrix <- sweep(rpk_matrix, 2, colSums(rpk_matrix), FUN = "/") * 1e6

# Format the data structure for cibersortx input
tpm_df <- CibersortX_Structure_Formatter(tpm_matrix)

print("Are the gene Ids unique?")
length(unique(tpm_df$Gene)) == nrow(tpm_df)

cat("\n==== Outputting TPM, CPM and Raw Count bulkRNA matrix versions ======\n")

write.table(tpm_df, file = "../../../data/outputs/bulkRNA_deconvolution_prep/corrected_TPM_counts_bulkRNA_CSX.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ======== Creating a CPM Matrix ======== #
cpm_matrix <- sweep(adjusted_counts, 2, colSums(adjusted_counts), FUN = "/") * 1e6
cpm_df <- CibersortX_Structure_Formatter(cpm_matrix)

write.table(cpm_df, file = "../../../data/outputs/bulkRNA_deconvolution_prep/corrected_CPM_counts_bulkRNA_CSX.txt",sep = "\t", quote = FALSE, row.names = FALSE)


# ======== Writing a Raw Count Version ======== #
raw_counts <- adjusted_counts
raw_counts_df <- CibersortX_Structure_Formatter(raw_counts)

write.table(raw_counts_df,file = "../../../data/outputs/bulkRNA_deconvolution_prep/corrected_raw_counts_bulkRNA_CSX.txt", sep = "\t", quote = FALSE, row.names = FALSE)
