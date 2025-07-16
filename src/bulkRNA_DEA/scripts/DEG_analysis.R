#Libraries
library(DESeq2)
library(ggplot2)
library(goseq)
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(car)

# Load files
batch_corrected_counts <- read.table("../../../data/outputs/bulkRNA_processed/batch_corrected_bulkRNA_matrix.txt", sep = "\t")
metadata <- read.table("../../../data/originals/bulkRNA/bulkRNA_sample_metadata/Meta_data_2.csv", header = TRUE, sep = ",")
FC_counts <- read.table("../../../data/outputs/bulkRNA_aligned_outputs/FC_combo.txt", header = TRUE, sep = "\t")
batch_related_genes <- readLines("../../../data/outputs/bulkRNA_processed/batch_top5pct_genes.txt")


# Clean Column names
colnames(batch_corrected_counts) <- sub("^X2_", "2_", colnames(batch_corrected_counts))
print(colnames(batch_corrected_counts))
print(length(colnames(batch_corrected_counts)))

# Remove erreonous samples
metadata <- metadata[!metadata$ID %in% c("2_BR10A", "2_BR13A", "BR18A"), ]

sample_names <- colnames(batch_corrected_counts)
metadata_df <- data.frame(
  row.names = sample_names, 
  Condition = metadata$condition,  
  Sex = metadata$Sex
)

# Creating a second metadata_df2 which will be used to for normalisation
metadata_df2 <- data.frame(
  row.names = metadata$ID,  
  Batch = metadata$batch,   
  Condition = metadata$condition, 
  Sex = metadata$Sex
)

metadata_df2$Batch <- as.factor(metadata_df2$Batch)
metadata_df2$Condition <- as.factor(metadata_df2$Condition)
metadata_df2$Sex <- as.factor(metadata_df2$Sex)
str(metadata_df2)

# Check if they are now properly aligned
cat("Check if rownames of metadata_df2 match colnames of the batch corrected counts: ", all(rownames(metadata_df2) == colnames(batch_corrected_counts)), "\n")

# Remove the genes that have over 8.3% of the variance explained by Batch
batch_corrected_counts <- batch_corrected_counts[ !rownames(batch_corrected_counts) %in% batch_related_genes, , drop = FALSE]

cat("\n======== Differential Gene Analysis ========\n\n")

dds <- DESeqDataSetFromMatrix(countData = batch_corrected_counts, colData = metadata_df2, design = ~ Batch + Sex + Condition)
dds <- estimateSizeFactors(dds) 
norm_counts <- counts(dds, normalized = TRUE)  
dds <- DESeq(dds)

# What is the order of comparison?
cat("The condition levels in the DESeqDataSet are:\n")
print(levels(dds$Condition))

# DEG Results
cat("Differential expression results computed using DESeq2.\n")
res <- results(dds)

upregulated_genes   <- subset(res, log2FoldChange >  0 & padj < 0.05) # Upregulated in Remission (relative)
downregulated_genes <- subset(res, log2FoldChange <  0 & padj < 0.05) # Upregulated in Relapse (relative)

cat("Upregulated Genes in Remission: \n")
cat("Number of Upregulated Genes in Remission:", nrow(upregulated_genes), "\n")
upregulated_genes

cat("Upregulated Genes in Relapse: \n")
cat("Number of Upregulated Genes in Relapse:", nrow(downregulated_genes), "\n")
downregulated_genes

upregulated_genes_df <- as.data.frame(upregulated_genes)
downregulated_genes_df <- as.data.frame(downregulated_genes)

cat("\n======== Go Analysis (Accounting for Read Length Bias) ========\n\n")

# All of the possible gene IDs extracted to a vector
all_gene_ids <- FC_counts$Geneid

# Get all of the gene IDs from the significant upregulated DEG results 
up_significant_gene_ids <- rownames(upregulated_genes_df)
down_significant_gene_ids <- rownames(downregulated_genes_df)


# Checks to see if each gene from the total genes list is in the significant gene Id list and returns the logical vector in a integer form (0 or 1)
up_gene_vector <- as.integer(all_gene_ids %in% up_significant_gene_ids)
down_gene_vector <- as.integer(all_gene_ids %in% down_significant_gene_ids)

names(up_gene_vector) <- all_gene_ids
names(down_gene_vector) <- all_gene_ids


# Check to see if the encoding has worked
cat("Preparing data for GOseq enrichment analysis...\n")
cat("Total number of genes (background):", length(all_gene_ids), "\n")

cat("Number of significantly upregulated genes (Remission):", nrow(upregulated_genes_df), "\n")
cat("Number of significantly downregulated genes (Relapse):", nrow(downregulated_genes_df), "\n")

cat("Number of upregulated genes correctly encoded in binary vector:", length(which(up_gene_vector == 1)), "\n")
cat("Number of downregulated genes correctly encoded in binary vector:", length(which(down_gene_vector == 1)), "\n")


# Create a gene_length dataframe
gene_length_df <- FC_counts[, c(1,6)]
colnames(gene_length_df) <- c("gene_id", "length")
head(gene_length_df)

# Filter out the batch related genes lengths and ensure the index order is exactly the same as the matrix's rownames
gene_length_df <- gene_length_df[!gene_length_df$gene_id %in% batch_related_genes, , drop = FALSE]

cat(
  "Row names identical between gene length dataframe and gene matrix counts? ",
  all(gene_length_df$gene_id == rownames(batch_corrected_counts)),
  "\n"
)

# Fetching of all the go terms for all genes
gene_ids <- rownames(batch_corrected_counts)
cleaned_gene_ids <- sub("\\.\\d+$", "", gene_ids)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113)

gene2go <- getBM(
  attributes = c("ensembl_gene_id", "go_id"),
  filters = "ensembl_gene_id",
  values = cleaned_gene_ids,
  mart = ensembl
)

head(gene2go)

# Remove all of the empty searches
gene2go <- gene2go[gene2go$go_id != "", ]

# Check: Checks to see how many gene_ids had go terms from the whole bulkRNA
cat("Number of unique Ensembl gene IDs in gene2go annotation:", length(unique(gene2go$ensembl_gene_id)), "\n")

cat("Number of cleaned gene IDs in your dataset:", length(cleaned_gene_ids), "\n")

cat("Proportion of your genes with GO annotation:",
    round(length(unique(gene2go$ensembl_gene_id)) / length(cleaned_gene_ids) * 100, 2), "%\n")

# Formatting gene2go for goseq standards
gene2go_list <- gene2go %>%
  group_by(ensembl_gene_id) %>%
  summarise(go_id = list(go_id)) 

gene2go_list <- setNames(gene2go_list$go_id, gene2go_list$ensembl_gene_id)
head(gene2go_list)


# Prepping for pwf to bias longer gene lengths
# Strip the version numbers of the ensembl ID
gene_length_df$gene_id <- sub("\\.\\d+$", "", gene_length_df$gene_id)

# Get the length of the total amount of unique ensemblID names
length(unique(gene_length_df$gene_id))

# Get the gene lengths as a vector 
gene_lengths <- gene_length_df$length

# Assign the gene ids are the name for each element in the vector 
names(gene_lengths) <- gene_length_df$gene_id
head (gene_lengths[1])

# Sanity checks
sum(is.na(gene_lengths))
head(names(up_gene_vector), 10)
head(names(down_gene_vector), 10)
head(names(gene_lengths), 10)

# Remove the version numbers in gene vectors (up and down)
names(up_gene_vector) <- sub("\\.\\d+$", "", names(up_gene_vector))
names(down_gene_vector) <- sub("\\.\\d+$", "", names(down_gene_vector))

# Important: reorder gene_lengths to match a gene_vector (either up or down) order
gene_lengths <- gene_lengths[names(up_gene_vector)]

pwf_up <- nullp(up_gene_vector, bias.data = gene_lengths)
pwf_down <- nullp(down_gene_vector, bias.data = gene_lengths)


# Running Checks
cat("Checking if gene names in upregulated vector match gene_lengths vector: ",
    all(names(up_gene_vector) == names(gene_lengths)), "\n")

cat("Checking if gene names in downregulated vector match gene_lengths vector: ",
    all(names(down_gene_vector) == names(gene_lengths)), "\n")


cat("Do all upregulated genes exist in GO list?: ",
    all(names(up_gene_vector) %in% names(gene2go_list)), "\n")

cat("Do all downregulated genes exist in GO list?: ",
    all(names(down_gene_vector) %in% names(gene2go_list)), "\n")



cat("First 10 gene IDs in upregulated gene vector:\n")
print(head(names(up_gene_vector), 10))

cat("\nFirst 10 gene IDs in downregulated gene vector:\n")
print(head(names(down_gene_vector), 10))

cat("\nGO terms associated with the first 10 upregulated gene IDs:\n")
print(gene2go_list[head(names(up_gene_vector), 10)])


# Running goseq

GO_up <- goseq(pwf_up, gene2cat = gene2go_list, use_genes_without_cat = TRUE)
GO_down <- goseq(pwf_down, gene2cat = gene2go_list, use_genes_without_cat = TRUE)

GO_up$over_represented_fdr <- p.adjust(GO_up$over_represented_pvalue, method = "BH")
GO_down$over_represented_fdr <- p.adjust(GO_down$over_represented_pvalue, method = "BH")


cat ("GO Terms for upregulated genes \n")
head(GO_up)
cat ("GO Terms for downregulated genes \n")
head(GO_down)

cat("Number of enriched GO terms for upregulated genes:", nrow(GO_up), "\n")
cat("Number of enriched GO terms for downregulated genes:", nrow(GO_down), "\n")

significant_GO_up <- GO_up[GO_up$over_represented_fdr < 0.05, ]
significant_GO_down <- GO_down[GO_down$over_represented_fdr < 0.05, ]

cat ("GO Terms for significant upregulated genes \n")
head(significant_GO_up)
cat ("GO Terms for significant downregulated genes \n")
head(significant_GO_down) 


cat("\n======== Go Analysis (ClusterProfiler) ========\n\n")

# Firstly fetch all ENTREZ IDs for DEGS
sig_upreg_gene_ids <- rownames(upregulated_genes_df)
sig_downreg_gene_ids <- rownames(downregulated_genes_df)

cleaned_sig_upreg_gene_ids <- sub("\\.\\d+$", "", sig_upreg_gene_ids)
cleaned_sig_downreg_gene_ids <- sub("\\.\\d+$", "", sig_downreg_gene_ids)

length(cleaned_sig_upreg_gene_ids)
length(cleaned_sig_downreg_gene_ids)

up_deg_annotations <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = cleaned_sig_upreg_gene_ids,
  mart = ensembl
)

down_deg_annotations <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = cleaned_sig_downreg_gene_ids,
  mart = ensembl
)

# Add a cleaned gene ID column to allow for easier merging of the deg annotations to the significant_genes
rownames(upregulated_genes_df) <- sub("\\.\\d+$", "", rownames(upregulated_genes_df))
rownames(downregulated_genes_df) <- sub("\\.\\d+$", "", rownames(downregulated_genes_df))


upregulated_genes_df$ensembl_id <- rownames(upregulated_genes_df)
downregulated_genes_df$ensembl_id <- rownames(downregulated_genes_df)


merged_up_genes <- merge(upregulated_genes_df, up_deg_annotations,
                                  by.x = "ensembl_id", 
                                  by.y = "ensembl_gene_id", 
                                  all.x = TRUE)


merged_down_genes <- merge(downregulated_genes_df, down_deg_annotations,
                                  by.x = "ensembl_id", 
                                  by.y = "ensembl_gene_id", 
                                  all.x = TRUE)



# Get a ENTREZ gene list for up/downregulated genes (Remission/Relapse):
up_deg_entrez <- unique(na.omit(merged_up_genes$entrezgene_id))

down_deg_entrez <- unique(na.omit(merged_down_genes$entrezgene_id))


##Go Analysis:
# A vector of all Ensembl IDs (no version suffix)
all_ensembl <- sub("\\.\\d+$", "", rownames(batch_corrected_counts))

up_go_enrich <- enrichGO(gene = cleaned_sig_upreg_gene_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   universe = all_ensembl,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)


down_go_enrich <- enrichGO(gene = cleaned_sig_downreg_gene_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   universe = all_ensembl,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)


########### Visualiation of GO Analysis
up_go_df   <- as.data.frame(up_go_enrich@result)
down_go_df <- as.data.frame(down_go_enrich@result)

if (nrow(up_go_df) > 0) {
  n_up <- nrow(up_go_df)
  # 0.3" per term (capped at 25) + 1" padding (Dynamic resizing of plots)
  h_up <- 0.3 * min(n_up, 25) + 1

  # Barplot
  bar_up <- barplot(up_go_enrich,
                    showCategory = 20,
                    title        = "Top 20 GO Enrichment (ALL Ontologies): Genes Upregulated in Remission",
                    font.size    = 12) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 10))
  ggsave("../outputs/visualisations/GO_analysis/GO_barplot_remission.png",
         plot   = bar_up,
         width  = 10, height = h_up, units = "in", dpi = 600)

  # Dotplot
  dot_up <- dotplot(up_go_enrich,
                    showCategory = 20,
                    title        = "Top 20 GO Enrichment (ALL Ontologies): Genes Upregulated in Remission") +
            theme_classic() +
            scale_size_continuous(range = c(3, 8))
  ggsave("../outputs/visualisations/GO_analysis/GO_dotplot_remission.png",
         plot   = dot_up,
         width  = 10, height = h_up, units = "in", dpi = 600)
}

# ------------------------------
# GO Barplots & Dotplots: Up in Relapse
# ------------------------------
if (nrow(down_go_df) > 0) {
  n_dn <- nrow(down_go_df)
  # same sizing logic
  h_dn <- 0.3 * max(n_dn, 25) + 1

  # Barplot
  bar_dn <- barplot(down_go_enrich,
                    showCategory = 15,
                    title        = "ALL GO Enrichment (ALL Ontologies): Genes Upregulated in Relapse",
                    font.size    = 12) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 10))
  ggsave("../outputs/visualisations/GO_analysis/GO_barplot_relapse.png",
         plot   = bar_dn,
         width  = 10, height = h_dn, units = "in", dpi = 600)

  # Dotplot
  dot_dn <- dotplot(down_go_enrich,
                    showCategory = 20,
                    title        = "ALL GO Enrichment (ALL Ontologies): Genes Upregulated in Relapse") +
            theme_classic() +
            scale_size_continuous(range = c(3, 8))
  ggsave("../outputs/visualisations/GO_analysis/GO_dotplot_relapse.png",
         plot   = dot_dn,
         width  = 10, height = h_dn, units = "in", dpi = 600)
}



# Retrieve background ENTREZ ID for KEGG analysis (Bioconductor 3.20 -> 

# Mappings were based on data provided by: Entrez Gene ftp://ftp.ncbi.nlm.nih.gov/gene/DATA With
# a date stamp from the source of: 2024-Sep20
# Map everything in one go
map_df <- bitr(all_ensembl,
               fromType = "ENSEMBL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)
               
# How many mapped?
cat("Number of genes successfully mapped in map_df:", nrow(map_df), "\n")

# Drop NAs (bitr() drops unmapped by default)
universe_entrez <- unique(map_df$ENTREZID)

########################
## KEGG Analysis:

up_kegg_enrich <- enrichKEGG(
  gene = up_deg_entrez,
  universe = universe_entrez,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

down_kegg_enrich <- enrichKEGG(
  gene = down_deg_entrez,
  universe = universe_entrez,
  organism = "hsa",
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# View the top enriched pathways:
head(as.data.frame(up_kegg_enrich))
head(as.data.frame(down_kegg_enrich))

up_kegg_df <- as.data.frame(up_kegg_enrich)
down_kegg_df <- as.data.frame(down_kegg_enrich)


# View top enriched pathways
cat("Top KEGG pathways for upregulated genes:\n")
print(head(up_kegg_df))

cat("\nTop KEGG pathways for downregulated genes:\n")
print(head(down_kegg_df))

# Get the vector of log fold changes in up/down regulated DEGs
# Upregulated FC vector (named by Entrez ID)
fc_up <- merged_up_genes$log2FoldChange
names(fc_up) <- merged_up_genes$entrezgene_id

# Downregulated FC vector (named by Entrez ID)
merged_down_genes_cleaned <- merged_down_genes[!is.na(merged_down_genes$entrezgene_id), ]
fc_down <- merged_down_genes_cleaned$log2FoldChange
names(fc_down) <- merged_down_genes_cleaned$entrezgene_id

################ Visualisation of KEGG Enrichment Analysis

# =========================================
# Upregulated in Remission
# =========================================
if (nrow(up_kegg_df) > 0) {
  # calculate height once
  n_up <- nrow(up_kegg_df)
  h_up <- 0.4 * min(n_up, 30) + 2
  
  # Dotplot
  dot_up <- dotplot(up_kegg_enrich, showCategory = 20) +
            ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Remission")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_dotplot_remission.png",
         plot   = dot_up,
         width  = 10, height = h_up, units = "in", dpi = 600)
  
  # Barplot
  bar_up <- barplot(up_kegg_enrich, showCategory = 20) +
            ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Remission")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_barplot_remission.png",
         plot   = bar_up,
         width  = 10, height = h_up, units = "in")
  
  # Cnetplot
  cnet_up <- cnetplot(up_kegg_enrich, showCategory = 20, foldChange = fc_up) +
             ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Remission") +
             labs(caption = "229 = ALDOB\n80201 = HKDC1") +
             theme_minimal(base_size = 12) +
             theme(
               plot.caption          = element_text(hjust = 0, size = 10, margin = margin(t = 10)),
               plot.caption.position = "plot"
             )
  ggsave(".../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_cnetplot_remission.png",
         plot   = cnet_up,
         width  = 10, height = h_up, units = "in", dpi = 600)
  
  # Heatplot
  heat_up <- heatplot(up_kegg_enrich, foldChange = fc_up, showCategory = 20) +
             ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Remission")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_heatplot_remission.png",
         plot   = heat_up,
         width  = 10, height = h_up, units = "in", dpi = 600)
}

# =========================================
# Upregulated in Relapse
# =========================================
if (nrow(down_kegg_df) > 0) {
  # calculate height once
  n_dn <- nrow(down_kegg_df)
  h_dn <- 0.4 * min(n_dn, 30) + 2
  
  # Dotplot
  dot_dn <- dotplot(down_kegg_enrich, showCategory = 15) +
            ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Relapse")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_dotplot_relapse.png",
         plot   = dot_dn,
         width  = 10, height = h_dn, units = "in", dpi = 600)
  
  # Barplot
  bar_dn <- barplot(down_kegg_enrich, showCategory = 20) +
            ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Relapse")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_barplot_relapse.png",
         plot   = bar_dn,
         width  = 10, height = h_dn, units = "in")
  
  # Cnetplot
  cnet_dn <- cnetplot(down_kegg_enrich, showCategory = 20, foldChange = fc_down) +
              ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Relapse") + 
             ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_cnetplot_relapse.png",
         plot   = cnet_dn,
         width  = 10, height = h_dn, units = "in", dpi = 600)
  
  # Heatplot
  heat_dn <- heatplot(down_kegg_enrich, foldChange = fc_down, showCategory = 20) +
             ggtitle("KEGG Pathway Enrichment: Genes Upregulated in Relapse")
  ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/KEGG_analysis/KEGG_heatplot_relapse.png",
         plot   = heat_dn,
         width  = 10, height = h_dn, units = "in", dpi = 600)
}


########### Biomarker plotting

# #### DESeq2 Plot Generation for ENHO

plotCounts(dds,
           gene    = "ENSG00000168913.7",
           intgroup = "Condition",
           returnData = FALSE)    # set TRUE if you just want the data frame

############ Manual Plot Generation of ENHO
# Compute the VST matrix (optional but to provide two normalised datasets)
vsd <- vst(dds)
vst_mat <- assay(vsd)


analyze_gene_set <- function(genes, norm_mat, vst_mat, metadata, outdir, width   = 6, height  = 5, dpi = 600) {
  
  # strip version suffix helper
  strip_version <- function(id) sub("\\.\\d+$", "", id)

  base_ids <- unique(strip_version(genes))
  gene_map <- getBM(
    attributes = c("ensembl_gene_id","hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = base_ids,
    mart       = ensembl
  )
  
  gene_map$base      <- strip_version(gene_map$ensembl_gene_id)
  sym_lookup         <- setNames(gene_map$hgnc_symbol, gene_map$base)
  
  for (gene in genes) {
    base    <- strip_version(gene)
    # find the exact row name in norm_mat
    pattern <- paste0("^", base, "(\\.[0-9]+)?$")
    matches <- grep(pattern, rownames(norm_mat), value = TRUE)
    if (length(matches) == 0) {
      warning("Skipping ", gene, ": not found in norm_mat")
      next
    }
    gene_key <- matches[1]
    gene_sym <- sym_lookup[base]
    if (is.na(gene_sym) || gene_sym == "") gene_sym <- base
    

    # build data.frame for a given matrix
    build_df <- function(mat) {
      df <- data.frame(
        Sample     = colnames(mat),
        Expression = mat[gene_key, ],
        stringsAsFactors = FALSE
      )
      df <- merge(df,
                  metadata[, c("ID", "condition", "batch", "Sex")],
                  by.x = "Sample", by.y = "ID")
      names(df)[names(df) == "condition"] <- "Condition"
      df$Condition <- factor(df$Condition,
                             levels = c("Relapse", "Remission"))
      return(df)
    }

    df_norm <- build_df(norm_mat)
    df_vst  <- build_df(vst_mat)

    
    run_and_plot <- function(df, suffix) {
      # Statistical testing assumptions for t-tests parametric and non-parametric
      sw1 <- shapiro.test(df$Expression[df$Condition == "Relapse"])$p.value
      sw2 <- shapiro.test(df$Expression[df$Condition == "Remission"])$p.value
      lev <- leveneTest(Expression ~ Condition, data = df)[["Pr(>F)"]][1]

      # Choose statistical test based on assumptions
      if (sw1 > 0.05 && sw2 > 0.05 && lev > 0.05) {
        tst   <- t.test(Expression ~ Condition, data = df, var.equal = TRUE)
        test_nm <- "Student's t-test"
      } else if (sw1 > 0.05 && sw2 > 0.05) {
        tst   <- t.test(Expression ~ Condition, data = df, var.equal = FALSE)
        test_nm <- "Welch's t-test"
      } else {
        tst   <- wilcox.test(Expression ~ Condition, data = df)
        test_nm <- "Wilcoxon rank-sum"
      }

      pval <- signif(tst$p.value, 3)
      ci   <- if (!is.null(tst$conf.int))
                sprintf("CI %.2f-%.2f", tst$conf.int[1], tst$conf.int[2])
              else NA

      message(sprintf("%s [%s]: %s p=%s %s", gene_sym, suffix, test_nm, pval, ci))

      csv_file <- file.path(outdir, sprintf("%s_%s_data.csv", gene_sym, suffix))
      write.csv(df, csv_file, row.names = FALSE)

      # Visualisation of Data
      p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
        scale_y_log10() +
        labs(
          title    = sprintf("Expression of %s: %s", gene_sym, gene_key),
          subtitle = sprintf("%s p = %s %s", test_nm, pval, ci),
          x        = NULL,
          y        = ifelse(suffix == "vst",
                            "VST counts (log10)",
                            "Norm counts (log10)")
        ) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5)
        )
      
      png_file <- file.path(outdir, sprintf("%s_%s_plot.png", gene_sym, suffix))
      ggsave(filename = png_file,
             plot     = p,
             width    = width,
             height   = height,
             units    = "in",
             dpi      = dpi,
             bg       = "white")
    }

    # run for both normalized and VST
    run_and_plot(df_norm, suffix = "norm")
    run_and_plot(df_vst,  suffix = "vst")
  }
  invisible(NULL)
}

analyze_gene_set(
  genes    = c("ENSG00000168913.7", "ENSG00000164379.7"),
  norm_mat = norm_counts,
  vst_mat  = vst_mat,
  metadata = metadata,
  outdir   = "../../../data/outputs/visualisations/bulkRNA_DEA/biomarkers",
  width    = 6,
  height   = 5,
  dpi      = 300
)

############# DEG Volano plot

res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
res_df$ensembl_id <- sub("\\.\\d+$", "", res_df$ensembl_id)

# Drop NAs (which are genes with too-low counts)
deg_clean <- res_df[!is.na(res_df$padj), ]

# Mapping to HGNC symbols 
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = deg_clean$ensembl_id,
  mart       = ensembl
)
deg_annot <- merge(
  deg_clean, gene_map,
  by.x = "ensembl_id", by.y = "ensembl_gene_id",
  all.x = TRUE
)

# Flagging significance of genes
lfc_cut <- 0.5
padj_cut <- 0.05
deg_annot$Significance <- "Not significant"
deg_annot$Significance[
  deg_annot$padj < padj_cut & deg_annot$log2FoldChange >  lfc_cut
] <- "Upregulated in Remission"
deg_annot$Significance[
  deg_annot$padj < padj_cut & deg_annot$log2FoldChange < -lfc_cut
] <- "Upregulated in Relapse"

# Compute –log10(padj) and cap outliers
deg_annot$negLog10Padj <- -log10(deg_annot$padj)
deg_annot$negLog10Padj_capped <- pmin(deg_annot$negLog10Padj, 20)

# Visalisation using ggplot2
p_volcano <- ggplot(deg_annot,
                    aes(x = log2FoldChange,
                        y = negLog10Padj_capped,
                        color = Significance)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut),
             linetype = "dashed") +
  geom_text_repel(
    data = subset(deg_annot, Significance != "Not significant"),
    aes(label = hgnc_symbol),
    size         = 2.8,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.alpha= 0.5
  ) +
  scale_color_manual(values = c(
    "Not significant" = "grey65",
    "Upregulated in Remission"     = "dodgerblue",
    "Upregulated in Relapse"   = "firebrick"
  )) +
  labs(
    title    = "Volcano Plot of Differential Gene Expression",
    subtitle = paste0("Cutoffs: |log₂FC| > ", lfc_cut, 
                      ", padj < ", padj_cut),
    x     = "log2(Fold Change)",
    y     = "-log10(Adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  coord_cartesian(
    ylim = c(0, max(deg_annot$negLog10Padj_capped, na.rm = TRUE) + 1)
  )

print(p_volcano)

ggsave("../../../data/outputs/visualisations/bulkRNA_DEA/DEG/volcano_plot.png", plot = p_volcano,
       width = 8, height = 6, dpi = 600)


























 















