library(compositions)
library(ggplot2)
library(car) 
library(effsize)
library(vegan)
library(DirichletReg)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Cell Proportion Estimations (Deconvolution)
tpm_fractions <- read.csv("../../../data/outputs/deconvolution/cibersortx/TPM_fractions_0.5/CIBERSORTx_Adjusted.txt", sep = '\t')
# Metadata
metadata <- read.csv("../../../data/originals/bulkRNA/bulkRNA_sample_metadata/Meta_data_2.csv")

# Remove duplicates and erreanously handled samples
metadata <- metadata[!metadata$ID %in% c("2_BR10A", "2_BR13A", "BR18A"), ]

# Filtering and Prepping data:
tpm_fractions <- tpm_fractions[tpm_fractions$Mixture %in% metadata$ID, ]

tpm_fractions <- merge(tpm_fractions, metadata[, c("ID", "condition")],
                       by.x = "Mixture", by.y = "ID")

# Computing the CLR after accounting for 0 values, by adding a small non-zero values to all samples
frac.cols <- setdiff(names(tpm_fractions), c("Mixture","P.value","condition","Correlation","RMSE"))
X <- as.matrix(tpm_fractions[ , frac.cols ])

# Compute a small epsilon relative to the smallest real fraction
min_nonzero <- min(X[X > 0])
eps <- min_nonzero * 1e-3

# Add eps to every entry & renormalize rows to sum=1
X2 <- X + eps
X2 <- X2 / rowSums(X2)

X2.clr <- clr(X2)
head(X2.clr)

#=============== Uni-varate testing to see if there is a difference in individual cell composition between conditions

cond <- factor(tpm_fractions$condition,
               levels = c("Remission","Relapse"))

results <- data.frame(
  CellType = colnames(X2.clr),
  Test     = character(ncol(X2.clr)),
  p.value  = numeric(ncol(X2.clr)),
  stringsAsFactors = FALSE
)

for (i in seq_along(results$CellType)) {
  vals <- X2.clr[, i]
  # Normality per group testing
  p1 <- shapiro.test(vals[cond=="Relapse"])$p.value
  p2 <- shapiro.test(vals[cond=="Remission"])$p.value
  # Variance homogeneity testing
  lev <- leveneTest(vals ~ cond)$"Pr(>F)"[1]
  
  # Conditional test applied based on assumptions meet or failed
  if (p1 > 0.05 && p2 > 0.05) {
    if (lev > 0.05) {
      # Student t-test
      res <- t.test(vals ~ cond, var.equal = TRUE)
      results$Test[i] <- "Student t-test"
    } else {
      # Welch’s t-test
      res <- t.test(vals ~ cond)
      results$Test[i] <- "Welch t-test"
    }
  } else {
    # Wilcoxon rank-sum
    res <- wilcox.test(vals ~ cond)      
    results$Test[i] <- "Wilcoxon rank-sum"
  }
  
  results$p.value[i] <- res$p.value
}

# Ensure adjusting for multiple testing
results$p.adj <- p.adjust(results$p.value, method = "BH")

print(results)

############
clr_means <- aggregate(X2.clr, by = list(Condition = cond), FUN = mean)
# transpose & compute difference
diffs <- clr_means[clr_means$Condition=="Relapse",-1] -
         clr_means[clr_means$Condition=="Remission",-1]
effect_sizes <- data.frame(
  CellType = names(diffs),
  MeanCLR_Diff = as.numeric(diffs)
)
print(effect_sizes)

d_vals <- sapply(colnames(X2.clr), function(ct) {
  cohen.d(X2.clr[,ct], cond, hedges.correction = TRUE)$estimate
})
effect_sizes$CohensD <- d_vals

volcano_df <- data.frame(
  CellType = results$CellType,
  logP     = -log10(results$p.value),
  meanDiff = effect_sizes$MeanCLR_Diff
)
mean_clr_plot <- ggplot(volcano_df, aes(x = meanDiff, y = logP, label = CellType)) +
  geom_point() +
  geom_text(vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(x = "Mean CLR (Relapse − Remission)",
       y = "-log10(raw p-value)")

ggsave("../../../data/outputs/visualisations/cibersortx_analysis/mean_clr_plot.png",
       plot   = mean_clr_plot,
       width  = 10,
       height = 6,
       dpi    = 800)

########### Multivariate Analysis ################
distA <- dist(X2.clr)   
adonis2(distA ~ condition, data = tpm_fractions)

### Checking homogeneity of dispersion
bd <- betadisper(distA, tpm_fractions$condition)
permutest(bd)         
boxplot(bd)    

df_bd <- data.frame(
  condition = tpm_fractions$condition,
  dist2cent = bd$distances
)

bd_disp <- ggplot(df_bd, aes(x = condition, y = dist2cent, fill = condition)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = c("Relapse"   = "#FC8D62",
                               "Remission" = "#66C2A5")) +
  labs(
    title = "Homogeneity of multivariate dispersion",
    x     = NULL,
    y     = "Distance to centroid\n(CLR-Euclidean)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position     = "none",
    plot.title          = element_text(hjust = 0.5)
  )

ggsave("../../../data/outputs/visualisations/cibersortx_analysis/bd_disp.png",
       plot   = bd_disp,
       width  = 6,
       height = 4,
       dpi    = 800)

pcoa  <- cmdscale(distA, eig = TRUE, k = 2)
scores <- data.frame(PC1 = pcoa$points[,1],
                     PC2 = pcoa$points[,2],
                     condition = tpm_fractions$condition)


pcoa_scores <- ggplot(scores, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2) +
  theme_minimal() +
  labs(
    x     = "PCoA1",
    y     = "PCoA2",
    color = "Condition"        
  ) +
  guides(
    color = guide_legend(
      title.theme = element_text(size = 16, face = "bold"),
      label.theme = element_text(size = 14)
    )
  )
ggsave("../../../data/outputs/visualisations/cibersortx_analysis/pcoa_scores.png",
       plot   = pcoa_scores,
       width  = 9,
       height = 6,
       dpi    = 800)

######## Compositional Analysis using Dirichlet Regression #######
tpm_fractions_DR <- tpm_fractions
tpm_fractions_DR$comp <- DR_data(X2)

mod0 <- DirichReg(comp ~ 1,         data = tpm_fractions_DR)
#mod1 <- DirichReg(comp ~ condition, data = tpm_fractions_DR)

warns <- character()

mod1 <- withCallingHandlers(
  DirichReg(comp ~ condition, data = tpm_fractions_DR),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("=== Captured Warnings ===\n")
cat(paste0("* ", warns), sep = "\n")
cat("\n=========================\n")

# per‐component effects on the full model:
summary(mod1)

# global likelihood‐ratio test:
anova(mod0, mod1)

sm <- summary(mod1)
coef_all <- sm$coefficients


pos <- grep(":conditionRemission$", names(coef_all))

CellType <- sub(":conditionRemission$", "", names(coef_all)[pos])

Estimate <- coef_all[pos]

SE    <- sm$coef.mat[pos, "Std. Error"]
Z     <- sm$coef.mat[pos, "z value"]
P_raw <- sm$coef.mat[pos, "Pr(>|z|)"]

P_adj <- p.adjust(P_raw, method = "BH")

FC    <- exp(Estimate)
CI_lo <- exp(Estimate - 1.96*SE)
CI_hi <- exp(Estimate + 1.96*SE)

dirichlet_results <- data.frame(
  CellType, Estimate, SE, Z, P_raw, P_adj, FC, CI_lo, CI_hi,
  stringsAsFactors = FALSE
)

dirichlet_results <- dirichlet_results[order(dirichlet_results$P_raw), ]
print(dirichlet_results)

########### Visualisation of Proportional Data #############
frac.cols <- setdiff(
  names(tpm_fractions),
  c("Mixture","P.value","condition","Correlation","RMSE")
)

long <- tpm_fractions %>%
  select(Mixture, condition, all_of(frac.cols)) %>%
  pivot_longer(
    cols      = all_of(frac.cols),
    names_to  = "CellType",
    values_to = "Fraction"
  ) %>%
  mutate(
    CellType  = gsub("_", " ", CellType),  
    Mixture   = factor(Mixture, levels = unique(Mixture)),
    condition = factor(condition, levels = c("Relapse", "Remission"))
  ) %>% 
    mutate(
    CellType = recode(CellType,
      "Cycling S G2M cell" = "Cycling (S/G2–M) cell"
    )
  )

n_ct    <- n_distinct(long$CellType)
base_pal <- brewer.pal(12, "Set3")  
full_pal <- colorRampPalette(base_pal, space = "Lab")(n_ct)
anchors    <- round(seq(1, n_ct, length.out = length(base_pal)))
full_pal[anchors] <- base_pal
names(full_pal) <- sort(unique(long$CellType))
full_pal["CD4 Naive T cell"] <- "#386cb0"

#––– Creating a stacked‐bar plot
p1 <- ggplot(long, aes(x = Mixture, y = Fraction, fill = CellType)) +
  geom_col(width = 0.8) +
  facet_grid(~ condition, scales = "free_x", space = "free") +
  scale_fill_manual(values = full_pal) +
  labs(
    title = "CIBERSORTx‐inferred Cell Fractions by Sample",
    x     = NULL,
    y     = "Estimated Cell Fraction",
    fill  = "Cell Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    strip.background   = element_rect(fill = "grey90", colour = NA),
    strip.text         = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

print(p1)
ggsave("../../../data/outputs/visualisations/cibersortx_analysis/cibersortx_stacked_bar.png",
       plot   = p1,
       width  = 10,
       height = 7,
       dpi    = 800)

tpm_fractions$condition <- factor(
  tpm_fractions$condition,
  levels = c("Relapse","Remission")
)


comp_list <- lapply(
  levels(tpm_fractions$condition), 
  function(g) {
    idx      <- which(tpm_fractions$condition == g)
    X_sub    <- X2[idx, , drop = FALSE]

    X_sub_acomp <- acomp(X_sub)
    center      <- mean(X_sub_acomp)
    sd_clr      <- apply(clr(X_sub_acomp), 2, sd)

    data.frame(
      condition = g,
      CellType  = names(center),
      Center    = as.numeric(center),
      SD_CLR    = sd_clr,
      stringsAsFactors = FALSE
    )
  }
)

# 2) Bind into one data.frame
comp_stats_by_cond <- bind_rows(comp_list)

# 3) Now it’s safe to group and arrange
comp_stats_by_cond <- comp_stats_by_cond %>%
  group_by(condition) %>%
  arrange(desc(Center), .by_group = TRUE)

# 4) Pick off top‐3 and bottom‐3 within each condition
top3_by_cond    <- comp_stats_by_cond %>% slice_head(n = 3)
bottom3_by_cond <- comp_stats_by_cond %>% slice_tail(n = 3)

# 5) Render tables
cat("**Top 3 cell types by compositional center, per condition**\n\n")
print(top3_by_cond)



