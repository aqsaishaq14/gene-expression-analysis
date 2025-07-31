#DEGS

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

# Load required libraries
library(Biobase)
library(limma)
library(stringr)
library(EnhancedVolcano)

# Step 1: Define group labels
disease <- ifelse(str_detect(pData(my_dataset_final)$Factor.Value..DISEASE.STATE., "none"), 
                  "Normal", "LC")

# Step 2: Design matrix
design <- model.matrix(~0 + factor(disease))
colnames(design) <- c("LC", "Normal")
rownames(design) <- rownames(pData(my_dataset_final))

# Save the design matrix (optional)
write.csv(design, "design_matrix.csv")

# Step 3: Fit the linear model
fit <- lmFit(my_dataset_final, design)

# Step 4: Set up contrast: Lung Cancer vs Normal
contrast_matrix <- makeContrasts(LC_vs_Normal = LC - Normal, levels = design)

# Step 5: Apply contrasts and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Step 6: Get DEGs (top table)
top_table <- topTable(fit2, number = Inf, adjust = "fdr", sort.by = "P")
write.csv(top_table, "top_table_final.csv")

# Step 7: Filter DEGs with adjusted p-value < 0.05
DEG_filtered <- subset(top_table, adj.P.Val < 0.05 & !is.na(SYMBOL))
write.csv(DEG_filtered, "filtered_DEGs.csv")

# Step 8: Volcano Plot
EnhancedVolcano(DEG_filtered,
                lab = DEG_filtered$SYMBOL,
                x = 'logFC',
                y = 'P.Value',
                ylim = c(0, -log10(10e-12)),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = "Lung Cancer vs Normal")
# Add Expression_Status column to DEG_filtered
DEG_filtered$Expression_Status <- ifelse(
  DEG_filtered$logFC > 1 & DEG_filtered$adj.P.Val < 0.05, "Upregulated",
  ifelse(
    DEG_filtered$logFC < 1 & DEG_filtered$adj.P.Val < 0.05, "Downregulated",
    "Not Significant"
  )
)

# Optional: save the new table
write.csv(DEG_filtered, "new_file.csv", row.names = FALSE)


