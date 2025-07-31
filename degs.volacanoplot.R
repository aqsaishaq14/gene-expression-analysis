#DEGS
#Setting data for linear models and contrasts
individual <-
  as.character(Biobase::pData(my_dataset_final)$Factor.Value..DISEASE.STATE.)
disease <- str_replace_all(Biobase::pData(my_dataset_final)$Factor.Value..DISEASE.STATE.,
                           "" , "")
disease <- ifelse(str_detect(Biobase::pData(my_dataset_final)$Factor.Value..DISEASE.STATE.,
                             "none"), "none" , "Lung cancer")
my_dataset_design <- model.matrix(~0 + disease)
colnames(my_dataset_design)[1:2] <- c("LC", "Normal")
rownames(my_dataset_design) <- individual
write.csv(my_dataset_design, "design matrix.csv")
#installing limma package for lmfit function
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")

library(limma)

#Contrasts and hypothesis tests
fit <- lmFit(my_dataset_design)
contrast_matrix <- makeContrasts('LC-Normal', levels=my_dataset_design )
#Building linear models considering contrast design
my_dataset_fit <- eBayes(contrasts.fit(lmFit(my_dataset_final,my_dataset_design), contrast_matrix))
#Extracting final results
table <- topTable(my_dataset_fit, number=Inf)
#chatgpt generated bcz my tottable obj does not found
library(limma)

# Create contrast matrix if needed (e.g., for two conditions)
contrast.matrix <- makeContrasts(Disease_vs_Control = Disease - Control, levels = design)

# Fit the linear model
fit <- lmFit(my_dataset_final, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top differentially expressed genes
top_table <- topTable(fit2, number = Inf, adjust = "fdr", sort.by = "P")

write.csv(top_table, "top_table_final.csv")


#Excluding genes with no symbols
table <- subset(table, !is.na(SYMBOL))

#Most significant differentially expressed genes
DEG_Norm_DS <- subset(table,adj.P.val < 0.05)
write.csv(DEG_Norm_DS, 'finalize_table.csv')

#Volcano Plot 
Enhancedvolcano(table,
                lab= table$SYMBOL,
                x= "LogFC", 
                y= "P.value",
                ylim= c(0, -log10(10e-12)),
                pCutoff= 0.05,
                FCcutoff=0.5,
                tile="Healthy vs Lung Cancer")
volcano_names <- ifelse(abs(my_dataset_fit$coefficients)>=0.5,
                        my_dataset_fit$genes$SYMBOL, NA)





#chatgpt generated correct code
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


