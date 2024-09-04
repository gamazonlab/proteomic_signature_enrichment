# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)

# Define paths
to_plot_path <- "/bettimj/gamazon/hs_sig_proteomics/GeneList.txt"
file_path <- "/bettimj/gamazon/protein_expression/gtex_protein_normalized_abundance.csv"
somascan_path <- "/bettimj/gamazon/hs_sig_proteomics/somascan_genes_7k.csv"

# Read and process the to_plot file
to_plot_file <- read.table(to_plot_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
to_plot_df <- as.data.frame(to_plot_file)

# Read the main data file
data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)

# Combine the first two rows to form the new column names
tissues <- data[2, 3:ncol(data)]
col_names <- data[3, ]
data <- data[4:nrow(data), ]
colnames(data) <- col_names
gene_ids <- data$gene.id

# Process data for plotting and tau index calculation
colnames_vector <- as.vector(unlist(col_names[3:length(col_names)]))
tissues <- as.vector(unlist(tissues))
data <- data[, c(3:ncol(data))]
t_data <- transpose(data)
t_data <- as.data.frame(t_data)
t_data_merge <- cbind(colnames_vector, tissues, t_data)
colnames(t_data_merge)[3:ncol(t_data_merge)] <- gene_ids

# Remove the reference row
df <- t_data_merge
df <- df[!(df$colnames_vector == "reference"),]
df <- df[,c(2:ncol(df))]

# Reshape the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "tissues", variable.name = "gene", value.name = "expression")

# Ensure the expression column is numeric
df_long$expression <- as.numeric(df_long$expression)

# Calculate mean expression for each tissue-gene combination, excluding NA values
df_mean <- df_long %>%
  group_by(tissues, gene) %>%
  summarize(mean_expression = ifelse(all(is.na(expression)), 0, mean(expression, na.rm = TRUE)), .groups = 'drop')

# Sort genes alphabetically
df_mean <- as.data.frame(df_mean)
df_mean$gene <- as.character(df_mean$gene)
df_mean <- df_mean %>% arrange(gene)

ref_df <- fread("~/MR-JTI/model_training/JTI/gencode.v32.GRCh37.txt", header = TRUE, sep = "\t", quote = "")
ref_df <- as.data.frame(ref_df)
ref_df <- ref_df[,c("geneid", "genename")]

merged_df <- merge(df_mean, ref_df, by.x = "gene", by.y = "geneid", all.x = TRUE)
merged_df[is.na(merged_df$genename),"genename"] <- merged_df[is.na(merged_df$genename),"gene"]

somascan_file <- read.csv(somascan_path, header = TRUE, stringsAsFactors = FALSE)
somascan_df <- as.data.frame(somascan_file)

merged_df <- merged_df[(merged_df$genename %in% somascan_df$EntrezGeneSymbol),]

# Aggregate data to remove duplicates (e.g., by taking the mean of expression values if there are duplicates)
df_aggregated <- merged_df %>%
  group_by(genename, tissues) %>%
  summarise(mean_expression = mean(mean_expression, na.rm = TRUE)) %>%
  ungroup()

# Reshape the data to wide format
df_wide <- df_aggregated %>%
  pivot_wider(names_from = tissues, values_from = mean_expression, values_fill = 0)

df_wide <- as.data.frame(df_wide)
rownames(df_wide) <- df_wide$genename

df_wide <- df_wide[,2:ncol(df_wide)]

library("SummarizedExperiment")

#Opened the list of genes of interest
genes_path <- "/bettimj/gamazon/hs_sig_proteomics/GeneList.txt"
genes_file <- read.table(genes_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes_df <- as.data.frame(genes_file)
genes_df <- unique(genes_df)
write.table(genes_df, file = "/bettimj/gamazon/hs_sig_proteomics/GeneList.unique.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

tissue_metadata <- data.frame(row.names = colnames(df_wide), Tissue = colnames(df_wide))

#Tissue-specific gene retrieval
se <- SummarizedExperiment(assays = SimpleList(as.matrix(df_wide)),rowData = row.names(df_wide),colData = colnames(df_wide))
output <- teGeneRetrieval(se)
head(assay(output))

library("TissueEnrich")
library("ggplot2")
genes <- "/bettimj/gamazon/hs_sig_proteomics/GeneList.unique.txt"
inputGenes<-scan(genes, character())
gs <- GeneSet(geneIds = inputGenes)
output2 <- teEnrichmentCustom(gs, output)
enrichmentOutput<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

enrichment_df <- as.data.frame(enrichmentOutput)
enrichment_df$p <- 10^-(enrichment_df$Log10PValue)
enrichment_df$tissue <- rownames(enrichment_df)
enrichment_df <- enrichment_df[,c("tissue", "Tissue.Specific.Genes", "fold.change", "p")]

write.table(enrichment_df, file = "/bettimj/gamazon/hs_sig_proteomics/gtex_protein_enrichment_hs_genes.somascan_7k.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

library(tibble)
# Convert row names to a column called "Tissue"
#enrichmentOutput <- rownames_to_column(enrichmentOutput, var = "Tissue")

pdf("gtex_protein_enrichment_hs_genes.somascan_7k.pdf")
ggplot(enrichmentOutput, aes(x = reorder(Tissue, -Log10PValue), y = Log10PValue, 
                             label = Tissue.Specific.Genes, fill = Tissue)) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = '-LOG10(P-Adjusted)') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#Plot only those with a non-one adjusted p-value
enrichmentOutput_nonzero <- enrichmentOutput[(enrichmentOutput$Log10PValue > 0),]

pdf("gtex_protein_enrichment_hs_genes_nonzero_log_p.somascan_7k.pdf")
ggplot(enrichmentOutput_nonzero, aes(x = reorder(Tissue, -Log10PValue), y = Log10PValue, 
                             label = Tissue.Specific.Genes, fill = Tissue)) +
  geom_bar(stat = 'identity') +
  labs(x = '', y = '-LOG10(P-Adjusted)') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()