#packages
install.packages("data.table")
install.packages('R.utils')
install.packages('tidyverse')
library(data.table)
library(tidyverse)
#reading the tables
umi_matrix <- fread('Ovarian_Cancer_raw_UMI_matrix.csv', sep='\t')
barcode_table <- read_tsv('Ovarian_Cancer_cell_annotation.tsv')

# Make a compatible umi table
row_names <- umi_matrix$Gene_symbol
#umi_matrix <- umi_matrix[, -1]  
rownames(umi_matrix) <- row_names
umi_matrix <- as.data.frame(umi_matrix)

#finding duplicates
barcode_duplicates <- barcode_table[duplicated(barcode_table$Index), ]
first_duplicates <- barcode_duplicates[,1]
first_duplicates <- as.list(first_duplicates$Index)

second_duplicates <- colnames(umi_matrix)[!colnames(umi_matrix) %in% barcode_table$Index]
second_duplicates <- as.list(second_duplicates)
duplicates <- c(first_duplicates, second_duplicates)

#excluding rows that we don't know where do they come from, we have duplicates in barcode table, while in umi matrix they are renamed 

umi_matrix <- umi_matrix[, !colnames(umi_matrix) %in% duplicates]

barcode_table <- barcode_table[!(barcode_table$Index %in% duplicates), ]

#excluding patients with early stages

umi_matrix$Index <- rownames(umi_matrix)
umi_wp <- merge.data.frame(barcode_table[,c(1,2)],umi_matrix, by.x = 'Index')
umi_wp <- umi_wp[umi_wp$Patient %in% c('ovCHA004','ovCHA034','ovCHA066', 'ovCHA070', 'ovCHA107', 'ovCHA110'), ]
umi_fs <- as.matrix(t(umi_wp[,-2]))

#making sparse matrix


colnames <- as.character(umi_fs[1,])
colnames
umi_data <- umi_fs[-1,]
umi_data[1,1]
rownames <- rownames(umi_data)
rownames
saveRDS(umi_data, file = "umi_data.rds")
umi_data <- readRDS("umi_data.rds")
umi_data[1:10,1:10]
saveRDS(umi_data, file = "umi_data.rds")


#cluster creation and converting matrix to numeric
library(parallel)
detectCores()
num_cores <- min(detectCores(), 8)  
cl <- makeCluster(num_cores)
print(cl)
umi_data <- as.matrix(parApply(cl, umi_data, 2, as.numeric))
stopCluster(cl)
str(umi_data)
str(umi_data[1:10,1:10])

#matrix sparsing
install.packages("Matrix")
library(Matrix)
sparse_umi_matrix <- Matrix(umi_data, sparse = TRUE)
rownames(sparse_umi_matrix) <- rownames
colnames(sparse_umi_matrix) <- colnames
saveRDS(sparse_umi_matrix, file = "sparse_umi_matrix.rds")
barcode_table <- as.data.frame(barcode_table)
rownames(barcode_table) <- barcode_table$Index
barcode_table <- barcode_table[,-1]
saveRDS(barcode_table, file = "barcode_table.rds")
barcode_table <- readRDS('barcode_table.rds')
sparse_umi_matrix <- readRDS('sparse_umi_matrix.rds')

#reading lung cancer table
lung_patients <- fread('meta.csv')
lung_patients <- lung_patients[lung_patients$PatientID %in% c('PS08','PS09', 'PA20', 'PA15', 'PA18', 'PA19'),]
rownames(lung_patients) <- lung_patients$Cells
saveRDS(lung_patients, file = "lung_patients.rds")
lung_patients <- readRDS('lung_patients.rds')
lung_barcodes <- c(lung_patients$Cells, "V1") #as V1 is a column name for genes in the assay
sc_lung <- fread('raw_data.csv', select = lung_barcodes, sep= ',', quote = "\"" )

sc_lung_rownames <- sc_lung$V1
sc_lung <- sc_lung[, -23185, with = FALSE]
rownames(sc_lung) <- sc_lung_rownames
saveRDS(sc_lung, 'sc_lung.rds')
sc_lung_matrix <- as.matrix(sc_lung)
sparse_sc_lung <- Matrix(sc_lung_matrix, sparse = TRUE)
sparse_sc_lung <- t(sparse_sc_lung)
saveRDS(sparse_sc_lung, file = "sparse_sc_lung.rds")
sparse_sc_lung <- readRDS('sparse_sc_lung.rds')


# Excluding non gene symbol rows:
classify_id <- function(id) {
  grepl("^[A-Z0-9]+(-[A-Z0-9]+)*$", id)
}

gene_ids_ovarian <- rownames(sparse_umi_matrix)
gene_ids_lung <- rownames(sparse_sc_lung)
valid_ids_ovarian <- sapply(gene_ids_ovarian, classify_id)
valid_ids_lung <- sapply(gene_ids_lung, classify_id)
ovarian_filtered <- sparse_umi_matrix[valid_ids_ovarian, , drop = FALSE]
lung_filtered <- sparse_sc_lung[valid_ids_lung, , drop = FALSE]
ovarian_id_types_filtered <- sapply(rownames(ovarian_filtered), classify_id)
table(ovarian_id_types_filtered)

#Seurat 
# packages
install.packages("devtools")
install.packages("Seurat")
install.packages('BiocManager')
BiocManager::install('glmGamPoi')

library(dplyr)
library(Seurat)
library(patchwork)

#creating the Seurat objects

ovarian <- CreateSeuratObject(counts = ovarian_filtered, min.cells = 3, min.features = 200, project = 'ovarian')
lung <- CreateSeuratObject(counts = lung_filtered, min.cells = 3, min.features = 200, project = 'lung')
saveRDS(ovarian, file = "ovarian.rds")



# Cutoff setting
ovarian[["percent.mt"]] <- PercentageFeatureSet(ovarian, pattern = "^MT-")
VlnPlot(ovarian, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ovarian <- subset(ovarian, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)

lung[["percent.mt"]] <- PercentageFeatureSet(lung, pattern = "^MT-")
VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
lung <- subset(lung, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

# Add metadata
ovarian[["patient"]] <- barcode_table$Patient
lung[["patient"]] <- lung_patients$PatientID
ovarian$dataset <- 'ovarian'
lung$dataset <- 'lung'

# Integration
options(future.globals.maxSize = 3000 * 1024^2)  # Set to 3000 MiB

lung_ovarian_combined <- merge(ovarian, lung, add.cell.ids = c("ovarian", "lung"))
lung_ovarian_combined <- SCTransform(lung_ovarian_combined, vst.flavor = "v2", conserve.memory = TRUE)
lung_ovarian_combined <- RunPCA(lung_ovarian_combined)
lung_ovarian_combined <- IntegrateLayers(
  object = lung_ovarian_combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT",
  verbose = FALSE
)
lung_ovarian_combined <- FindNeighbors(lung_ovarian_combined, reduction = "integrated.rpca", dims = 1:30)
lung_ovarian_combined <- FindClusters(lung_ovarian_combined, resolution = 2, cluster.name = "rpca_clusters")
lung_ovarian_combined <- RunUMAP(lung_ovarian_combined, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
nk_seurat <- subset(lung_ovarian_combined, lung_ovarian_annotations == "NK cells")

DimPlot(lung_ovarian_combined, reduction = "umap.rpca", group.by = "lung_ovarian_annotations", label = TRUE, label.size = 3)+
  theme(legend.position = "none")

#Cluster profiler
BiocManager::install("clusterProfiler")
library(clusterProfiler)
keyType <- "SYMBOL"
gene_list <- c('KLRC3', 'CD3A', 'CD160', 'KLRF1', 'CD244 ', 'KLRC1 ', 'GNLY', 'PRF1', 'XCL2', 'IL18RAP')
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = keyType,
                ont = "BP",  # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE) 
#BlueprintEncodeData
install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
library(celldex)
library(SingleR)
reference <- celldex::BlueprintEncodeData()
reference_expression <- assay(reference, "logcounts")
lung_ovarian_counts <- GetAssayData(lung_ovarian_combined, layer = "counts")

main_labels <- reference@colData$label.main
# Perform SingleR cell annotation
lung_ovarian_annotations <- SingleR(
  test = lung_ovarian_counts, 
  ref = reference_expression,
  labels = main_labels
)
lung_ovarian_annotations
str(reference)
lung_ovarian_combined[['lung_ovarian_annotations']] <- lung_ovarian_annotations$labels
confidence_scores <- lung_ovarian_annotations$scores
hist(confidence_scores, breaks = 50, main = "Distribution of Confidence Scores", xlab = "Confidence Score")


DimPlot(lung_ovarian_combined, reduction = "umap.rpca", group.by = "rpca_clusters")
DimPlot(nk_seurat, reduction = "umap.rpca", split.by = "dataset")
FeaturePlot(lung_ovarian_combined, features = c("NKG7"), min.cutoff = "q10", max.cutoff = "q90",
            cols = c("lightblue", "darkblue"), label = TRUE, pt.size = 0.5) +
  ggplot2::ggtitle("Expression of NKG7 on UMAP") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
nk_cells_based_on_clusters <- subset(lung_ovarian_combined, subset = rpca_clusters %in% c('11', '27', '26'))
DimPlot(nk_cells_based_on_clusters, reduction = "umap.rpca", split.by = "dataset", group.by = 'lung_ovarian_annotations')
FeaturePlot(nk_cells_based_on_clusters, features = c("NKG7"), min.cutoff = "q10", max.cutoff = "q90",
            cols = c("lightblue", "darkblue"), label = TRUE, pt.size = 0.5) +
  ggplot2::ggtitle("Expression of NKG7 on UMAP") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
# Count NK cells
ovarian_cells <- lung_ovarian_combined[, lung_ovarian_combined$dataset == 'ovarian']
total_cells_ovarian <- ncol(ovarian_cells)
total_cells_ovarian
nk_cell_count <- sum(ovarian_cells$lung_ovarian_annotations == "NK cells")
nk_cell_count

lung_cells <- lung_ovarian_combined[, lung_ovarian_combined$dataset == 'lung']
total_cells_lung <- ncol(lung_cells)
total_cells_lung
nk_cell_lung <- sum(lung_cells$lung_ovarian_annotations == "NK cells")
nk_cell_lung

##1st DE for all cells
#Pseudobulk DE 
pseudobulk <- AggregateExpression(lung_ovarian_combined, assays = "RNA", return.seurat = T, group.by = c( "patient", "lung_ovarian_annotations", "dataset"))
tail(Cells(pseudobulk))
pseudobulk$celltype.dataset <- paste(pseudobulk$lung_ovarian_annotations, pseudobulk$dataset, sep = "_")
Idents(pseudobulk) <- "celltype.dataset"
nk_subset <- subset(pseudobulk, idents = c("NK cells_ovarian", "NK cells_lung"))
table(Idents(nk_subset))
table(pseudobulk$celltype.dataset)

# Now perform the DE analysis on this NK cell subset
bulk.nk.de <- FindMarkers(nk_subset, 
                          ident.1 = "NK cells_ovarian", 
                          ident.2 = "NK cells_lung",
                          test.use = "DESeq2")

head(bulk.nk.de, n = 50)


#Heatmap
library(ComplexHeatmap)
library(circlize)

bulk.nk.de.filtered <- bulk.nk.de %>%
  filter(p_val_adj <= 0.05)

#Top genes based on adjusted p-values
top_genes <- rownames(head(bulk.nk.de.filtered[order(bulk.nk.de.filtered$avg_log2FC), ], 40))

# Prepare data for the heatmap

data_for_heatmap <- GetAssayData(nk_subset, layer = "data")[top_genes, ]

# Create the heatmap
Heatmap(as.matrix(data_for_heatmap),
        name = "expression",
        column_title = "Expression of Top 40 DE Genes in NK cells",
        row_names_side = "left",
        cluster_rows = TRUE,
        cluster_columns = TRUE
        )

#Volcanoplot
library(ggplot2)
library(ggrepel)
# Sorting by absolute log2 fold change to find the top 30 genes
top_genes <- bulk.nk.de %>%
  mutate(gene = rownames(bulk.nk.de)) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice(1:30)

# Create a volcano plot with gene names annotated
volcano_plot <- ggplot(bulk.nk.de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj <= 0.05), alpha = 0.4) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3,
                  box.padding = 0.35, point.padding = 0.5,
                  max.overlaps = 50) +
  labs(x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Volcano Plot with Top 30 Genes Labeled") +
  theme_minimal() +
  theme(legend.position = "none")

# Print the plot without specified xlim and ylim
print(volcano_plot)


##2nd DE for NK subset only
#Pseudobulk DE 
pseudobulk_nk <- AggregateExpression(nk_cells_based_on_clusters, assays = "RNA", return.seurat = T, group.by = c( "patient", "dataset"))
tail(Cells(pseudobulk_nk))
Idents(pseudobulk_nk) <- "dataset"
# Now perform the DE analysis on this NK cell subset
bulk.nk.de_2 <- FindMarkers(pseudobulk_nk, 
                          ident.1 = "ovarian", 
                          ident.2 = "lung",
                          test.use = "DESeq2")

head(bulk.nk.de_2, n = 50)


#Heatmap
bulk.nk.de_2.filtered <- bulk.nk.de_2 %>%
  filter(p_val_adj <= 0.05)

#Top genes based on adjusted p-values
top_genes_2 <- rownames(head(bulk.nk.de_2.filtered[order(bulk.nk.de_2.filtered$avg_log2FC), ], 40))

# Prepare data for the heatmap

data_for_heatmap_2 <- GetAssayData(pseudobulk_nk, layer = "data")[top_genes_2, ]

# Create the heatmap
Heatmap(as.matrix(data_for_heatmap_2),
        name = "expression",
        column_title = "Expression of Top 40 DE Genes in NK cells",
        row_names_side = "left",
        cluster_rows = TRUE,
        cluster_columns = TRUE
)

#Volcanoplot
# Sorting by absolute log2 fold change to find the top 30 genes
top_genes_2 <- bulk.nk.de_2 %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice(1:30) %>%
  mutate(gene = rownames(.))
# Create a volcano plot with gene names annotated
volcano_plot_2 <- ggplot(bulk.nk.de_2, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj <= 0.05), alpha = 0.4) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(data = top_genes_2, aes(label = gene), size = 3,
                  box.padding = 0.35, point.padding = 0.5,
                  max.overlaps = 50) +
  labs(x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Volcano Plot with Top 30 Genes Labeled") +
  theme_minimal() +
  theme(legend.position = "none")

# Print the plot without specified xlim and ylim
print(volcano_plot_2)
