
#overall set of call functions to implement seurat workflow

#1 package setup



# Load necessary libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(EnhancedVolcano)
# Define file paths for matrices right now I am excluding the first patient in this paper because the basic read quality metrics were way off
matrix_files <- list(
  SRR12621863 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621863_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621864 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621864_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621865 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621865_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621862 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621866_S1_Solo.out/Gene/filtered/matrix.mtx"
)

# Load data from each matrix and create Seurat objects
allCells<- lapply(names(matrix_files), function(sample) {
  matrix_path <- matrix_files[[sample]]
  matrix <- readMM(matrix_path)
  
  # Load barcodes and genes
  barcodes <- readLines(gsub("matrix.mtx", "barcodes.tsv", matrix_path))
  features <- readLines(gsub("matrix.mtx", "features.tsv", matrix_path))
  
  # Set row and column names for the matrix
  rownames(matrix) <- features
  colnames(matrix) <- barcodes
  
  # Create Seurat object
  allCells <- CreateSeuratObject(
    counts = matrix,
    project = sample,
    min.cells = 3,        # Genes must be expressed in at least 3 cells
    min.features = 200    # Cells must express at least 200 genes
  )
  #this is based on other papers, I think it's reasonable
  return(allCells)
})

# Merge Seurat objects into one
allCells <- Reduce(function(x, y) merge(x, y), allCells)


# Add mitochondrial transcript info
allCells[["percent.mt"]] <- PercentageFeatureSet(allCells, pattern = "\tMT-")



# Apply filters Here I am removing cells with less than 200 or greater than 3000 individual transcripts (should get rid of many doublets and other non-descriptive low count cells)
allCells <- subset(
  allCells,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15
)



#normalizing data (this is all the default currently): for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
allCells <- NormalizeData(allCells, normalization.method = "LogNormalize", scale.factor = 10000)



# the gene names were output in a weird tab-delimted format, fixing it for plotting
gene_symbols <- sapply(strsplit(rownames(allCells), "\t"), function(x) x[2])
rownames(allCells) <- make.unique(gene_symbols)



# let's start by just nabbing highly differentially expressed genes
allCells <- FindVariableFeatures(allCells, selection.method = "vst", nfeatures = 2000)



# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(allCells), 20)


# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(allCells)
#plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#plot1 + plot2


## ok now onto clustering

# Scale the data
allCells <- ScaleData(allCells)

#run PCA
allCells <- RunPCA(allCells, features = VariableFeatures(object = allCells))


#checking out loadings 
VizDimLoadings(allCells, dims = 1:4, reduction = "pca")

#basic PCA plot w/ 2000 most diff expressed genes
DimPlot(allCells, reduction = "pca") + NoLegend()

#now looking at cells from across the spectrum to really check out how things are divided

DimHeatmap(allCells, dims = 1:5, cells = 500,  balanced = TRUE)

#now an elbow plot to see which PCAs get fed into UMAP
ElbowPlot(allCells)


# based on elbow plot I will use PCA 1:14 

allCells <- FindNeighbors(allCells, dims = 1:14) #nearest nearest neighbor analysis
allCells <- FindClusters(allCells, resolution = .5) # find cluster, I didn't mess with the resolution much... just using what seems standard in literature

allCells <- RunUMAP(allCells, dims = 1:14, n_neighbors=30, min_dist=.3) #run unmap

DimPlot(allCells, reduction = "umap", label='true') #plot it all

allCells <- JoinLayers(allCells, assay = "RNA") # got this step from chatgpt, not totally clear on how this seurat object is 

#let's do a quick screen to confirm what our clusters are
FeaturePlot(allCells, features = c("GAD2","SLC17A6","GRIK1", "TH", "CADPS2")) 

FeaturePlot(allCells, features = c("VCAN","MOBP","AQP4","FOXJ1", "PDGFRB","CLDN5","CD74","TTR"))
FeaturePlot(allCells, features = c("CD74","CD163"))


allCells <- JoinLayers(allCells, assay = "RNA") # got this step from chatgpt, not totally clear on how this seurat object is designed, kinda weird to work w/

cluster0.markers <- FindMarkers(allCells, ident.1 = 2, ident.2 = c(13))
head(cluster0.markers, n = 50)
# what are the differences w/ this freak population, ie should we remove it?
# answer yes, top transcript is MOG, yikes, yuck, some kind of oligo thing? I will name accordingly




#with these parameters 
#clusters 0 1 3 4 and 12 are oligodendrocytes
#cluster 2 is microglia
# cluster 5 is oligodendrocyte pre-cursor cells
#cluster 6 and 9 are astrocytes (not sure what differentiates them)
#cluster 7 is endothelial cells
#cluster 8 is neuronal cells
#cluster 10 is pericytes
#cluster 11 is ependymal cells
#13 is CSF1R expressing oligodendrocytes
# 14 is gabaergic neurons







new.cluster.ids <- c("Oligo 1", "Oligo 2", "Microglia", "Oligo 3", "Oligo 4", "OPCs", "Astrocyte1","Endothelial Cells", "Excitatory Neurons", "Astrocyte2",  
                     "Pericytes", "Ependymal cells", "Oligo 5", "CSF1R expressing Oligo", "Gabaergic Neurons")
names(new.cluster.ids) <- levels(allCells)
allCells <- RenameIdents(allCells, new.cluster.ids)
DimPlot(allCells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()





#let's look at the proportion of GPNMB expressing cells


# Define a detection threshold 
gpnmb_expression <- FetchData(allCells, vars = "GPNMB")

# Convert to binary
gpnmb_detected <- gpnmb_expression > 0


# Add GPNMB detection as metadata to Seurat object
allCells$GPNMB_Expressed <- gpnmb_detected


allCells$cluster_names <- Idents(allCells)



# Fetch GPNMB expression and convert to binary (expressed or not)
gpnmb_expression <- FetchData(allCells, vars = "GPNMB")
gpnmb_detected <- gpnmb_expression > 1

# Add the binary GPNMB expression to the metadata
allCells$GPNMB_Expressed <- gpnmb_detected

# Calculate proportions based on renamed cluster identities
gpnmb_proportion <- allCells@meta.data %>%
  dplyr::group_by(cluster_names) %>%  # Use the renamed cluster names
  dplyr::summarize(Proportion = mean(GPNMB_Expressed))

# View the proportions per cluster
print(gpnmb_proportion)

library(ggplot2)
ggplot(gpnmb_proportion, aes(x = cluster_names, y = Proportion, fill = cluster_names)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of GPNMB+ Cells Across  Identified Cell Types",
       x = "Cluster",
       y = "Proportion of GPNMB+ Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#ok, now let's just subcluster our microglia and go from there

microGliaOnly <- subset(allCells, idents = "Microglia")


FeaturePlot(microGliaOnly, features = c("P2RY12","CX3CR1","TMEM119", "CD74","CD163","GPNMB","APOE","LPL","PTPRG","SPP1", "CLEC7A","TREM2"))
FeaturePlot(microGliaOnly, features = c("AIF1"))


FeaturePlot(microGliaOnly, features = c("P2RY12","CX3CR1","GPNMB","APOE",  "HSP90AA1","IL1B"))


FeaturePlot(seurat_object, features = c("GPNMB", "CD163", "SAT1", "TALAM1", "TREM2","APOE",  "CD74","P2RY12","CX3CR1"))

FeaturePlot(microGliaOnly, features = c("GPNMB", "CD163", "CD74", "ATP1A1", "ATP1A3","KCNMA1"))


#VlnPlot(seurat_object, features = c("GPNMB","P2RY12","TMEM119","IBA1", "TREM2","LPL","CLEC7A","APOE"))


FeaturePlot(microGliaOnly, features = c("VIC","GAPDH","PPIA","ACTB",  "NED", "JUN","UBC","HMBS","TBP"))


FeaturePlot(microGliaOnly, features = c("Iba1","TREM2","CXC3R1","P2RY12",  "TMEM119", "CD43","CD41","CD235a","CD44"))

FeaturePlot(microGliaOnly, features = c("ADAM10","Syndecan-4","OAZ1","GAPDH",  "ACTB", "RPL27","CST3","HEXB","CD33"))

FeaturePlot(microGliaOnly, features = c("ADAM10","Syndecan-4","OAZ1","GAPDH",  "ACTB", "RPL27","CST3","HEXB","CD33"))

FeaturePlot(microGliaOnly, features = c("ADAM10","Syndecan-4","OAZ1","GAPDH",  "ACTB", "RPL27","CST3","HEXB","CD33"))


FeaturePlot(microGliaOnly, features = c("P2RY13","PTPRC-4","PU.1","RUNX1",  "APOE", "CLEC7A","LPL","ABCA7","GPR141", "COLEC12","IL6","IL1B","TNFA"))

FeaturePlot(microGliaOnly, features= c("AIF1","RPL27"))

#re-run clustering within the microglia subset 
microGliaOnly <- FindNeighbors(microGliaOnly, dims = 1:14)
microGliaOnly <- FindClusters(microGliaOnly, resolution = 0.25)

# Visualize the clusters within microglia
DimPlot(microGliaOnly, reduction = "umap", label = TRUE)

new.cluster.ids <- c("Homeostatic", "Transition","DAM","Homeostatic2")
names(new.cluster.ids) <- levels(microGliaOnly)
microGliaOnly <- RenameIdents(microGliaOnly, new.cluster.ids)




#looks good, ok now let's look at the primary difference between putative DAM microglia and the others


DAMvsHomeo.markers <- FindMarkers(microGliaOnly, ident.1 = 2, ident.2 = c(0) )
head(DAMvsHomeo.markers, n = 100)

DAMvsMid.markers <- FindMarkers(microGliaOnly, ident.1 = 2, ident.2 = c(1) )
head(DAMvsMid.markers, n = 200)

MidvsHomeo.markers <- FindMarkers(microGliaOnly, ident.1 = 1, ident.2 = c(0) )
head(MidvsHomeo.markers, n = 100)

weirdvsHomeo.markers <- FindMarkers(microGliaOnly, ident.1 = 3, ident.2 = c(0) )
head(weirdvsHomeo.markers, n = 100)

#seperating based on top hits in each cluster 


# find markers for every cluster compared to all remaining cells, report only the positive
# ones

microGliaOnly <- ScaleData(microGliaOnly, features = rownames(microGliaOnly))

microGliaOnly.markers <- FindAllMarkers(microGliaOnly)
microGliaOnly.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# 




microGliaOnly.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(microGliaOnly, features = top10$gene) + NoLegend()

# let's also make some volcano plots

EnhancedVolcano(DAMvsHomeo.markers , 
                rownames(DAMvsHomeo.markers ),
                x ="avg_log2FC", 
                y ="p_val",
                title = 'Putative DAM cluster vs Putative Homeostatic Cluster',
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)


EnhancedVolcano(DAMvsMid.markers , 
                rownames(DAMvsMid.markers ),
                x ="avg_log2FC", 
                y ="p_val",
                title = 'Putative DAM cluster vs Putative Transition Cluster',
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)


EnhancedVolcano(MidvsHomeo.markers , 
                rownames(MidvsHomeo.markers ),
                x ="avg_log2FC", 
                y ="p_val",
                title = 'Putative Transition cluster vs Homeostatic Transition Cluster',
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)






## on to frank pseudotime, woohoo!

# Convert Seurat object to Monocle3 CDS (Cell Data Set)
#removing the homeostatic 2 cluster (it messes with trajectory analysis)
microGliaOnly <- subset(microGliaOnly, idents = c("DAM", "Homeostatic", "Transition"))



cds <- as.cell_data_set(microGliaOnly)
cds
#reorienting data correctly, thank you to bionformagician who covers this explicitly in her video https://www.youtube.com/watch?v=iq4T_uzMFcY

fData(cds)$gene_short_name <- rownames(fData(cds) ) 


#retrieve clustering information so we can integrate seurat w/ monocle3

#1 add partition
recreate.partition <-c(rep(1,length(cds@colData@rownames))) #making an empty vector (we want everything to be one partition since we are only looking at microglia)
names(recreate.partition) <- cds@colData@rownames
recreate.partition<- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <-recreate.partition
#add cluster info
list_cluster <- microGliaOnly@active.ident
cds@clusters$UMAP$clusters <-list_cluster
#add umap coordinates

cds@int_colData@listData$reducedDims$UMAP <- microGliaOnly@reductions$umap@cell.embeddings



# Learn the trajectory graph
cds <- learn_graph(cds, use_partition=FALSE)
plot_cells(cds, color_cells_by= 'cluster')

# Set root cells based on the "Homeostatic" cluster


cds_subset <- choose_cells(cds) #here I subset because there was a tiny cluster away from the primary homeostatic bunch that heavily weighted psuedotime output

cds <- order_cells(cds_subset, reduction_method= 'UMAP', root_cells = colnames(cds_subset[,clusters(cds_subset) == 'Homeostatic'] ))

# Visualize the trajectory
plot_cells(cds,
           color_cells_by ='pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves= FALSE)


#find genes that are differentially expressed across the trajectory

deg_micro <- graph_test(cds, neighbor_graph = 'principal_graph', cores=4 )



library(ggplot2)

deg_micro %>%
  arrange(q_value) %>%
  filter(status =='OK') %>%
  head(200)

# Filter for genes with Moran's I > 0.1 and rank by q-value
top_genes <- deg_micro %>%
  filter(morans_I > 0.1) %>%   # Select genes with strong spatial correlation
  arrange(q_value) %>%         # Rank by q-value (ascending order)
  head(50)                    # Select top 50 genes

# Display the table
print(top_genes)

# Create bar plot of Moran's I values for top genes
ggplot(top_genes, aes(x = reorder(gene_short_name, morans_I), y = morans_I, fill = -log10(q_value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Genes by Moran's I",
    x = "Gene",
    y = "Moran's I Value",
    fill = "-log10(q-value)"
  ) +
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient for significance
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8), 
    plot.title = element_text(hjust = 0.5)
  )


FeaturePlot(microGliaOnly, features = c("RAB42","CD163","MALAT1","CD74","TALAM1", "NEAT1","KCNMA1","GPNMB", "SLC11A1","HLA-DRA","PLCG2","PLIN2","NUPR1", "VIM","SORL1","HAMP"))


############################### this all failed for unknown reasons... sigh

#rerruning, eventually these will all be individual functions, I can sometimes get mixed. up as I just walk through the code
# Run graph test to identify significantly changing genes along the trajectory

#this is a weird bug with this function thanks to "fjrosello" for the fix https://github.com/cole-trapnell-lab/monocle3/issues/623


deg_results <- graph_test(cds, neighbor_graph = 'knn', cores = 4)

# Filter genes with significant q-values (adjusted p-value threshold)
sig_genes <- deg_results %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)  # Rank genes by significance

# View top genes


# Select the significant genes for clustering
gene_module_df <- find_gene_modules(
  cds[valid_genes, ], 
  resolution = 1e-3  
)



#########################

genes_of_interest <- c("GPNMB")  

plot_genes_in_pseudotime(cds[genes_of_interest, ], min_expr = 0.00001) +
  theme_minimal() +
  labs(title = "Expression of Selected Genes Across Pseudotime")


genes_of_interest <- c("RAB42")  

plot_genes_in_pseudotime(cds[genes_of_interest, ], min_expr = 0.000001) +
  theme_minimal() +
  labs(title = "Expression of Selected Genes Across Pseudotime")




# get and plot pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, ident, fill = ident)) +
  geom_boxplot()

# Plot GPNMB expression over pseudotime
plot_genes_in_pseudotime(cds[c("LPL", "APOE", "GPNMB","SPP1","TREM2")] )
plot_genes_in_pseudotime(cds[c("LPL", "APOE", "GPNMB","TREM2")] )

plot_genes_in_pseudotime(cds[c("RAB42", "CD74", "CD163","GPNMB")] )
plot_genes_in_pseudotime(cds[c("P2RY12", "CX3CR1", "TMEM119", "CD74")] )

plot_genes_in_pseudotime(cds["GPNMB"] )
plot_genes_in_pseudotime(cds["APOE"] )
plot_genes_in_pseudotime(cds["LPL"] )
plot_genes_in_pseudotime(cds["SPP1"] )
plot_genes_in_pseudotime(cds[c("CD74","GPNMB","RAB42","APOE")] )


plot_genes_in_pseudotime(cds[c("AIF1")] )
plot_genes_in_pseudotime(cds[c("RPL27")] )



plot_genes_in_pseudotime(cds["TREM2"] )

genes_of_interest <- c("GPNMB", "APOE", "TREM2", "CX3CR1","RAB42","CD163")  # Replace with your genes


gene_pseudotime_cor <- apply(exprs(cds[genes_of_interest, ]), 1, function(gene_expr) {
  cor(pseudotime(cds), gene_expr, method = "pearson", use = "complete.obs")
})

print(gene_pseudotime_cor)


genes_of_interest <- c("GPNMB")  # Replace with your genes

pheatmap::pheatmap(
  exprs(cds[genes_of_interest, ]),
  cluster_cols = FALSE,
  scale = "row",
  main = "Expression of Selected Genes Across Pseudotime"
)


library(ggplot2)

library(ggplot2)

pseudotime_vals <- pseudotime(cds)

for (gene in genes_of_interest) {
  data <- data.frame(
    pseudotime = pseudotime_vals,
    expression = exprs(cds[gene, ]),
    gene = gene
  )
  
  ggplot(data, aes(x = pseudotime, y = expression)) +
    geom_point() +
    geom_smooth(method = "loess") +
    ggtitle(paste("Expression of", gene, "over Pseudotime")) +
    theme_minimal()
}



# Aggregate expression data for the custom module
custom_module_expr <- aggregate_gene_expression(cds, "CustomModule")

# Plot the expression of the module over pseudotime
plot_cells(cds, color_cells_by = "CustomModule") +
  theme_minimal() +
  labs(title = "Custom Module Expression Along Pseudotime")






#let's be a bit more specific with our evals











#quick and dirty evaluation of GPNMB correlations across clusters

# Step 1: Extract Expression Data and Cluster Assignments
expr_data <- microGliaOnly[["RNA"]]@layers[["data"]]

expr_data <- expr_data[, Idents(microGliaOnly) == "Homeostatic"]


# Assign row names if missing
if (is.null(rownames(expr_data))) {
  rownames(expr_data) <- rownames(microGliaOnly[["RNA"]])
}

# Verify row names
head(rownames(expr_data))

# Extract cell cluster identities
cluster_assignments <- Idents(microGliaOnly)

# Ensure cell names match with expression data
if (is.null(names(cluster_assignments))) {
  names(cluster_assignments) <- colnames(expr_data)
}

# Step 2: Extract GPNMB Expression
# Find GPNMB row name (case insensitive)
gpnmb_row_name <- grep("GPNMB", rownames(expr_data), value = TRUE, ignore.case = TRUE)[1]

# Extract the expression data for GPNMB
gpnmb_expr <- expr_data[gpnmb_row_name, ]


# Step 3: Compute Correlations with P-Values
correlation_results <- apply(expr_data, 1, function(gene_expr) {
  cor_test <- cor.test(as.numeric(gpnmb_expr), as.numeric(gene_expr), method = "pearson", use = "complete.obs")
  c(correlation = cor_test$estimate, p_value = cor_test$p.value)
})

# Convert results to a dataframe
cor_results <- data.frame(
  gene = rownames(expr_data),
  correlation = correlation_results[1, ],
  p_value = correlation_results[2, ]
)

# Adjust p-values using the Benjamini-Hochberg method
cor_results$adjusted_p_value <- p.adjust(cor_results$p_value, method = "BH")

# Remove GPNMB from the results
cor_results <- cor_results[cor_results$gene != gpnmb_row_name, ]


# Step 4: Save and View Results
#write.csv(cor_results, "GPNMB_correlations.csv", row.names = FALSE)

# View top correlated and anticorrelated genes

# Sort by descending to get positively correlated genes
top_correlated_genes <- cor_results %>%
  filter(adjusted_p_value < 0.05) %>%
  arrange(desc(correlation)) %>%
  head(50)

# Sort by ascending to get negatively correlated genes
top_anticorrelated_genes <- cor_results %>%
  filter(adjusted_p_value < 0.05) %>%
  arrange(correlation) %>%  # Sorting by correlation in ascending order
  head(50)

# View the top results
head(top_correlated_genes, 20)
head(top_anticorrelated_genes, 20)

# Step 5: Visualization
library(ggplot2)

ggplot(top_correlated_genes, aes(x = reorder(gene, -correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top GPNMB Correlated Genes", x = "Gene", y = "Correlation") +
  theme_minimal()


ggplot(top_anticorrelated_genes, aes(x = reorder(gene, -correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top GPNMB antiCorrelated Genes", x = "Gene", y = "Correlation") +
  theme_minimal()








# Step 1: Extract Expression Data and Cluster Assignments
expr_data <- microGliaOnly[["RNA"]]@layers[["data"]]

expr_data <- expr_data[, Idents(microGliaOnly) == "Homeostatic"]

# Assign row names if missing
if (is.null(rownames(expr_data))) {
  rownames(expr_data) <- rownames(microGliaOnly[["RNA"]])
}

# Verify row names
head(rownames(expr_data))

# Extract cell cluster identities
cluster_assignments <- Idents(microGliaOnly)

# Ensure cell names match with expression data
if (is.null(names(cluster_assignments))) {
  names(cluster_assignments) <- colnames(expr_data)
}

# Step 2: Extract GPNMB Expression
# Find GPNMB row name (case insensitive)
gpnmb_row_name <- grep("GPNMB", rownames(expr_data), value = TRUE, ignore.case = TRUE)[1]

# Extract the expression data for GPNMB
gpnmb_expr <- expr_data[gpnmb_row_name, ]

# Step 3: Compute Correlations with P-Values
correlation_results <- apply(expr_data, 1, function(gene_expr) {
  cor_test <- cor.test(as.numeric(gpnmb_expr), as.numeric(gene_expr), method = "pearson", use = "complete.obs")
  c(correlation = cor_test$estimate, p_value = cor_test$p.value)
})

# Convert results to a dataframe
cor_results <- data.frame(
  gene = rownames(expr_data),
  correlation = correlation_results[1, ],
  p_value = correlation_results[2, ]
)

# Adjust p-values using the Benjamini-Hochberg method
cor_results$adjusted_p_value <- p.adjust(cor_results$p_value, method = "BH")

# Remove GPNMB from the results
cor_results <- cor_results[cor_results$gene != gpnmb_row_name, ]

# Step 4: Save and View Results
#write.csv(cor_results, "GPNMB_correlations.csv", row.names = FALSE)

# Ensure new variable names for each run
timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")

assign(paste0("top_correlated_genes_", timestamp), cor_results %>%
         filter(adjusted_p_value < 0.05) %>%
         arrange(desc(correlation)) %>%
         head(50))

assign(paste0("top_anticorrelated_genes_", timestamp), cor_results %>%
         filter(adjusted_p_value < 0.05) %>%
         arrange(correlation) %>%
         head())

head(get(paste0("top_correlated_genes_", timestamp)), 20)
head(get(paste0("top_anticorrelated_genes_", timestamp)), 20)

# Step 5: Visualization
library(ggplot2)

ggplot(get(paste0("top_correlated_genes_", timestamp)), aes(x = reorder(gene, correlation), y = correlation)) +
  geom_bar(stat = "identity",fill='darkred') +
  coord_flip() +
  labs(title = "Top GPNMB Correlated Genes Homeostatic Cells", x = "Gene", y = "Correlation") +
  theme_minimal()

ggplot(get(paste0("top_anticorrelated_genes_", timestamp)), aes(x = reorder(gene, -correlation), y = correlation)) +
  geom_bar(stat = "identity",fill='darkgreen') +
  coord_flip() +
  labs(title = "Top GPNMB antiCorrelated Genes Homeostatic Cells", x= ' ', y = "Correlation") +
  theme_minimal()
 









# Visualize QC metrics
print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", plot.cor= "false")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor= "false")
plot1 + plot2




# just for a lil sanity check let's re-image qc after these changes

# Visualize QC metrics
print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", plot.cor= "false")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor= "false")
plot1 + plot2



## there is a difference in baseline feature counts between conditions let's see if participant ID dominates variability


# Assign participant markers based on layers
seurat_object$participant_id <- NA  # Initialize participant metadata

# Loop through layers and assign  IDs
layers <- Layers(seurat_object)
participant_mapping <- c(
  "counts.SRR12621863.SeuratProject.SeuratProject" = "pt1",
  "counts.SRR12621864.SeuratProject.SeuratProject" = "pt2",
  "counts.SRR12621865.SeuratProject" = "pt3",
  "counts.SRR12621866" = "pt4"
)

# Loop through the layers and assign  IDs
for (layer_name in names(participant_mapping)) {
  # Extract cell names for the current layer
  cells_in_layer <- Cells(seurat_object, layer = layer_name)
  
  # Assign the corresponding  ID
  seurat_object$participant_id[cells_in_layer] <- participant_mapping[layer_name]
}

seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
# normally this would be 2000 but rn I'm just screening for patient effects
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Scale the data
seurat_object <- ScaleData(seurat_object)
# Run PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
DimPlot(seurat_object, reduction = "pca", group.by = "participant_id", label = TRUE)
# ok great! Participant status doesn't seem to influence clustering at all




##



seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 250)

# Create a variable feature plot
variable_feature_plot <- VariableFeaturePlot(seurat_object)

# Add labels for the top 10 variable features
variable_feature_plot <- LabelPoints(plot = variable_feature_plot, points = top10, repel = TRUE)

# Display the plot
print(variable_feature_plot)


#regress out unwanted features
seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")


## ok let's do some specific evaluations of GPNMB

# Ensure GPNMB is in your dataset
if (!"ENSG00000136235\tGPNMB\tGene Expression" %in% rownames(seurat_object)) {
  stop("tGPNMB is not found in the dataset. Check gene names.")
}


grep("GPNMB", rownames(seurat_object), value = TRUE)


# Combine normalized data from all layers
combined_data <- do.call(cbind, lapply(grep("^data", Layers(seurat_object[["RNA"]]), value = TRUE), function(layer) {
  GetAssayData(seurat_object, assay = "RNA", layer = layer)
}))



AddLayer(seurat_object, assay = "RNA", layer = "lognorm_combined", data = combined_data)




gpnmb_expression <- GetAssayData(seurat_object, assay = "RNA", layer = "lognorm_combined")["ENSG00000136235\tGPNMB\tGene Expression", ]





# Save the Seurat object for future use
saveRDS(seurat_object, file = "/Users/pennerc/Desktop/starSoloCountMatrices/seurat_object.rds")




