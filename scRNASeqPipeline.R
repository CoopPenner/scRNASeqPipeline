
#overall set of call functions to implement seurat workflow

#1 package setup




# Load necessary libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)

# Define file paths for matrices right now I am excluding the first patient in this paper because the basic read quality metrics were way off
matrix_files <- list(
  SRR12621863 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621863_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621864 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621864_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621865 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621865_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621866 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621866_S1_Solo.out/Gene/filtered/matrix.mtx"
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

allCells <- RunUMAP(allCells, dims = 1:14) #run unmap

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


FeaturePlot(microGliaOnly, features = c("P2RY12","CX3CR1", "CD74","CD163","MIF","HSP90AA1", "GPNMB","TREM2","APOE","LPL","CLEC7A"))
FeaturePlot(microGliaOnly, features = c("CD163"))


FeaturePlot(microGliaOnly, features = c("GPNMB","TREM2","APOE",  "CD74","MIF","P2RY12","CX3CR1","HSP90AA1"))


FeaturePlot(seurat_object, features = c("GPNMB", "CD163", "SAT1", "TALAM1", "TREM2","APOE",  "CD74","P2RY12","CX3CR1"))

FeaturePlot(microGliaOnly, features = c("GPNMB", "CD163", "CD74", "ATP1A1", "ATP1A3","KCNMA1"))


#VlnPlot(seurat_object, features = c("GPNMB","P2RY12","TMEM119","IBA1", "TREM2","LPL","CLEC7A","APOE"))




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

## on to frank pseudotime, woohoo!
#





library(monocle3)

# Convert Seurat object to Monocle3 CDS (Cell Data Set)
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

#plot to make sure the transfer worked

plot_cells(cds, color_cells_by='cluster',
           label_groups_by_cluster = FALSE,
           group_label_size=5) +
  theme(legend.position="right")


cluster.names<- plot_cells(cds, color_cells_by='cluster',
           label_groups_by_cluster = FALSE,
           group_label_size=5) +
  scale_color_manual(values= c('red', 'blue','green','maroon')) +
  theme(legend.position="right")



# Learn the trajectory graph
cds <- learn_graph(cds, use_partition=FALSE)
plot_cells(cds, color_cells_by= 'cluster')

# Set root cells based on the "Homeostatic" cluster
cds <- order_cells(cds, reduction_method= 'UMAP', root_cells = colnames(cds[,clusters(cds) == 2] ))

# Visualize the trajectory
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = TRUE)


# Plot GPNMB expression over pseudotime
plot_genes_in_pseudotime(cds, genes = c("GPNMB"))

# plot variable features 
plot1 <- VariableFeaturePlot(cluster3.markers)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


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




