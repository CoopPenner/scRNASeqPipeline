
#overall set of call functions to implement seurat workflow

#1 get all relevant pacakges 


# Install if not already installed
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
  install.packages("leiden")
  install.packages("monocle3")
}



# Load necessary libraries
library(Seurat)
library(Matrix)
library(dplyr)

# Define file paths for matrices right now I am excluding the first patient in this paper because the basic read quality metrics were way off
matrix_files <- list(
  SRR12621863 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621863_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621864 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621864_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621865 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621865_S1_Solo.out/Gene/filtered/matrix.mtx",
  SRR12621866 = "/Users/pennerc/Desktop/starSoloCountMatrices/SRR12621866_S1_Solo.out/Gene/filtered/matrix.mtx"
)

# Load data from each matrix and create Seurat objects
seurat_objects <- lapply(names(matrix_files), function(sample) {
  # Load the matrix
  matrix_path <- matrix_files[[sample]]
  matrix <- readMM(matrix_path)
  
  # Load barcodes and genes
  barcodes <- readLines(gsub("matrix.mtx", "barcodes.tsv", matrix_path))
  features <- readLines(gsub("matrix.mtx", "features.tsv", matrix_path))
  
  # Set row and column names for the matrix
  rownames(matrix) <- features
  colnames(matrix) <- barcodes
  
  # CreateSeurat object
  seurat_obj <- CreateSeuratObject(
    counts = matrix,
    project = sample,
    min.cells = 3,        # Genes must be expressed in at least 3 cells
    min.features = 200    # Cells must express at least 200 genes
  )
  
  return(seurat_obj)
})

# Merge all Seurat objects into one
seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)


# Add mitochondrial transcript info
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "\tMT-")



# Apply filters Here I am removing cells with less than 200 or greater than 3000 individual transcripts (should get rid of many doublets and other non-descriptive low count cells)
seurat_object <- subset(
  seurat_object,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15
)



#normalizing data (this is all the default currently): for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)



# the gene names were output in a weird tab-delimted format, fixing it for plotting
gene_symbols <- sapply(strsplit(rownames(seurat_object), "\t"), function(x) x[2])
rownames(seurat_object) <- make.unique(gene_symbols)



# let's start by just nabbing highly differentially expressed genes
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)



# Identify the 200 most highly variable genes
top20 <- head(VariableFeatures(seurat_object), 20)






# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


## ok now onto PCA stuff
# Scale the data
seurat_object <- ScaleData(seurat_object)

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))


#checking out loadings (the irritating way gene names are input makes this a lil hard to look at rn 1/15/25)
VizDimLoadings(seurat_object, dims = 1:4, reduction = "pca")

#basic PCA plot w/ 2000 most diff expressed genes
DimPlot(seurat_object, reduction = "pca") + NoLegend()

#now looking at cells from across the spectrum to really check out how things are divided

DimHeatmap(seurat_object, dims = 1:5, cells = 500,  balanced = TRUE)


ElbowPlot(seurat_object)


#now finding nearest neighbors based on PCA space, based on elbow plot I will use PCA 1:14 

seurat_object <- FindNeighbors(seurat_object, dims = 1:14)
seurat_object <- FindClusters(seurat_object, resolution = .5)

seurat_object <- RunUMAP(seurat_object, dims = 1:14)

DimPlot(seurat_object, reduction = "umap", label='true')

seurat_object <- JoinLayers(seurat_object, assay = "RNA")

#let's do a quick screen to confirm what our clusters are
FeaturePlot(seurat_object, features = c("GAD2","SLC17A6","GRIK1", "TH", "CADPS2")) 

FeaturePlot(seurat_object, features = c("VCAN","MOBP","AQP4","FOXJ1", "PDGFRB","CLDN5","CD74","TTR"))


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
names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



cluster0.markers <- FindMarkers(seurat_object, ident.1 = 2, ident.2 = c(13))
head(cluster0.markers, n = 50)
# what are the differences


#let's look at the proportion of GPNMB expressing cells


# Define a detection threshold 
gpnmb_expression <- FetchData(seurat_object, vars = "GPNMB")

# Convert to binary
gpnmb_detected <- gpnmb_expression > 0


# Add GPNMB detection as metadata to Seurat object
seurat_object$GPNMB_Expressed <- gpnmb_detected


seurat_object$cluster_names <- Idents(seurat_object)



# Fetch GPNMB expression and convert to binary (expressed or not)
gpnmb_expression <- FetchData(seurat_object, vars = "GPNMB")
gpnmb_detected <- gpnmb_expression > 2

# Add the binary GPNMB expression to the metadata
seurat_object$GPNMB_Expressed <- gpnmb_detected

# Calculate proportions based on renamed cluster identities
gpnmb_proportion <- seurat_object@meta.data %>%
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



#ok, now let's just subcluster 2

microGliaOnly <- subset(seurat_object, idents = "Microglia")


FeaturePlot(microGliaOnly, features = c("GPNMB","TREM2","APOE",  "CD74","MIF","P2RY12","CX3CR1","HSP90AA1"))


FeaturePlot(seurat_object, features = c("GPNMB", "CD163", "SAT1", "TALAM1", "TREM2","APOE",  "CD74","P2RY12","CX3CR1"))

FeaturePlot(microGliaOnly, features = c("GPNMB", "CD163", "CD74", "ATP1A1", "ATP1A3","KCNMA1"))


#VlnPlot(seurat_object, features = c("GPNMB","P2RY12","TMEM119","IBA1", "TREM2","LPL","CLEC7A","APOE"))




#re-run clustering within the microglia subset 
microGliaOnly <- FindNeighbors(microGliaOnly, dims = 1:14)
microGliaOnly <- FindClusters(microGliaOnly, resolution = 0.25)

# Visualize the clusters within microglia
DimPlot(microGliaOnly, reduction = "umap", label = TRUE)



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

new.cluster.ids <- c("Homeostatic", "Transition","DAM","Homeostatic2")
names(new.cluster.ids) <- levels(microGliaOnly)
microGliaOnly <- RenameIdents(microGliaOnly, new.cluster.ids)


microGliaOnly.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 4) %>%
  ungroup() -> top10
DoHeatmap(microGliaOnly, features = top10$gene) + NoLegend()

## on to frank pseudotime, woohoo!

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")



install.packages("monocle3")
install.packages("SeuratWrappers")

library(monocle3)
library(SeuratWrappers)

# Convert Seurat object to Monocle3 CDS (Cell Data Set)
cds <- as.cell_data_set(microglia_cells)

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)

# Reduce dimensions
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Visualize cells
plot_cells(cds, color_cells_by = "ident")






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




