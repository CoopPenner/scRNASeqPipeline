
#overall set of call functions to implement seurat workflow

#1 get all relevant pacakges 


# Install if not already installed
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
  install.packages("leiden")
  
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

# Visualize QC metrics
print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", plot.cor= "false")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor= "false")
plot1 + plot2


# Apply filters Here I am removing cells with less than 200 or greater than 3000 individual transcripts (should get rid of many doublets and other non-descriptive low count cells)
seurat_object <- subset(
  seurat_object,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15
)




#normalizing data (this is all the default currently): for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# just for a lil sanity check let's re-image qc after these changes

# Visualize QC metrics
print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", plot.cor= "false")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor= "false")
plot1 + plot2





#ok, nice, looks good, 











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




