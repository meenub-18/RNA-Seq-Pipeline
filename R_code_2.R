##--------------------------- Load packages ---------------------------

#install.packages("Seurat")
#BiocManager::install("Matrix")

devtools::install_github('satijalab/seurat-data')
options(timeout = 300)
BiocManager::install("SeuratData")

library(Seurat)
library(SeuratData)
library(ggplot2)

##---------------------------- check for avialable data in Seurat Data package ------------------------
AvailableData()

##---------------------------- Install one demo data i.e stxBrain data -------------------------------
#InstallData("stxBrain")
InstallData("stxKidney")

# Load the data
brain <- LoadData("stxBrain", type = "anterior1")
brain

#Just like the brain dataset, we can load the kidney detaset using below command. However in kidney detaset there is no specification like anterior and posterior ( similar to the brain dataset). It will take kidney detaset as default.
#kidney <- LoadData("stxKidney")

# check the metadata
View(brain@meta.data)

# check the rownames and colnames of brain obj
View(colnames(brain))
View(rownames(brain))

# qc plots
plot_a <- VlnPlot(brain, features = "nFeature_Spatial", pt.size = 1.5)
plot_a


plot_b <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 1.5)

plot_a + plot_b

SpatialFeaturePlot(brain, features = "nFeature_Spatial")



# load the data from 10X dataset
getwd() # get working dir

list.files() # list out all the files for the current working directory

brain_data <- Load10X_Spatial(data.dir = "H:/decode_session/new_data/",
                              filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_raw_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "slice1")
brain_data


View(brain_data@meta.data)

# Dimension of data
dim(x = brain_data)

# number of rows
nrow(x = brain_data)

# see the staring few line contents
head(x = rownames(brain_data), n = 15)

head(x = colnames(brain_data), n = 15)



# filtering
brain_subset <- subset(brain, subset = nFeature_Spatial < 7500)
brain_subset

print(paste("Filter out", ncol(brain) - ncol(brain_subset), 
            "samples because of the outlier QC metrics, with", ncol(brain_subset),
            "samples left."))

names(brain)



# Normalization
brain_norm <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
names(brain_norm)

dim(brain_norm@assays$SCT@scale.data)

brain_obj <- RunPCA(brain_norm, assay = "SCT", verbose = FALSE)
brain_obj <- FindNeighbors(brain_obj, reduction = "pca", dims = 1:30)
brain_obj <- FindClusters(brain_obj, verbose = TRUE)

# Runing umap
brain_obj <- RunUMAP(brain_obj, reduction = "pca", dims = 1:30)

# generate umap plot
plot_umap <- DimPlot(brain_obj, reduction = "umap", label = TRUE) + NoLegend()
plot_umap

plot_spatial <- SpatialDimPlot(brain_obj, label = TRUE, label.size = 2)
plot_spatial




packageVersion("Seurat")
packageVersion("SeuratData")
