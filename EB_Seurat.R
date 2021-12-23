library(Seurat)
library(tidyverse)
library(patchwork)

#Read and create a Seurat object
eb.data <- Read10X(data.dir = "EB-scRNAseq/T2_3B")
eb <- CreateSeuratObject(counts = eb.data, project = "embryoid_body",
                          min.cells = 3, min.features = 200)
eb

#Pre-processing
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
eb[["percent.mt"]] <- PercentageFeatureSet(eb, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(eb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <-FeatureScatter(eb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <-FeatureScatter(eb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

eb <-subset(eb, subset = nFeature_RNA >200 & nFeature_RNA <3500 & percent.mt<5)

#Normalizing the data
#normalizes the feature expression measurements for each cell by the total expression, 
#multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
#Normalized values are stored in pbmc[["RNA"]]@data
eb <-NormalizeData(eb, normalization.method = "LogNormalize", scale.factor = 10000)
#same as
#pbmc <- NormalizeData(pbmc)

#Identifying highly variable features. Return 2,000 features per dataset
eb <- FindVariableFeatures(eb, selection.method = "vst", nfeatures = 2000)
#top 10 highly variable genes
top10 <- eb %>% VariableFeatures() %>% head(.,10)
#same as
#top10 <- head(VariableFeatures(pbmc),10)
#plot variable feautres
plot1 <- VariableFeaturePlot(eb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data, aka linear transformation
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(eb)
eb <- ScaleData(eb, features = all.genes)

#Linear dimension reduction and visualization
eb<- RunPCA(eb, features = VariableFeatures(object = eb))
VizDimLoadings(eb, dims = 1:2, reduction = 'pca')
DimPlot(eb, reduction = 'pca') #dimplot
DimHeatmap(eb, dims = 1, cells = 500, balanced = TRUE) #heatmap


#Determine the dimensionality. Processing time can be long
# eb <- JackStraw(eb, num.replicate = 100)
# eb <- ScoreJackStraw(eb, dims = 1:20)

#Significant’ PCs will show a strong enrichment of features with low p-values 
#(solid curve above the dashed line)
# JackStrawPlot(eb, dims = 1:15) #error max dimension is 0

#Alternative - Elbow plot
#Ranking of PC based on the percentage of variance explained by each one 
ElbowPlot(eb)


#Cluster the cells. Resolution sets granularity of the downstream clustering. Bigger more clusters.
#0.4-1.2 typically returns good results for sc datasets around 3k cells
eb <- FindNeighbors(eb, dims = 1:20)
eb <- FindClusters(eb, resolution = 0.5)
eb %>% Idents(.) %>% head(.,5) #show first 5 cell IDs

#Run non-linear dimensional reduction
#UMAP
# eb_umap<-eb
eb_umap <- RunUMAP(eb, dims = 1:20)
DimPlot(eb_umap, reduction = 'umap', label = TRUE)
eb <- eb_umap #save umap result to eb object
#T-SNE
# eb_tsne <-eb
eb_tsne <- RunTSNE(eb, dims = 1:20)
DimPlot(eb_tsne, reduction = 'tsne', label = TRUE)

eb<-eb_tsne #save tsne result while preserving umap result


#Save processed file
# saveRDS(eb, file = 'eb01.rds')

