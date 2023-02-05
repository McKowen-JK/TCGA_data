#use script from https://github.com/cpreid2/gdc-rnaseq-tool to retrieve data from TCGA manifest
libpath = export PATH="/home/keller/R/x86_64-pc-linux-gnu-library/4.1/:$PATH"
library(dplyr)
library(Seurat)
library(patchwork)


# Initialize the Seurat object with the raw (non-normalized data).
counts = read.table(file = "/home/keller/Documents/Mac/cptac.tsv", sep = "\t", header = TRUE, row.names = 1)
matrix = as.matrix(counts)
cptac <- CreateSeuratObject(matrix, project = "TCGA-cptac", min.cells = 20, min.features = 2000)

#normalize
cptac <- NormalizeData(cptac)
cptac <- FindVariableFeatures(cptac, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cptac), 10)

# plot variable features with and without labels
png(file="/home/keller/Documents/Mac/cptac_variable_genes.png",
width=600, height=350)
VariableFeaturePlot(cptac)
dev.off()

#plot1 <- VariableFeaturePlot(cptac)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

#scaling data
all.genes <- rownames(cptac)
cptac <- ScaleData(cptac, features = all.genes)

cptac <- RunPCA(cptac, features = VariableFeatures(object = cptac))

#pca plots
png(file="/home/keller/Documents/Mac/cptac_pca_x.png",
width=600, height=400)
DimPlot(cptac, dims = c(1, 2), reduction = "pca")
dev.off()


#Find significant PCAs
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
cptac <- JackStraw(cptac, num.replicate = 100)
cptac <- ScoreJackStraw(cptac, dims = 1:20)

png(file="/home/keller/Documents/Mac/cptac_PCA_evaluationx.png",
width=1000, height=600)
JackStrawPlot(cptac, dims = 1:20)
dev.off()

#14 significant PCAs!

#determine parameters

cptac <- FindNeighbors(cptac, dims = 1:12)
cptac <- FindClusters(cptac, resolution = 1.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
cptac <- RunUMAP(cptac, dims = 1:12, umap.method = 'umap-learn', metric = 'correlation')

png(file="/home/keller/Documents/Mac/cptac_UMAPx.png",
width=600, height=400)
DimPlot(cptac, reduction = "umap")
dev.off()

saveRDS(pbmc, file = "../output/cptac_mac.rds")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(cptac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

