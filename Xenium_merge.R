library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
path1 <- "~/R/Xenium/catalyst_release_SQ_UA_Apr23/0027365_1_f"
path2 <- "~/R/Xenium/catalyst_release_SQ_UA_Apr23/0027365_2_e"
path3 <- "~/R/Xenium/catalyst_release_SQ_UA_Apr23/0027365_3_d"
path4 <- "~/R/Xenium/catalyst_release_SQ_UA_Apr23/0027365_4_c"

raw_all = list()
samples = c('x1','x2','x3','x4')

raw_all[['x1']] = LoadXenium(path1, fov = "fov.1")
raw_all[['x2']] = LoadXenium(path2, fov = "fov.2")
raw_all[['x3']] = LoadXenium(path3, fov = "fov.3")
raw_all[['x4']] = LoadXenium(path4, fov = "fov.4")

for (i in samples) {
  raw_all[[i]] <- subset(raw_all[[i]], subset = nCount_Xenium > 0)
  Idents(raw_all[[i]]) <- i
  raw_all[[i]]$sample <- i
  raw_all[[i]] <- RenameCells(raw_all[[i]], add.cell.id = i)}

merge1 <- merge(raw_all[['x1']], y = c(raw_all[['x2']]))
merge2 <- merge(raw_all[['x3']], y = c(raw_all[['x4']]))
x.merge <- merge(merge1, merge2)


x.merge <- SCTransform(x.merge, assay = "Xenium")
x.merge  <- RunPCA(x.merge , npcs = 30, features = rownames(x.merge))
x.merge  <- RunUMAP(x.merge , dims = 1:30)
x.merge <- FindNeighbors(x.merge, reduction = "pca", dims = 1:30)
x.merge  <- FindClusters(x.merge, resolution = 0.7)
DimPlot(x.merge, reduction = "umap", group.by = c("ident", "sample"))
DimPlot(x.merge, reduction = "umap", group.by = c("seurat_clusters"))

ImageDimPlot(x.merge,fov = c('fov.1', 'fov.2', 'fov.3', 'fov.4'), molecules = c("Vwf", "Ager"))

ImageDimPlot(raw_all[['x1']],fov = c('fov.1'), cols = "polychrome", size = 0.75)
ImageDimPlot(x.merge,fov = c('fov.1', 'fov.2', 'fov.3', 'fov.4'), cols = "polychrome", size = 0.75)
ImageDimPlot(x.merge,fov = c('fov.1'), cols = "polychrome", size = 0.75)

Idents(x.merge) <- "sample"
x1.subset <- subset(x.merge, idents = ("x1"))
x2.subset <- subset(x.merge, idents = ("x2"))
x3.subset <- subset(x.merge, idents = ("x3"))
x4.subset <- subset(x.merge, idents = ("x4"))

library(spacexr)

#Reference
library(SeuratDisk)
Convert("/Users/zhiyudai/pythonProjects/LungMAP_MouseLung_CellRef.v1.1.h5ad", dest = "h5seurat", overwrite = TRUE)
ref <- LoadH5Seurat("/Users/zhiyudai/pythonProjects/LungMAP_MouseLung_CellRef.v1.1.h5seurat")
Idents(ref) <- "celltype_level1"
ref <- subset(ref, idents= c("Chondrocyte"), invert = TRUE)
counts <- GetAssayData(ref, assay = "RNA", layer = "counts")
cluster <- as.factor(ref$celltype_level1)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# For x1 subset
query.counts1 <- GetAssayData(x1.subset, assay = "Xenium", slot = "counts")[, Cells(x1.subset[["fov.1"]])]
coords1 <- GetTissueCoordinates(x1.subset[["fov.1"]], which = "centroids")

# Ensure row names match
rownames(coords1) <- coords1$cell
coords1$cell <- NULL

# Create SpatialRNA object
query1 <- SpatialRNA(coords1, query.counts1, colSums(query.counts1))

# Create and run RCTD
RCTD1 <- create.RCTD(query1, reference, max_cores = 8)
RCTD1 <- run.RCTD(RCTD1, doublet_mode = "doublet")

# Extract and assign annotations
annotations.df1 <- RCTD1@results$results_df
annotations1 <- annotations.df1$first_type
names(annotations1) <- rownames(annotations.df1)
x1.subset$predicted.celltype <- annotations1

# For x2 subset
query.counts2 <- GetAssayData(x2.subset, assay = "Xenium", slot = "counts")[, Cells(x2.subset[["fov.2"]])]
coords2 <- GetTissueCoordinates(x2.subset[["fov.2"]], which = "centroids")

# Ensure row names match
rownames(coords2) <- coords2$cell
coords2$cell <- NULL

# Create SpatialRNA object
query2 <- SpatialRNA(coords2, query.counts2, colSums(query.counts2))

# Create and run RCTD
RCTD2 <- create.RCTD(query2, reference, max_cores = 8)
RCTD2 <- run.RCTD(RCTD2, doublet_mode = "doublet")

# Extract and assign annotations
annotations.df2 <- RCTD2@results$results_df
annotations2 <- annotations.df2$first_type
names(annotations2) <- rownames(annotations.df2)
x2.subset$predicted.celltype <- annotations2

# For x3 subset
query.counts3 <- GetAssayData(x3.subset, assay = "Xenium", slot = "counts")[, Cells(x3.subset[["fov.3"]])]
coords3 <- GetTissueCoordinates(x3.subset[["fov.3"]], which = "centroids")

# Ensure row names match
rownames(coords3) <- coords3$cell
coords3$cell <- NULL

# Create SpatialRNA object
query3 <- SpatialRNA(coords3, query.counts3, colSums(query.counts3))

# Create and run RCTD
RCTD3 <- create.RCTD(query3, reference, max_cores = 8)
RCTD3 <- run.RCTD(RCTD3, doublet_mode = "doublet")

# Extract and assign annotations
annotations.df3 <- RCTD3@results$results_df
annotations3 <- annotations.df3$first_type
names(annotations3) <- rownames(annotations.df3)
x3.subset$predicted.celltype <- annotations3

# For x4 subset
query.counts4 <- GetAssayData(x4.subset, assay = "Xenium", slot = "counts")[, Cells(x4.subset[["fov.4"]])]
coords4 <- GetTissueCoordinates(x4.subset[["fov.4"]], which = "centroids")

# Ensure row names match
rownames(coords4) <- coords4$cell
coords4$cell <- NULL

# Create SpatialRNA object
query4 <- SpatialRNA(coords4, query.counts4, colSums(query.counts4))

# Create and run RCTD
RCTD4 <- create.RCTD(query4, reference, max_cores = 8)
RCTD4 <- run.RCTD(RCTD4, doublet_mode = "doublet")

# Extract and assign annotations
annotations.df4 <- RCTD4@results$results_df
annotations4 <- annotations.df4$first_type
names(annotations4) <- rownames(annotations.df4)
x4.subset$predicted.celltype <- annotations4


sc.cols = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A", "#FF7A5C", "#53377A", "#FF8E00", "#B32851", 
            "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16","#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6","#2E00FF", "#9408F7", "#371377", 
            "#7700FF", "#9E0142", "#FF0080", "#F8766D", "#00A9FF", "#E68613", "#8494FF","#232C16","#FF7A5C", "#F6768E", "#A6BDD7")

p1 <- DimPlot(x1.subset, reduction = "umap", group.by = "predicted.celltype", label = TRUE, cols = sc.cols)
p2 <- DimPlot(x1.subset, reduction = "umap", group.by = "seurat_clusters", cols = sc.cols)
p1 | p2
keep.cells <- Cells(x1.subset)[!is.na(x2.subset$predicted.celltype)]
Idents(x1.subset) <- "predicted.celltype"
x1_predict.markers <- FindAllMarkers(x1.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(x1_predict.markers, file = "~/R/Xenium/x1_predict.markers.csv")

x1.subset_v2 <- subset(x1.subset, cells = keep.cells)
DimPlot(x1.subset_v2, reduction = "umap", group.by = "predicted.celltype")

table(x1.subset_v2$predicted.celltype)

x.merge_an <- merge(x1.subset, x2.subset)
x.merge_an2 <- merge(x3.subset, x4.subset)
x.merge_anf <- merge(x.merge_an, x.merge_an2)
keep.cells <- Cells(x.merge_anf)[!is.na(x.merge_anf$predicted.celltype)]

x.merge_anf2 <- subset(x.merge_anf, cells = keep.cells)

x.merge_anf2  <- RunPCA(x.merge_anf2, npcs = 30, features = rownames(x.merge_anf2))
x.merge_anf2  <- RunUMAP(x.merge_anf2, dims = 1:30)
x.merge_anf2 <- FindNeighbors(x.merge_anf2, reduction = "pca", dims = 1:30)
x.merge_anf2  <- FindClusters(x.merge_anf2, resolution = 0.5)

p1 <- DimPlot(x.merge_anf2, reduction = "umap", group.by = "predicted.celltype", split.by = "sample", label = TRUE, cols = sc.cols)
p2 <- DimPlot(x.merge_anf2, reduction = "umap", group.by = "seurat_clusters", cols = sc.cols)
p1 | p2
DimPlot(x.merge_anf2, reduction = "umap", group.by = "predicted.celltype", split.by = "sample", label = TRUE, cols = sc.cols)
DimPlot(x.merge_anf2, reduction = "umap", group.by = "predicted.celltype", split.by = "sample", cols = sc.cols)

ImageDimPlot(x.merge_anf2,fov = c('fov.1', 'fov.2', 'fov.3', 'fov.4'), group.by = "predicted.celltype", molecules = c("Vwf", "Ager"))
p1 <- ImageDimPlot(x.merge_anf2,fov = c('fov.1'), group.by = "predicted.celltype", cols = "polychrome")
ImageDimPlot(x.merge_anf2,fov = c('fov.1','fov.3'), group.by = "predicted.celltype", cols = "polychrome",coord.fixed = FALSE)

Idents(x.merge_anf2) <- "predicted.celltype"
cell.prop<-as.data.frame(prop.table(table(Idents(x.merge_anf2), x.merge_anf2$sample)))
colnames(cell.prop)<-c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster)) + geom_bar(stat="identity",position="fill") + ggtitle("") + theme_bw() + theme(axis.ticks.length=unit(0.5,'cm')) + guides(fill=guide_legend(title=NULL)) + scale_fill_manual(values = sc.cols) 

x.merge_anf2@meta.data$group <- gsub(pattern = "x1", replacement = "control", x = x.merge_anf2@meta.data$sample)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$group <- gsub(pattern = "x2", replacement = "control", x = x.merge_anf2@meta.data$group)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$group <- gsub(pattern = "x3", replacement = "KO", x = x.merge_anf2@meta.data$group)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$group <- gsub(pattern = "x4", replacement = "KO", x = x.merge_anf2@meta.data$group)# changeST033r to ST033 to newid
DimPlot(x.merge_anf2, reduction = "umap", group.by = "predicted.celltype", split.by = "group", label = TRUE, cols = sc.cols)

cell.prop<-as.data.frame(prop.table(table(Idents(x.merge_anf2), x.merge_anf2$group)))
colnames(cell.prop)<-c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster)) + geom_bar(stat="identity",position="fill") + ggtitle("") + theme_bw() + theme(axis.ticks.length=unit(0.5,'cm')) + guides(fill=guide_legend(title=NULL)) + scale_fill_manual(values = sc.cols) 


x.merge_anf2@meta.data$sample <- gsub(pattern = "x1", replacement = "control1", x = x.merge_anf2@meta.data$sample)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$sample <- gsub(pattern = "x2", replacement = "control2", x = x.merge_anf2@meta.data$sample)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$sample <- gsub(pattern = "x3", replacement = "KO1", x = x.merge_anf2@meta.data$sample)# changeST033r to ST033 to newid
x.merge_anf2@meta.data$sample <- gsub(pattern = "x4", replacement = "KO2", x = x.merge_anf2@meta.data$sample)# changeST033r to ST033 to newid
DimPlot(x.merge_anf2, reduction = "umap", group.by = "predicted.celltype", split.by = "sample",  cols = sc.cols)

library("scProportionTest")
prop_test <- sc_utils(x.merge_anf2)
prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype",
  sample_1 = "control", sample_2 = "KO",
  sample_identity = "group"
)
Cellproportion_group <- permutation_plot(prop_test)
saveRDS(x.merge_anf2, file = "~/R/Xenium/x.merge_anf2.rds")
#DEG
#DEGs for all cell types, D1 control vs Stimulate
Idents(x.merge_anf2) <- "group"
x.merge_anf2$celltype.test <- paste(Idents(x.merge_anf2), x.merge_anf2$predicted.celltype, sep = "_")
Idents(x.merge_anf2) <- "celltype.test"

# Get unique cell types from x.merge_anf2$MergeCellType
cell_types <- unique(x.merge_anf2$predicted.celltype)

# Create an empty list to store results for each cell type
results_list <- list()

slot(object = x.merge_anf2@assays$SCT@SCTModel.list[[3]], name="umi.assay") <-"Xenium" # fix PrepSCTFindMarkers error, [[1]], 1,2,3,4 need to run

x.merge_anf2 <- PrepSCTFindMarkers(x.merge_anf2)
# Loop through each cell type and perform differential expression analysis
for (cell_type in cell_types) {
  celltype_control <- paste("control", cell_type, sep = "_")
  celltype_PH <- paste("PH", cell_type, sep = "_")
  
  # Perform differential expression analysis
  results <- FindMarkers(
    x.merge_anf2, 
    ident.1 =  celltype_PH, 
    ident.2 = celltype_control, 
  )
  
  # Store the results in the list with the cell type as the name
  results_list[[cell_type]] <- results
}

# Now, results_list contains DE analysis for each cell type without mitochondrial genes

library(scRNAtoolVis)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


results_list <- lapply(results_list, function(cell_type) {
  cell_type$gene <- rownames(cell_type)
  return(cell_type)
}) # put the rownames (gene symbol) to a new column gene to avoid the change of gene name after combination into one data frame.

# Ensure cell_types array is defined and matches the names in results_list
cell_types <- names(results_list)
# Add a cell type column to each data frame and combine them into one
degs_combined <- do.call(rbind, lapply(cell_types, function(cell_type) {
  df <- results_list[[cell_type]]
  df$CellType <- cell_type
  return(df)
}))

degs_combined$cluster <- degs_combined$CellType # need cluster column
write.csv(degs_combined, file = "~/R/Xenium/Xenium_DEG.csv")

# Assuming degs_combined is your data frame and it has a column named 'cluster'
cell_types <- unique(degs_combined$cluster)

# Split cell types into two groups (first 6 and next 6)
EC_cell_types <- cell_types[c(1, 4, 13, 15,26)]

# Subset data for the first six cell types
degs_EC <- subset(degs_combined, cluster %in% EC_cell_types)

# Plot using jjVolcano for the first six cell types
mygene <- c("Des", "Mfap4", "Pi16", "Tmem100", "Hpgd")

jjVolcano(diffData = degs_EC, fontface = 'italic', log2FC.cutoff = 0.5, myMarkers = mygene)


# Split cell types into two groups (first 6 and next 6)
Mensenchy_cell_types <- cell_types[c(3, 8, 11, 25,33,32)]
# Subset data for the first six cell types
degs_Mensenchy <- subset(degs_combined, cluster %in% Mensenchy_cell_types)
# Plot using jjVolcano for the first six cell types
mygene <- c("Des", "Mfap4", "Pi16", "Tagln", "Myh11")
jjVolcano(diffData = degs_Mensenchy, fontface = 'italic', log2FC.cutoff = 0.50, myMarkers = mygene)

# Split cell types into two groups (first 6 and next 6)
AT_cell_types <- cell_types[c(2, 5, 10)]
# Subset data for the first six cell types
degs_AT <- subset(degs_combined, cluster %in% AT_cell_types)
# Plot using jjVolcano for the first six cell types
mygene <- c("Krt8", "Krt19")
jjVolcano(diffData = degs_AT, fontface = 'italic', log2FC.cutoff = 0.50, myMarkers = mygene)


#second_9_cell_types <- cell_types[10:18]
# Subset data for the first six cell types
#degs_second_9 <- subset(degs_combined, cluster %in% second_9_cell_types)
# Plot using jjVolcano for the first six cell types
#jjVolcano(diffData = degs_second_9, fontface = 'italic', log2FC.cutoff = 0.585)
#jjVolcano(diffData = degs_second_9, log2FC.cutoff = 0.585, col.type = "adjustP")

#Figure 5
subset_data_AEC <- subset(x.merge_anf2, predicted.celltype == c("AEC", "EPC"))
subset_data_CAP <- subset(x.merge_anf2, predicted.celltype == c("CAP1","CAP2"))
subset_data_EC <- subset(x.merge_anf2, predicted.celltype == c("CAP1","CAP2", "AEC", "EPC", "VEC"))

p1 <- ImageDimPlot(
  object = subset_data_AEC,
  fov = c('fov.1', 'fov.3'),
  split.by = "group",
  cols = c("red", "yellow"), coord.fixed = FALSE
)

p2 <- ImageDimPlot(
  object = subset_data_CAP,
  fov = c('fov.1', 'fov.3'),
  split.by = "group",
  cols = c("orange", "green"), coord.fixed = FALSE
)
EC_markers_featureplot <- p1 | p2 | p3

ggsave("~/R/Xenium/EC_markers_featureplot.pdf", plot =EC_markers_featureplot, width = 20, height = 4, device = "pdf")
ggsave("~/R/Xenium/AEC_featureplot.pdf", plot =p3, width = 6, height = 3, device = "pdf")
ggsave("~/R/Xenium/CAP1_featureplot.pdf", plot =p1, width = 6, height = 3, device = "pdf")
ggsave("~/R/Xenium/CAP2_featureplot.pdf", plot =p2, width = 6, height = 3, device = "pdf")

p3 <- ImageFeaturePlot(subset_data_EC, features = c("Sox17"),fov = c('fov.1', 'fov.3'), size = 0.75, cols = c("white", "red"),coord.fixed = FALSE, min.cutoff =0, dark.background = FALSE,scale = c("all"))
p1 <- ImageFeaturePlot(subset_data_EC, features = c("Plvap"),fov = c('fov.1', 'fov.3'), size = 0.75, cols = c("white", "red"),coord.fixed = FALSE, min.cutoff =1, dark.background = FALSE,scale = c("all"))
p2 <- ImageFeaturePlot(subset_data_EC, features = c("Prx"),fov = c('fov.1', 'fov.3'), size = 0.75, cols = c("white", "red"),coord.fixed = FALSE, min.cutoff =1, dark.background = FALSE,scale = c("all"))


