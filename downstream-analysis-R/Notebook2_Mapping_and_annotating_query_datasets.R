
library(Seurat)

heart <- readRDS("~/heart-int-FND-filtered.Rds")
table(heart$orig.ident)

Idents(heart) <- heart$orig.ident
FY <- subset(heart,subset = orig.ident != "DCM")
table(FY$orig.ident)

FY$orig.ident <- factor(FY$orig.ident, levels = c("Fetal","ND"))

DefaultAssay(heart) <- "RNA"
FY.list <- SplitObject(heart, split.by = "biorep")

for (i in 1:length(FY.list)) {
  FY.list[[i]] <- NormalizeData(FY.list[[i]], verbose = FALSE)
  FY.list[[i]] <- FindVariableFeatures(FY.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}

FY.anchors <- FindIntegrationAnchors(object.list = FY.list, dims = 1:30)
FY.integrated <- IntegrateData(anchorset = FY.anchors, dims = 1:30)
saveRDS(FY.integrated, file="~/FY.integrated.Rds")
FY.integrated <- readRDS(file="~/FY.integrated.Rds")

library(ggplot2)
library(cowplot)
library(patchwork)

DefaultAssay(FY.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
FY.integrated <- ScaleData(FY.integrated, verbose = FALSE)
FY.integrated <- RunPCA(FY.integrated, npcs = 30, verbose = FALSE)
FY.integrated <- RunUMAP(FY.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(FY.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(FY.integrated, reduction = "umap", group.by = "Broad_celltype", label = TRUE, repel = TRUE)
p1 + p2

#-----20240316-smallPatch
heart <- readRDS("~/heart-int-FND-filtered.Rds")
F2 <- subset(heart,subset = orig.ident == "")
DefaultAssay(heart) <- "RNA"
heart <- NormalizeData(heart)
heart <- FindVariableFeatures(heart)
heart <- ScaleData(heart)


DefaultAssay(data) <- "RNA"
query <- NormalizeData(data)
query <- FindVariableFeatures(query)
query <- ScaleData(query)


FY.anchors <- FindTransferAnchors(reference = heart, query = query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = FY.anchors, refdata = heart$Broad_celltype,
                            dims = 1:30)
query <- AddMetaData(query, metadata = predictions)

#query$prediction.match <- query$predicted.id == query$Broad_celltype

table(query$predicted.id)
VlnPlot(query, c("MYH7", "ANLN", "PECAM1", "NRXN1", "PDGFRB", "MYH11"), group.by = "predicted.id")

ctrl <- query
F2 <- RunUMAP(F2, dims = 1:30, reduction = "pca", return.model = TRUE)
ctrl <- MapQuery(anchorset = FY.anchors, reference = F2, query = ctrl,
                           refdata = list(celltype = "Broad_celltype"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(F2, reduction = "umap", group.by = "Broad_celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(ctrl, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)+ ggtitle("Query transferred labels")


p1 + p2

new_media$cellType <- ctrl$predicted.celltype




#----


#let's get prediction.score.Fib-[sum(prediction.score.pericyte+prediction.score.neuron+etc.)]>0.3
table(predictions$predicted.id)

#CM      CM(Prlf)     Endo      Fib   Neuron   Pericyte      Smc 
#2001      343         69       113        1       20        7 

cellbarcodes <- data.frame(matrix(NA, ncol = 1, nrow = dim(ctrl)[2]))
colnames(cellbarcodes) <- "trueCell"
rownames(cellbarcodes) <- rownames(predictions)

index <- which(predictions$predicted.id == "Smc")
for (i in 1:length(index)){
print(predictions$prediction.score.Smc[index[i]]-sum(predictions[index[i],c(2:8,10)]))
print(index[i])
cellbarcodes$trueCell[index[i]] <- (predictions$prediction.score.Smc[index[i]]-sum(predictions[index[i],c(2:8,10)]))> 0.3
}

index <- match(rownames(cellbarcodes)[cellbarcodes$trueCell=="TRUE"],colnames(ctrl))
index <- na.omit(index)
ctrl <- ctrl[,index]

saveRDS(ctrl,file= "~/ctrl-filtered.Rds")
ctrl <- readRDS("~/Endo-Data-Cleaning/ctrl.Rds")
ctrl$Broad_celltype <- ctrl$predicted.celltype
Idents(ctrl) <- ctrl$Broad_celltype

markers <- FindAllMarkers(ctrl, assay = "RNA", slot = "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers %>% filter(markers$p_val_adj < 0.05)
DimPlot(ctrl, reduction = "ref.umap", label = TRUE,group.by = "predicted.celltype")
        
  #"CD93","CDH5","ERG","NOTCH4","PTPRB","TIE1","TSPAN11")      
FeaturePlot(ctrl, features = c("PTPRB"), reduction = "ref.umap")
VlnPlot(ctrl,assay = "RNA", slot = "data", features = c("FLI1"))

barplot(table(Idents(ctrl)),ylab="Number of nuclei",xlab="Clusters", ylim = c(0,2050))
title("Number of nuclei in each cluster")

#CM      CM(Prlf)     Endo      Fib 
#1946      320       16       73