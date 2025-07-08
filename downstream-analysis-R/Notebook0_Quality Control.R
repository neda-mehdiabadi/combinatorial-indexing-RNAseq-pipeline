# title: "Quality control of the cells"
# author: "Neda R. Mehdiabadi"
# date: "22/01/2022"


# Introduction
# The first step of the analysis is to perform quality control of the cells to make sure that low quality cells are removed prior to further analysis

# Load libraries and functions
setwd("/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/code/sci-plex-master/large_screen")
library(Seurat)
library(dplyr)
library(Homo.sapiens)
library(edgeR)
source("../../../../GENERAL_CODES/ggplotColors.R")
source("../../../../GENERAL_CODES/normCounts.R")
source("../../../../GENERAL_CODES/findModes.R")

# Read in data (seurat object) generated in "./Notebook3_preprocess_data.R"
alldata <- readRDS("../../../output/combo/all15098-seurat.Rds")
#data <- readRDS("../../../output/second-batch/all33817-seurat.Rds")
alldata$replicate <- factor(alldata$replicate, levels = paste("rep",1:8, sep = ""))
alldata$orig.ident <- factor(alldata$orig.ident, levels = c("WT_TTN","patient_het.50","patient_hom.42","patient_hom.59"))

#extract a batch 
data <- subset(alldata,subset = batch == "batch2" & orig.ident %in% c("patient_120i.3","patient_120i.4","patient_120i.6","WT_TPM1") & media =="MM")
data <- subset(alldata,subset = batch == "batch2" & orig.ident %in% c("patient_DSPKI","WT_DSP") & media =="MM")
data$orig.ident <- factor(data$orig.ident, levels = c("patient_DSPKI","WT_DSP"))
data$orig.ident <- factor(data$orig.ident, levels = c("patient_120i.3", "patient_120i.4", "patient_120i.6","WT_TPM1"))

list <- list()
for (cell_line in levels(factor(data$orig.ident))){
  for (density in levels(factor(data$density[which(data$orig.ident==cell_line)]))){
    for (media in levels(factor(data$media[which(data$density==density)]))){
      for (seed in levels(factor(data$seed[which(data$media==media)]))){
    list[[paste(cell_line,density,media, seed,replicate, sep = "_")]] <- data[,which(data$orig.ident==cell_line & data$density == density & data$media==media & data$seed == seed & data$replicate == replicate)] 
  }
    
  }
}
}

#list <- SplitObject(data, split.by = "orig.ident")
data$split <- paste(data$orig.ident,data$media, data$seed,data$replicate, sep = "_")
list <- SplitObject(data, split.by = "split")
all <- list()
for ( i in 1:length(list)){
  all[[i]] <- list[[i]]@assays$RNA@counts
  names(all)[[i]] <- names(list)[[i]]
}

for (i in 1:length(all)){
  #Take out mitochondrial, ribosomal and genes with no annotation
  ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all[[i]]),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
  m <- match(rownames(all[[i]]),ann$ENSEMBL)
  ann <- ann[m,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #20633
  mito <- grep("mitochondrial",ann$GENENAME)
  length(mito)
  #[1] 229
  ribo <- grep("ribosomal",ann$GENENAME)
  length(ribo)
  #[1] 201
  missingEZID <- which(is.na(ann$ENTREZID))
  length(missingEZID)
  #[1] 1366
  
  # Filter out non-informative genes
  chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
  length(chuck)
  all[[i]] <- all[[i]][-chuck,]
  ann <- ann[-chuck,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
#  num_cells_expressed <- rowSums(all[[i]] !=0)
#  fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
#  keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
#  table(keep.genes)
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  #num_cells_expressed <- rowSums(all[[i]] !=0)
  numzero.genes <- rowSums(all[[i]] ==0)
  #fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
  #keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
  table(numzero.genes > (ncol(all[[i]])-20))
  keep.genes <- numzero.genes < (ncol(all[[i]])-20)
  all[[i]] <- all[[i]][keep.genes,]
  dim(all[[i]])
  ann <- ann[keep.genes,]
  
  #remove low quality nuclei
  keep.nuclei <- colSums(all[[i]] !=0 ) >= 200
  all[[i]]  <- all[[i]][,keep.nuclei]
  
  # Take out sex chromosome genes
  xy <- ann$TXCHROM %in% c("chrX","chrY")
  table(xy)
  xy
  #FALSE  TRUE 
  #13571   516
  all[[i]] <- all[[i]][!xy,]
  ann <- ann[!xy,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #13571
}  


saveRDS(all,"../../../output/combo/batch2-dgCMatrix-filter-TPM1-MM-orig.ident.Rds")
all <- readRDS("../../../output/combo/all-firstbatch-dgCMatrix-filter.Rds")

for (i in 1:length(all)){
  index <- which(dim(all[[i]])[2]>50)
}

orig.ident <- list()
for (i in 1:length(all)){
  orig.ident[[i]] <- rep(levels(factor(list[[i]]@meta.data$orig.ident)),ncol(all[[i]]))
  names(orig.ident)[[i]] <- names(all)[[i]]
}

replicate <- list()
for ( i in 1:length(all)){
  index <- match(colnames(all[[i]]),colnames(list[[i]]))
  replicate[[i]] <- list[[i]]$replicate[index]
  names(replicate)[[i]] <- names(all)[[i]]
}

treatment <- list()
for ( i in 1:length(all)){
  treatment[[i]] <- rep(levels(factor(list[[i]]@meta.data$treatment)),ncol(all[[i]]))
  names(treatment)[[i]] <- names(all)[[i]]
}

density <- list()
for ( i in 1:length(all)){
  density[[i]] <- rep(levels(factor(list[[i]]@meta.data$density)),ncol(all[[i]]))
  names(density)[[i]] <- names(all)[[i]]
}
seed <- list()
for ( i in 1:length(all)){
  index <- match(colnames(all[[i]]),colnames(list[[i]]))
  seed[[i]] <- list[[i]]$seed[index]
  names(seed)[[i]] <- names(all)[[i]]
}
media <- list()
for ( i in 1:length(all)){
  index <- match(colnames(all[[i]]),colnames(list[[i]]))
  media[[i]] <- list[[i]]$media[index]
  names(media)[[i]] <- names(all)[[i]]
}

sciplex <- list()
for (i in 1:length(all)){
  sciplex[[i]] <- CreateSeuratObject(counts = all[[i]], project = names(all)[[i]])
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = orig.ident[[i]], col.name="orig.ident")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = replicate[[i]], col.name="replicate")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = density[[i]], col.name="density")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = seed[[i]], col.name="seed")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = media[[i]], col.name="media")
  #sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = treatment[[i]], col.name="treatment")
  names(sciplex)[[i]] <- names(all)[[i]]
}

for (i in 1:length(sciplex)){
  sciplex[[i]] <- NormalizeData(sciplex[[i]],scale.factor = 10000)
  sciplex[[i]] <- FindVariableFeatures(sciplex[[i]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(sciplex[[i]])
  sciplex[[i]] <- ScaleData(sciplex[[i]], features = all.genes)
  sciplex[[i]] <- RunPCA(sciplex[[i]], features = VariableFeatures(object = sciplex[[i]]))
  sciplex[[i]] <- RunUMAP(sciplex[[i]], reduction = "pca", dims = 1:15)
  sciplex[[i]] <- FindNeighbors(sciplex[[i]], dims = 1:15)
  sciplex[[i]] <- FindClusters(sciplex[[i]], resolution = 0.1)
}

sciplex <- lapply(X = sciplex, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
for (i in 1:length(sciplex)){
  if(dim(sciplex[[i]])[2]<50){
    print(i)
  }
}
sciplex <- sciplex[c(-100,-107,-125,-128)]
min <- min(sapply(sciplex, ncol))
features <- SelectIntegrationFeatures(object.list = sciplex)
heart.anchors <- FindIntegrationAnchors(object.list = sciplex, anchor.features = features, k.filter = min)
heart.combined  <- IntegrateData(anchorset = heart.anchors,dims=1:25)

sciplex[[12]] <- FindClusters(sciplex[[12]], resolution = 0.2)
sciplex[[12]] <- RunUMAP(sciplex[[12]], reduction = "pca", dims = 1:25)
DimPlot(sciplex[[1]], reduction = "umap",label=T,label.size = 5)

sciplex[[1]]$Broad_celltype <- "NA"
index <- which(sciplex[[3]]$RNA_snn_res.0.1=="2")
sciplex[[3]]$Broad_celltype[index] <- "Fib"
saveRDS(sciplex,file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-3rd-29.04.2022/output/sciplex-annotatedCellType.Rds")
Idents(sciplex[[1]]) <- sciplex[[1]]$Broad_celltype


sciplex.markers <- Markers(sciplex[[1]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(sciplex[[1]]),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(sciplex[[1]]),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(sciplex[[1]]))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])
write.csv(sciplex.markers,"../../../output/docs/WTmutation_flash_res0.2_markers.csv")
filter <- sciplex.markers %>% filter(sciplex.markers$cluster ==3)
                                       #MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2           "TNNT2"           #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(sciplex[[1]], features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)
FeaturePlot(sciplex[[1]], features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3, label = T)

VlnPlot(sciplex[[1]], features="ENSG00000198626")
FeaturePlot(sciplex[[1]], features = c("ENSG00000198626"), ncol = 3, label = T)
                                        #MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2          #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(sciplex[[1]], features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)


sciplex.sub <- merge(x = sciplex[[1]],y = sciplex[[14]])
sciplex.sub <- subset(sciplex.sub, subset=Broad_celltype != "NA")

heart.list <- SplitObject(sciplex.sub, split.by = "orig.ident")

#sciplex.list <- list()
#for ( i in 7:9){
#  sciplex.list[[i]] <- sciplex[[i+3]]
#  names(sciplex.list)[[i]] <- names(sciplex)[[i+3]]
#}
sciplex.list <- sciplex
min_filter <- min(sapply(sciplex.list, ncol))
for (i in 1:length(sciplex.list)) {
  sciplex.list[[i]] <- SCTransform(sciplex.list[[i]], verbose = FALSE)
}

for (i in 1:length(sciplex.list)) {
  sciplex.list[[i]] <- NormalizeData(sciplex.list[[i]], verbose = FALSE)
  sciplex.list[[i]] <- FindVariableFeatures(sciplex.list[[i]], selection.method = "vst", nfeatures = 2000,
                                     verbose = FALSE)
}

###extra https://satijalab.org/seurat/articles/integration_introduction.html
#Take out mitochondrial, ribosomal and genes with no annotation
data <- subset(alldata,subset = orig.ident %in% c("WT_TTN","patient_het.50") & media =="RPMI")
data <- subset(alldata,subset = orig.ident %in% c("WT_TTN","patient_hom.42","patient_hom.59") & media =="RPMI")
#data$orig.ident <- factor(data$orig.ident, levels = c("patient_120i.3", "patient_120i.4", "patient_120i.6","WT_TPM1"))
#data <- subset(alldata,subset = batch == "batch2" & orig.ident %in% c("patient_DSPKI","WT_DSP") & media =="MM")
#data$orig.ident <- factor(data$orig.ident, levels = c("patient_DSPKI","WT_DSP"))

ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 229
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 1366

# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
data <- data[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(data))

#remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
#  num_cells_expressed <- rowSums(all[[i]] !=0)
#  fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
#  keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
#  table(keep.genes)

#remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
#num_cells_expressed <- rowSums(all[[i]] !=0)
numzero.genes <- rowSums(data@assays$RNA@counts ==0)
#fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
#keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
table(numzero.genes > (ncol(data@assays$RNA@counts)-20))
keep.genes <- numzero.genes < (ncol(data@assays$RNA@counts)-20)
data <- data[keep.genes,]
dim(data)
ann <- ann[keep.genes,]

#remove low quality nuclei
keep.nuclei <- colSums(data@assays$RNA@counts !=0 ) >= 200
data  <- data[,keep.nuclei]

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
xy
#FALSE  TRUE 
#12891   489 
data <- data[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(data))
#TRUE 
#11185

data$stim <- NA
index <- which(data$orig.ident %in% c("patient_120i.3", "patient_120i.4", "patient_120i.6"))
#index <- which(data$orig.ident %in% c("WT_TPM1"))
data$stim[index] <- "MUT"

##https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
data <- NormalizeData(data,scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
ElbowPlot(data,ndims=50)
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.1)
data <- RunUMAP(data, reduction = "pca", dims = 1:15)
pdf(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/output/combo/docs/UMAP-DSPKI-no-batch.pdf",width =7,height = 3)
DimPlot(data, reduction = "umap", split.by ="orig.ident")
dev.off()
Idents(data) <- data$RNA_snn_res.0.1
pdf(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/output/combo/docs/UMAP-DSPKI-COL1A1-RYR2-MYH7-MKI67-TOP2A.pdf",width = 5,height = 15)
FeaturePlot(data, features = c("ENSG00000092054"), split.by = "orig.ident", reduction = "umap")
dev.off()
data$Broad_celltype <- NA
index <- which(data$RNA_snn_res.0.1=="1")
data$Broad_celltype[index] <- "CM(Prlf)"
Idents(data) <- data$Broad_celltype


sciplex.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])
write.csv(sciplex.markers,"../../../output/docs/WTmutation_flash_res0.2_markers.csv")
filter <- sciplex.markers %>% filter(sciplex.markers$cluster ==3)
#MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2           "TNNT2"           #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(invitro, features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)
VlnPlot(invitro, features = c("ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245"),ncol=3) 
FeaturePlot(heart.integrated, features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426"),split.by = "orig.ident", ncol = 3, label = T)

VlnPlot(sciplex[[1]], features="ENSG00000164692")
#MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2          #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(sciplex[[1]], features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)


data$split <- NA
data$split <- paste(data$orig.ident,data$density,data$seed,sep = "_")
data$split <- factor(data$split)
data.list <- SplitObject(data, split.by = "split")
pseudobulk <- matrix(NA,ncol=length(data.list),nrow=nrow(data.list[[1]]))
colnames(pseudobulk) <- names(data.list)
rownames(pseudobulk) <- rownames(data.list[[1]])
for (i in 1:length(data.list))
  pseudobulk[,i] <- rowSums(data.list[[i]])

y <- DGEList(pseudobulk)

ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data))
table(rownames(y) %in% ann$ENSEMBL)

m <- match(rownames(y),ann$ENSEMBL)
y$genes <- ann[m,]
#saveRDS(y,file="./data/pseudobulk.Rds")
keep <- rowSums(y$counts)>8
#keep <- rowSums(cpm.DGEList(y)>=0.5)>=8
table(keep)
y$samples$group <- c(rep("patient_120i.3_400_100",1),rep("patient_120i.3_400_50",1),rep("patient_120i.4_300_100",1),rep("patient_120i.4_300_50",1),rep("patient_120i.6_250_100",1),rep("patient_120i.6_250_50",1),rep("WT_TPM1_400_100",1),rep("WT_TPM1_400_50",1))
y$samples$group <- factor(y$samples$group, levels = c("patient_120i.3_400_100","patient_120i.3_400_50","patient_120i.4_300_100","patient_120i.4_300_50","patient_120i.6_250_100","patient_120i.6_250_50","WT_TPM1_400_100","WT_TPM1_400_50"))

y.keep <- y[keep,]
y.keep <- calcNormFactors(y.keep)
par(mfrow=c(1,1))
par(mar=c(4,4,2,2))
plotMDS(y.keep, top = 500,dim=c(1,2), cex = 2, gene.selection = "common") 
legend("bottom",fill=ggplotColors(4),legend=levels(factor(targets$Group)),bty="n")

keep <- rowSums(cpm.DGEList(y)> 2)>= 1  #median(y$samples$lib.size) is ~135,694. A CPM of 74 is used as it corresponds to a count of 10 for this median library size (~159,988) in this data set. 8 = minimum number of replicates per condition.
table(keep)
#keep
#FALSE  TRUE 
#13880  4818
y <- y[keep,]
y <- calcNormFactors(y)
dim(y)
#[1] 4818  192

### for first-batch
y$samples$group <- c(rep("MCHTB11.4_700",16),rep("MCHTB11.4_800",16),rep("MCHTB11.4_900",8),rep("MCHTB11.5_700",16),rep("MCHTB11.5_800",16),rep("MCHTB11.5_900",8),rep("MCHTB11.9_600",16),rep("MCHTB11.9_700",16),rep("MCHTB11.9_800",8),rep("MCHTB11.9 H3_700",16),rep("MCHTB11.9 H3_800",16),rep("MCHTB11.9 S2_700",16),rep("MCHTB11.9 S2_800",16),rep("MCHTB11.9 S2_900",8))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_700","MCHTB11.4_800","MCHTB11.4_900","MCHTB11.5_700","MCHTB11.5_800","MCHTB11.5_900","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800","MCHTB11.9 S2_700","MCHTB11.9 S2_800","MCHTB11.9 S2_900"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = c(brewer.pal(name="Paired", n = 12), c("chocolate","violetred4"))
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 5)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),14)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomright", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("MCHTB11.4_700","MCHTB11.4_800","MCHTB11.4_900","MCHTB11.5_700","MCHTB11.5_800","MCHTB11.5_900","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800","MCHTB11.9 S2_700","MCHTB11.9 S2_800","MCHTB11.9 S2_900"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)
####
features <- SelectIntegrationFeatures(object.list = sciplex.list)
heart.anchors <- FindIntegrationAnchors(object.list = sciplex.list, dims=1:30,anchor.features = features, k.filter = min_filter)
heart.integrated <- IntegrateData(anchorset = heart.anchors,dims=1:30)
DefaultAssay(object = heart.integrated) <- "integrated"
heart.integrated <- ScaleData(heart.integrated, verbose = FALSE)
saveRDS(data,file="/group/card3/ACTIVE/Neda/sci-RNA-seq3-7th-TTNtv/output/combo/int-vitro-TTNhom-RPMI.Rds")
heart.integrated <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/output/combo/int-vitro-batch1-TPM1-RPMI.Rds")

set.seed(10)
heart.integrated <- RunPCA(heart.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(heart.integrated,ndims=50)
heart.integrated <- FindNeighbors(heart.integrated, dims = 1:15)
heart.integrated <- FindClusters(heart.integrated, resolution = 0.1)
set.seed(10)
heart.integrated <- RunUMAP(heart.integrated, reduction = "pca", dims = 1:15)

##annotate cell types
heart.integrated$Broad_celltype <- NA
index <- which(heart.integrated$integrated_snn_res.0.1=="2")
heart.integrated$Broad_celltype[index] <- "Fib"
Idents(heart.integrated) <- heart.integrated$Broad_celltype

Idents(heart.integrated) <- heart.integrated$Broad_celltype
pdf(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/output/combo/docs/UMAP.pdf",width = 10,height = 3)
DimPlot(heart.integrated, reduction = "umap",label.size = 5, split.by = "orig.ident")
dev.off()

sciplex.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])

DefaultAssay(heart.integrated) <- "SCT"
FeaturePlot(heart.integrated, features = c("ENSG00000092054"), split.by = "orig.ident", ncol = 3, label = T)
VlnPlot(heart.integrated, features = c("ENSG00000092054"), split.by = "orig.ident", ncol = 3)

y.fetal <- DGEList(heart.integrated@assays$RNA@counts)
y.fetal$genes <- rownames(y.fetal)

logcounts <- normCounts(y.fetal,log=TRUE,prior.count=0.5)

DefaultAssay(heart.integrated) <- "RNA"
maxclust <- length(levels(Idents(heart.integrated)))-1

grp <- paste("c",Idents(heart.integrated),sep = "")
grp <- factor(grp,levels = paste("c",0:maxclust,sep=""))

design <- model.matrix(~0+grp)
colnames(design) <- levels(grp)

mycont <- matrix(NA,ncol=length(levels(grp)),nrow=length(levels(grp)))
rownames(mycont)<-colnames(mycont)<-levels(grp)
diag(mycont)<-1
mycont[upper.tri(mycont)]<- -1/(length(levels(factor(grp)))-1)
mycont[lower.tri(mycont)]<- -1/(length(levels(factor(grp)))-1)

fit <- lmFit(logcounts,design)
fit.cont <- contrasts.fit(fit,contrasts=mycont)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

fit.cont$genes <- rownames(y.fetal)

summary(decideTests(fit.cont))

treat <- treat(fit.cont,lfc=0.25)

dt<-decideTests(treat)

summary(dt)

contnames <- colnames(mycont)


for(i in 1:length(contnames)){
  topsig <- topTreat(treat,coef=i,n=Inf,p.value=0.05)
  ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(topsig),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
  m <- match(rownames(topsig),ann$ENSEMBL)
  ann <- ann[m,]
  table(ann$ENSEMBL==rownames(topsig))
  topsig$SYMBOL <- ann$SYMBOL
  write.csv(topsig[topsig$logFC>0,],file=paste("../../../output/second-batch/Up-ClusterTPM-",contnames[i],".csv",sep=""))
  #write.csv(topGO(goana(de=topsig$ENTREZID[topsig$logFC>0],universe=treat$genes$ENTREZID,species="Hs"),number=50),
  #          file=paste("../../../output/first-batch/GO-Cluster-",contnames[i],".csv",sep=""))
  
}

FeaturePlot(heart.integrated, features = c("ENSG00000092054"), split.by = "orig.ident", ncol = 3, label = T)
VlnPlot(heart.integrated, features = c("ENSG00000092054"), split.by = "orig.ident", ncol = 3)


heart.integrated$Broad_celltype <- NA

index <- which(heart.integrated$integrated_snn_res.0.1 =="4")
index <- na.omit(index)
heart.integrated$Broad_celltype[index] <- "Fib"
Idents(heart.integrated) <- heart.integrated$Broad_celltype
DefaultAssay(heart.integrated) <- "integrated"
DimPlot(heart.integrated, reduction = "umap",label=TRUE,label.size = 6,pt.size = 0.5,split.by = "orig.ident")+NoLegend()
DimPlot(invitro, reduction = "umap",label=TRUE,label.size = 6,pt.size = 0.5,split.by = "density")+NoLegend()
saveRDS(heart.integrated,file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-5th-DSP-TPM/output/second-batch/TPM-int-vitro.Rds")
invitro <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-5th-DSP-TPM/output/first-batch/int-vitro.Rds")
index <- which(invitro$orig.ident=="WT_S2")
index <- na.omit(index)
invitro$orig.ident[index] <- "MCHTB11.9 S2"

sciplex.markers <- FindAllMarkers(data.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data.int),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data.int),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data.int))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])
write.csv(sciplex.markers,"../../../output/docs/WTmutation_flash_res0.2_markers.csv")
filter <- sciplex.markers %>% filter(sciplex.markers$cluster ==3)
                                   #MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2           "TNNT2"           #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(invitro, features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)
VlnPlot(invitro, features = c("ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000118194","ENSG00000111245"),ncol=3) 
FeaturePlot(heart.integrated, features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426"),split.by = "orig.ident", ncol = 3, label = T)

VlnPlot(sciplex[[1]], features="ENSG00000164692")
                                        #MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN            #MYBPC3         #MYH6              #RYR2          #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
VlnPlot(sciplex[[1]], features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)


###Make MDS plot
data <- readRDS("../../../../sci-RNA-seq3-5th-DSP-TPM/output/second-batch/all33817-seurat.Rds")
data$replicate <- factor(data$replicate, levels = paste("rep",1:8, sep = ""))

data.int <- readRDS("../../../../sci-RNA-seq3-5th-DSP-TPM/output/second-batch/int-vitro.Rds")
#index <- match(colnames(data.int),colnames(data))
#index <- na.omit(index)
#data.int$density <- paste(data$orig.ident[index],data$density[index],sep = "_")

DefaultAssay(data) <- "RNA"
DimPlot(data, reduction = "umap",label=FALSE,label.size = 6,pt.size = 0.2,split.by = "density")
DefaultAssay(data) <- "SCT"    
FeaturePlot(heart.integrated, features = c("ENSG00000092054"),split.by = "orig.ident",
            pt.size = 0.1,ncol = 3, label = F)
data$Broad_celltype <- factor(data.int$Broad_celltype, levels = c("CM","CM(Prlf)","Fib","unknown"))
VlnPlot(data.int, features=c("ENSG00000175206","ENSG00000148677"))
data <- data.frame(table(data.int$Broad_celltype, data.int$orig.ident))
colnames(data) <- c("CellType", "Group","Percentage")
for ( i in 1:3){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[1]*100
}

for ( i in 4:6){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[2]*100
}
for ( i in 7:9){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[3]*100
}
for ( i in 13:16){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[4]*100
}

for ( i in 17:20){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[5]*100
}
for ( i in 16:18){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[6]*100
}
for ( i in 19:21){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$orig.ident)[7]*100
}

for ( i in 22:24){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[8]*100
}
for ( i in 25:27){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[9]*100
}
for ( i in 28:30){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[10]*100
}

for ( i in 31:33){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[11]*100
}
for ( i in 34:36){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[12]*100
}
for ( i in 37:39){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[13]*100
}

for ( i in 40:42){
  data$RelativePercentage[i] <- data$Percentage[i]/table(data.int$density)[14]*100
}

options(digits=2)
library(ggplot2)
ggplot(data, aes(fill=CellType, y=RelativePercentage, x=Group)) +
  geom_bar(position="stack", stat="identity")+scale_fill_brewer(palette = "Paired")+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+labs(y="Relative Percentage [%]", x = "Group")


DefaultAssay(data.int) <- "RNA"
data.int <- subset(data.int,subset = orig.ident !="MCHTB11.9 S2")
cm <- subset(data.int,subset = Broad_celltype =="Fib")
dim(cm)
#9994 42951
index <-   index <- match(colnames(cm),colnames(data))
index <- na.omit(index)
data.second <- data[,index]
cm$split <- NA
cm$split <- paste(cm$density,"second",sep = "_")
cm.list <- SplitObject(cm,split.by = "split")
backup <- cm.list
list <- c(backup,cm.list)
cm.list <- list
data <- cbind(data.first@assays$RNA@counts,data.second@assays$RNA@counts)
numzero.genes <- rowSums(data==0)
table(numzero.genes > (ncol(data)-20))
keep.genes <- numzero.genes < (ncol(data)-20)
table(keep.genes)
data <- data[keep.genes,]
dim(data)

ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(data),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(data),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(data))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 231
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 775

# Filter out non-informative genes
#ENSG00000251562 stands for MALAT1
chuck <- unique(c(mito,ribo,missingEZID,which(rownames(data)=="ENSG00000251562")))
length(chuck)
#[1] 1127

data <- data[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(data))
#TRUE 
#19506

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
#xy
#FALSE  TRUE 
#18698   808
data <- data[!xy,]
ann <- ann[!xy,]


condition <- as.factor(c(rep("mut",120),rep("wt",32),rep("mut",32),rep("wt",32)))
batch = as.factor(c(rep("1",72),rep("2",64)))
design0 <- model.matrix(~0+condition+batch)
dim(design0)
rownames(design0) <- colnames(y)

y <- DGEList(data)
keep <- rowSums(cpm.DGEList(y)> 1)>=8  #median(y$samples$lib.size) is ~135,694. A CPM of 74 is used as it corresponds to a count of 10 for this median library size (~159,988) in this data set. 8 = minimum number of replicates per condition.
table(keep)
#keep
#FALSE  TRUE 
#13880  4818
y <- y[keep,]
y <- calcNormFactors(y)
dim(y)

for (i in 1:length(cm.list)){
  index <- match(colnames(cm.list[[i]]),colnames(data))
  index <- na.omit(index)
  cm.list[[i]] <- data[,index]
}

# MDS plot of all samples
#To get a high-level idea of the overall sources of variability in the dataset, I have summed the counts over all cells within a sample to obtain a “pseudobulk” sample and made MDS plots using functions in edgeR.
for (i in 1:length(cm.list)){
  #remove low quality nuclei
  keep.nuclei <- colSums(cm.list[[i]] !=0 ) >= 200
  cm.list[[i]]  <- cm.list[[i]][,keep.nuclei]
}  

pseudobulk <- matrix(NA,ncol=length(cm.list),nrow=nrow(cm.list[[1]]))
colnames(pseudobulk) <- names(cm.list)
rownames(pseudobulk) <- rownames(cm.list[[1]])
for (i in 1:length(cm.list))
  pseudobulk[,i] <- rowSums(cm.list[[i]])

#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(pseudobulk),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(pseudobulk),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(pseudobulk))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 231
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 775

# Filter out non-informative genes
#ENSG00000251562 stands for MALAT1
chuck <- unique(c(mito,ribo,missingEZID,which(rownames(pseudobulk)=="ENSG00000251562")))
length(chuck)
#[1] 1127

pseudobulk <- pseudobulk[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(pseudobulk))
#TRUE 
#19506

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
#xy
#FALSE  TRUE 
#18698   808
pseudobulk <- pseudobulk[!xy,]
ann <- ann[!xy,]

#Filter lowly expressed genes
y <- DGEList(pseudobulk)
keep <- rowSums(cpm.DGEList(y)> 2)>= 1  #median(y$samples$lib.size) is ~135,694. A CPM of 74 is used as it corresponds to a count of 10 for this median library size (~159,988) in this data set. 8 = minimum number of replicates per condition.
table(keep)
#keep
#FALSE  TRUE 
#13880  4818
y <- y[keep,]
y <- calcNormFactors(y)
dim(y)
#[1] 4818  192

### for first-batch
y$samples$group <- c(rep("MCHTB11.4_700",16),rep("MCHTB11.4_800",16),rep("MCHTB11.4_900",8),rep("MCHTB11.5_700",16),rep("MCHTB11.5_800",16),rep("MCHTB11.5_900",8),rep("MCHTB11.9_600",16),rep("MCHTB11.9_700",16),rep("MCHTB11.9_800",8),rep("MCHTB11.9 H3_700",16),rep("MCHTB11.9 H3_800",16),rep("MCHTB11.9 S2_700",16),rep("MCHTB11.9 S2_800",16),rep("MCHTB11.9 S2_900",8))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_700","MCHTB11.4_800","MCHTB11.4_900","MCHTB11.5_700","MCHTB11.5_800","MCHTB11.5_900","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800","MCHTB11.9 S2_700","MCHTB11.9 S2_800","MCHTB11.9 S2_900"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = c(brewer.pal(name="Paired", n = 12), c("chocolate","violetred4"))
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 5)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),14)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomright", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("MCHTB11.4_700","MCHTB11.4_800","MCHTB11.4_900","MCHTB11.5_700","MCHTB11.5_800","MCHTB11.5_900","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800","MCHTB11.9 S2_700","MCHTB11.9 S2_800","MCHTB11.9 S2_900"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)


### for first-batch
y$samples$group <- c(rep("MCHTB11.9_600",16),rep("MCHTB11.9_700",16),rep("MCHTB11.9_800",8),rep("MCHTB11.9 H3_700",16),rep("MCHTB11.9 H3_800",16))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = c(brewer.pal(name="Paired", n = 3), c("chocolate","violetred4"))
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 5)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),5)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("topright", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9_800","MCHTB11.9 H3_700","MCHTB11.9 H3_800"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)

pca <- prcomp(y)
biplot(pca)
pc1 <- order(abs(pca$x[,1]), decreasing=TRUE)[1:500]
genes <- rownames(pca$x)[pc1]
ann1 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann1$ENSEMBL)
ann1 <- ann1[m,]
table(ann1$ENSEMBL==genes)

pc2 <- order(abs(pca$x[,2]), decreasing=TRUE)[1:500]
genes <- rownames(pca$x)[pc2]
ann2 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann2$ENSEMBL)
ann2 <- ann2[m,]
table(ann2$ENSEMBL==genes)

round(pc$rotation[1:10,1:2], 3)
dev.off()


###for second batch

y$samples$group <- c(rep("MCHTB11.4_600",16),rep("MCHTB11.4_700",16),rep("MCHTB11.5_600",16),rep("MCHTB11.5_700",16),rep("MCHTB11.9_600",16),rep("MCHTB11.9_700",16),rep("MCHTB11.9 H3_800",16),rep("MCHTB11.9 H3_900",16),rep("MCHTB11.9 S2_900",16))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_600","MCHTB11.4_700","MCHTB11.5_600","MCHTB11.5_700","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9 H3_800","MCHTB11.9 H3_900","MCHTB11.9 S2_900"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = brewer.pal(name="Paired", n =9)
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 15)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),9)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomright", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("MCHTB11.4_600","MCHTB11.4_700","MCHTB11.5_600","MCHTB11.5_700","MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9 H3_800","MCHTB11.9 H3_900","MCHTB11.9 S2_900"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)
dev.off()

###for second batch splitby density
y$samples$group <- c("MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9 H3_800","MCHTB11.9 H3_900")
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9 H3_800","MCHTB11.9 H3_900"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = brewer.pal(name="Paired", n =4)
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 15)
plotMDS(y, top = 400,dim=c(1,2), pch=rep(c(21,21,21,21),4)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomright", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("MCHTB11.9_600","MCHTB11.9_700","MCHTB11.9 H3_800","MCHTB11.9 H3_900"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)
dev.off()

pca <- prcomp(y)
biplot(pca)
pc1 <- order(abs(pca$x[,1]), decreasing=TRUE)[1:500]
genes <- rownames(pca$x)[pc1]
ann1 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann1$ENSEMBL)
ann1 <- ann1[m,]
table(ann1$ENSEMBL==genes)

pc2 <- order(abs(pca$x[,2]), decreasing=TRUE)[1:500]
genes <- rownames(pca$x)[pc2]
ann2 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann2$ENSEMBL)
ann2 <- ann2[m,]
table(ann2$ENSEMBL==genes)

##both first and second batch
### for first-batch
#Batch effect removal

batch = as.factor(c(rep("1",152),rep("2",128)))
batch = as.factor(c(rep("1",11),rep("2",8)))
logFC <- predFC(y,design,prior.count = 1, dispersion = 0.05)
cor(logFC[,1:6])
logCPM <- cpm(y,log =TRUE, prior.count =20)
condition <- as.factor(c(rep("mut",120),rep("wt",32),rep("mut",96),rep("wt",32)))
condition <- as.factor(c(rep("mut",9),rep("wt",2),rep("mut",6),rep("wt",2)))
design0 <- model.matrix(~0+condition)
logCPM <- removeBatchEffect(logCPM,batch =batch, design = design0)

y$samples$group <- c(rep("MCHTB11.4_700_first",16),rep("MCHTB11.4_800_first",16),rep("MCHTB11.4_900_first",8),rep("MCHTB11.5_700_first",16),rep("MCHTB11.5_800_first",16),rep("MCHTB11.5_900_first",8),rep("MCHTB11.9_600_first",16),rep("MCHTB11.9_700_first",16),rep("MCHTB11.9_800_first",8),rep("MCHTB11.9 H3_700_first",16),rep("MCHTB11.9 H3_800_first",16),rep("MCHTB11.9 S2_700_first",16),rep("MCHTB11.9 S2_800_first",16),rep("MCHTB11.9 S2_900_first",8), rep("MCHTB11.4_600_second",16),rep("MCHTB11.4_700_second",16),rep("MCHTB11.5_600_second",16),rep("MCHTB11.5_700_second",16),rep("MCHTB11.9_600_second",16),rep("MCHTB11.9_700_second",16),rep("MCHTB11.9 H3_800_second",16),rep("MCHTB11.9 H3_900_second",16),rep("MCHTB11.9 S2_900_second",16))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_700_first","MCHTB11.4_800_first","MCHTB11.4_900_first","MCHTB11.5_700_first","MCHTB11.5_800_first","MCHTB11.5_900_first","MCHTB11.9_600_first","MCHTB11.9_700_first","MCHTB11.9_800_first","MCHTB11.9 H3_700_first","MCHTB11.9 H3_800_first","MCHTB11.9 S2_700_first","MCHTB11.9 S2_800_first","MCHTB11.9 S2_900_first","MCHTB11.4_600_second","MCHTB11.4_700_second","MCHTB11.5_600_second","MCHTB11.5_700_second","MCHTB11.9_600_second","MCHTB11.9_700_second","MCHTB11.9 H3_800_second","MCHTB11.9 H3_900_second","MCHTB11.9 S2_900_second"))
y$samples$group <- c("MCHTB11.4_700_first","MCHTB11.4_800_first","MCHTB11.4_900_first","MCHTB11.5_700_first","MCHTB11.5_800_first","MCHTB11.5_900_first","MCHTB11.9_600_first","MCHTB11.9_700_first","MCHTB11.9_800_first","MCHTB11.9 H3_700_first","MCHTB11.9 H3_800_first","MCHTB11.9 S2_700_first","MCHTB11.9 S2_800_first","MCHTB11.9 S2_900_first","MCHTB11.4_600_second","MCHTB11.4_700_second","MCHTB11.5_600_second","MCHTB11.5_700_second","MCHTB11.9_600_second","MCHTB11.9_700_second","MCHTB11.9 H3_800_second","MCHTB11.9 H3_900_second","MCHTB11.9 S2_900_second")
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_700_first","MCHTB11.4_800_first","MCHTB11.4_900_first","MCHTB11.5_700_first","MCHTB11.5_800_first","MCHTB11.5_900_first","MCHTB11.9_600_first","MCHTB11.9_700_first","MCHTB11.9_800_first","MCHTB11.9 H3_700_first","MCHTB11.9 H3_800_first","MCHTB11.9 S2_700_first","MCHTB11.9 S2_800_first","MCHTB11.9 S2_900_first","MCHTB11.4_600_second","MCHTB11.4_700_second","MCHTB11.5_600_second","MCHTB11.5_700_second","MCHTB11.9_600_second","MCHTB11.9_700_second","MCHTB11.9 H3_800_second","MCHTB11.9 H3_900_second","MCHTB11.9 S2_900_second"))
y$samples$group <- c(rep("MCHTB11.4",72),rep("MCHTB11.5",72),rep("MCHTB11.9",72),rep("MCHTB11.9 H3",64))
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4","MCHTB11.5","MCHTB11.9","MCHTB11.9 H3"))
y$samples$group <- c("MCHTB11.4_700_first","MCHTB11.4_800_first","MCHTB11.4_900_first","MCHTB11.5_700_first","MCHTB11.5_800_first","MCHTB11.5_900_first","MCHTB11.9_600_first","MCHTB11.9_700_first","MCHTB11.9_800_first","MCHTB11.9 H3_700_first","MCHTB11.9 H3_800_first","MCHTB11.4_600_second","MCHTB11.4_700_second","MCHTB11.5_600_second","MCHTB11.5_700_second","MCHTB11.9_600_second","MCHTB11.9_700_second","MCHTB11.9 H3_800_second","MCHTB11.9 H3_900_second")
y$samples$group <- factor(y$samples$group, levels = c("MCHTB11.4_700_first","MCHTB11.4_800_first","MCHTB11.4_900_first","MCHTB11.5_700_first","MCHTB11.5_800_first","MCHTB11.5_900_first","MCHTB11.9_600_first","MCHTB11.9_700_first","MCHTB11.9_800_first","MCHTB11.9 H3_700_first","MCHTB11.9 H3_800_first","MCHTB11.4_600_second","MCHTB11.4_700_second","MCHTB11.5_600_second","MCHTB11.5_700_second","MCHTB11.9_600_second","MCHTB11.9_700_second","MCHTB11.9 H3_800_second","MCHTB11.9 H3_900_second"))

#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = c(brewer.pal(name="Spectral", n = 11), c(brewer.pal(name="BrBG", n = 8),c("blue")))
mycolors = c(brewer.pal(name="Spectral", n = 5))
y$samples$colour = rep(mycolors,19)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plot.pdf",width = 10,height = 15)
plotMDS(logCPM, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),23)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
plotMDS(logCPM, cex = 0.5, gene.selection = "common")
legend("bottomright", title="Condition", # << THIS IS THE HACKISH PART
       legend=levels(y$samples$group), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)
dev.off()
#prcomp(y)

pca <- prcomp(y)
biplot(pca)
pc1 <- order(abs(pca$x[,1]), decreasing=TRUE)[1:1000]
genes <- rownames(pca$x)[pc1]
ann1 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann1$ENSEMBL)
ann1 <- ann1[m,]
table(ann1$ENSEMBL==genes)

pc2 <- order(abs(pca$x[,2]), decreasing=TRUE)[1:1000]
genes <- rownames(pca$x)[pc2]
ann2 <- AnnotationDbi:::select(Homo.sapiens,keys= genes ,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(genes,ann2$ENSEMBL)
ann2 <- ann2[m,]
table(ann2$ENSEMBL==genes)

### for second-batch TPM
y$samples$group <- c(rep("patient_i3_900",16),rep("patient_i6_800",8),rep("patient_i6_850",8),rep("WT_TPM1_950",16))
y$samples$group <- factor(y$samples$group, levels = c("patient_i3_900","patient_i6_800","patient_i6_850","WT_TPM1_950"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
library(RColorBrewer)
mycolors = c(brewer.pal(name="Paired", n = 4))
y$samples$colour = rep(mycolors,4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 5)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),4)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomleft", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("patient_i3_900","patient_i6_800","patient_i6_850","WT_TPM1_950"), ncol=1, fill = mycolors, box.lwd =1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)

#prcomp(y)
dev.off()


###DEG 26.09.2024
#link:https://support.bioconductor.org/p/121073/
#link: https://support.bioconductor.org/p/86588/
#https://evangelynsim.github.io/2021_UoM_Yap_shRNA_nuclei_RNAseq_ATACseq/03.RNAseq_EdgeR_and_ScatterPlot.html
condition <- as.factor(c(rep("mut",40),rep("wt",32),rep("mut",32),rep("wt",32)))
batch = as.factor(c(rep("1",72),rep("2",64)))
design0 <- model.matrix(~0+condition+batch)
rownames(design0) <- colnames(y)


data.int <- readRDS("/Volumes/CARD3/Neda/20240307-sciRNAseq3-7th-TTNtv/output/combo/int-vitro-TTNhom-RPMI.Rds")
col = c("CM"="#7570B3","Fib"="#E6AB02", "CM(Prlf)" = "#CF9FFF")
DimPlot(data.int, reduction = "umap",label=FALSE,label.size = 6,pt.size = 0.5,split.by = "orig.ident", cols = col)
FeaturePlot(data.int,features = "ENSG00000107796", split.by = "orig.ident")
VlnPlot(data.int,features = "ENSG00000107796", split.by = "orig.ident")
data.int$Broad_celltype <- factor(data.int$Broad_celltype, levels = c("CM","CM(Prlf)","Fib"))
DefaultAssay(data.int) <- "RNA"
cm <- subset(data.int,subset = Broad_celltype=="CM")
cm$pt <- paste0(cm$orig.ident,cm$density)
cm$ptrep <- paste0(cm$pt,cm$replicate)
cm.subset <- cm
dim(cm.subset)
#9923 3290

all <- readRDS("/Volumes/CARD3/Neda/20240307-sciRNAseq3-7th-TTNtv/output/combo/all15098-seurat.Rds")
#all$replicate <- factor(all$replicate, levels = paste("rep",8, sep = ""))
index <- match(colnames(cm.subset), colnames(all))
index <- na.omit(index)
all <- all[,index]

#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(org.Hs.eg.db,keys=rownames(all),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","CHR"),keytype = "ENSEMBL")
m <- match(rownames(all),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(all))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 232
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 200
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1]  690
  
# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
#1044
all <- all[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(all))
#numzero.genes <- rowSums(all@assays$RNA@counts==0)
numzero.genes <- apply(all@assays$RNA@counts, 1, function(x) sum(x == 0))
table(numzero.genes < (ncol(all)-20))
#FALSE  TRUE 
#9779  9810 
keep.genes <- numzero.genes < (ncol(all)-20)
all <- all[keep.genes,]
dim(all)
ann <- ann[keep.genes,]
  
# Take out sex chromosome genes
xy <- ann$CHR %in% c("X","Y")
table(xy)

all <- all[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(all))

shared.genes <- Reduce(intersect, list(rownames(all), ann$ENSEMBL, rownames(cm.subset)))
cm.subset <- cm.subset[shared.genes,]
all <- all[shared.genes,]
ann <- ann[ann$ENSEMBL %in% shared.genes,]
table(rownames(cm.subset)==rownames(all))
table(rownames(cm.subset)==ann$ENSEMBL)
table(rownames(all)==ann$ENSEMBL)
#TRUE 
#9453 

broadct <- "CM"
cm.subset$biorep <- "NA"
index <- which(cm.subset$pt =="patient_hom.59500")
index <- na.omit(index)
#cm.subset$biorep[index] <- paste("pthet400",sep = "")
cm.subset$biorep[index] <- paste("pthom59500",cm.subset$replicate[index],sep = "")
index <- which(cm.subset$pt =="WT_TTN200")
index <- na.omit(index)
#cm.subset$biorep[index] <- paste("cor900",sep = "")
cm.subset$biorep[index] <- paste("cor",cm.subset$replicate[index],sep = "")

#batch 1
sam <- factor(cm.subset$biorep,levels=c(paste(rep("pthom42rep",8),1:8,sep = ""),paste(rep("pthom59400rep",8),1:8,sep = ""),paste(rep("pthom59500rep",8),1:8,sep = ""),paste(rep("correp",8),1:8,sep = "")))

sam <- factor(cm.subset$biorep,levels=c(paste(rep("pthet400rep",8),1:8,sep = ""),paste(rep("pthet500rep",8),1:8,sep = ""),paste(rep("correp",8),1:8,sep = "")))

newgrp <- paste(broadct,sam,sep=".")
newgrp <- factor(newgrp,levels=paste("CM",levels(sam),sep="."))
table(newgrp)

des <- model.matrix(~0+newgrp)
colnames(des) <- levels(newgrp)
dim(des)
#3290   24
dim(all)
#9453 3290
pb <- all@assays$RNA@counts %*% des
y.pb <- DGEList(pb)
keep <- rowSums(cpm(y.pb)>1)>=8  #median(y$samples$lib.size) is ~231,091. A CPM of 43 is used as it corresponds to a count of 10 for this median library size (~231091) in this data set. 8 = minimum number of replicates per condition.
table(keep)          #43 batch 1 TPM1 MM CM
#keep                
#TRUE 
#9453
y.pb <- y.pb[keep,]
ann <- ann[keep,]

saminfo <- matrix(unlist(strsplit(colnames(y.pb$counts),split="[.]")),ncol=2,byrow=TRUE)
bct <- factor(saminfo[,1])
indiv <- factor(saminfo[,2])
group <- rep(NA,ncol(y.pb))
group[grep("cor",indiv)] <- "cor"
group[grep("pt",indiv)] <- "pt"
group <- factor(group,levels=c("cor","pt"))
table(rownames(y.pb)==ann$ENSEMBL)
#TRUE 
#9453 

y.pb$genes <- ann
bct2 <- as.character(bct)
bct2[bct2 == "CM"] <- "CM"
bct2[bct2 == "Fib"] <- "Fib"



newgrp <- paste(bct2,group,sep=".")
newgrp <- factor(newgrp)

design <- model.matrix(~0+newgrp)
colnames(design) <- levels(newgrp)

par(mfrow=c(1,1))
y.pb <- calcNormFactors(y.pb)
plotMDS.DGEList(y.pb,gene.selection="common",dim=c(1,2))  #,pch=rep(c(21,21,21,21),4)[factor(y.pb$samples$group)]
v <- voom(y.pb,design,plot=TRUE)
fit <- lmFit(v,design)
cont.cardio <- makeContrasts(ptvscorCM = CM.pt-CM.cor,
                             #ptvscorFib = Fib.pt-Fib.cor, 
                             #ptvscor = 0.5*(CM.pt+Fib.pt)-0.5*(CM.cor+Fib.cor),
                             levels=design)
fit.cardio <- contrasts.fit(fit,contrasts = cont.cardio)
fit.cardio <- eBayes(fit.cardio,robust=TRUE)

summary(decideTests(fit.cardio))
treat.cardio <- treat(fit.cardio,lfc=0)

dt.cardio<-decideTests(treat.cardio)

summary(dt.cardio)

options(digits=3)
topTreat(treat.cardio,coef=1,n=20,p.value=0.05)[,-c(1:3)]

fdr <- apply(treat.cardio$p.value, 2, function(x) p.adjust(x, method="BH"))
output <- data.frame(treat.cardio$genes,LogFC=treat.cardio$coefficients,AveExp=treat.cardio$Amean,tstat=treat.cardio$t, pvalue=treat.cardio$p.value, fdr=fdr)
library(dplyr)
colnames(output)[6] <- "LogFC"
colnames(output)[8] <- "tstat"
colnames(output)[9] <- "p.value"
colnames(output)[10] <- "fdr"

write.csv(output, "~/Desktop/paper-layout/Fig3/CM-TTNhom-MM.csv")

VlnPlot(cm.subset, features = c("ENSG00000133107"), split.by = "pt")


output$rank <- sign(output$LogFC)*-1*log10(output$fdr)
output <- output[order(output$rank,decreasing = T),]
write.csv(output, "../../../../sci-RNA-seq3-7th-TTNtv/output/combo/docs/CM-TTNhom-MM-8.csv")
output <- read.csv("../../../../sci-RNA-seq3-7th-TTNtv/output/combo/docs/CM-TTNhet-RPMI-8.csv")
index <- which(output$fdr < 0.05 )
index2 <- which(output$fdr < 0.05 & output$LogFC>0)
output <- output[index2,]
View(output[index,])
View(output.first[index,])
shared.up <- intersect(output.first$SYMBOL[index],output$SYMBOL[index2])
shared.dn <- intersect(output.first$SYMBOL[index],output$SYMBOL[index2])
shared.ns <- intersect(output.first$SYMBOL[index],output$SYMBOL[index2])
shared <- c(shared.up,shared.dn,shared.ns)
output.shared <- output.first[output.first$SYMBOL %in% shared,]
write.csv(output.shared, "../../../../../../../Fib-both-updnns.csv")
library(ggrepel)

output.shared$col <- NA
index <- which(output.shared$p_value < 0.05 & output.shared$score>0)
output.shared$col[index] <- "up"
index <- which(output.shared$p_value < 0.05 & output.shared$score < 0)
output.shared$col[index] <- "down"
#output.shared <- output.shared[order(output.shared$rank,decreasing = T),]
index <- which((output.shared$score > 3 | output.shared$score < -3) & output.shared$p_value < 0.05)
output.shared$labels <- NA
output.shared$labels[index] <- TRUE
pdf(file="../../../../../../../ggplotnew4.pdf",width = 10,height = 5)
ggplot(data=output.shared, aes(x=score, y=-log10(p_value), colour = col, label=ifelse(labels == TRUE, as.character(source),""))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_classic()
dev.off()

fgp <- read.csv("/group/card2/Neda/MCRI_LAB/Publications_MCRI_Enzo Porrello & David Elliott/DCM-fetal-gene-program-manuscript/FGP-CMup.csv")
index <- which(output$fdr < 0.05 & output$LogFC>0)
outputv2 <- output[index,] 
index <- match(fgp$fgpUP,outputv2$SYMBOL)
index <- na.omit(index)
outputv2$SYMBOL[index]

load("/group/card2/Neda/MCRI_LAB/scRNAseq-ES/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.c2,treat.cardio$genes$ENTREZID)
reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
cardio.camera <- cameraPR(treat.cardio$t[,1], reactome.id)
cardio.camera$f <- -log10(cardio.camera$FDR)
cardio.camera.up <- cardio.camera[cardio.camera[,2]=="Up",]
library(dplyr)
cardio.camera.up <- cardio.camera.up %>% filter(cardio.camera.up$FDR < 0.05)

cardio.camera.dn <- cardio.camera[cardio.camera[,2]=="Down",]
library(dplyr)
cardio.camera.dn <- cardio.camera.dn %>% filter(cardio.camera.dn$FDR < 0.05)

cardio.camera <- cardio.camera %>% filter(cardio.camera$FDR < 0.05)
write.csv(cardio.camera,"/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-6th-TPM1-DSPKI/output/combo/docs/DSP-fib-batch2-reactom.csv")

index <- match(rownames(cardio.camera.up), names(reactome.id))
names(reactome.id[index])
entrezID <- reactome.id[["REACTOME_GLYCOLYSIS"]]
index <- match(entrezID,output$ENTREZID)
index <- na.omit(index)
output[index,]

load("/group/card2/Neda/MCRI_LAB/scRNAseq-ES/RDataObjects/human_c5_v5p2.rdata")
c5.id <- ids2indices(Hs.c5,output$ENTREZID)
go.cardio <- cameraPR(treat.cardio$t[,1],c5.id)
go.cardio$names <- rownames(go.cardio)
go.cardio.up <- go.cardio[go.cardio[,2]=="Up",]
library(dplyr)
go.cardio.up <- go.cardio.up %>% filter(go.cardio.up$FDR < 0.05)
go.cardio.dn <- go.cardio[go.cardio[,2]=="Down",]
library(dplyr)
go.cardio.dn <- go.cardio.dn %>% filter(go.cardio.dn$FDR < 0.05)


#goSeq
genes <- output.shared$fdr <.05 & output.shared$LogFC<0
names(genes)=output.shared$ENTREZID
print(head(genes))
library(goseq)
pwf <- nullp(genes, "hg19","knownGene")
go.results <- goseq(pwf, "hg19","knownGene", test.cats=c("GO:BP"),use_genes_without_cat=TRUE)
go.results$adjP <- p.adjust(go.results$over_represented_pvalue, method="BH")
go.results %>% 
  top_n(20, wt=-adjP) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=adjP, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="Adj.p value", size="Count")
go <- go.results %>% filter(go.results$adjP < 0.05)


  
#https://github.com/stephenturner/msigdf/blob/master/data-raw/msigdf.R
#used this code for fetal gene program Ch4 22.05
output.shared <- read.csv("../../../../../../../CM-both-updnns.csv")
fetal <- read.csv("../../../../../../../155FGP-CM.csv")

index <- match(fetal$total_genes[fetal$dir == "Up-regulated"],output.shared$SYMBOL[output.shared$LogFC>0 & output.shared$fdr < 0.05])
index <- na.omit(index)
shared <- output.shared[index,]
DB <- paste("org", "Hs", "eg", "db", sep = ".")
require(DB, character.only = TRUE)
GO2ALLEGS <- paste("org", "Hs", "egGO2ALLEGS", sep = ".")
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]

index <- which(shared$LogFC > 0 & shared$fdr < 0.05)
sig <- shared[index,]
topgo <- goana(de=sig$ENTREZID,universe=output.shared$ENTREZID,species="Hs")
topgo$adjP <- p.adjust(topgo$P.DE, method="BH")
topgo <- topgo %>% filter(topgo$P.DE < 0.05 & topgo$Ont=="BP") 
#topgo <- topgo %>% filter(topgo$P.DE < 0.05) 
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))
de.by.go <- lapply(de.by.go, FUN=function(x) { x[x %in% sig$ENTREZID] })
result <- data.frame(matrix(NA,nrow = dim(topgo)[1],ncol = 5))
for (i in 1:dim(topgo)[1]){
  genename <- sig[sig$ENTREZID %in% de.by.go[[paste(rownames(topgo)[i],topgo$Ont[i],sep = ".")]],]
  result[i,1] <- topgo$Term[i]
  result[i,2] <- topgo$Ont[i]
  result[i,3] <- topgo$P.DE[i]
  result[i,4] <- topgo$DE[i]
  result[i,5] <- topgo$N[i]
  result[i,6] <- topgo$DE[i]*100/topgo$N[i]
  vec <- c(genename$SYMBOL)
  fvec <- shQuote(vec, type = "cmd")
  comma_vec <- paste(fvec, collapse = ", ")
  result[i,7] <- comma_vec
}
colnames(result) <- c("Term","Ont","P","noDE","N","hits","DE")
rownames(result) <- rownames(topgo)


#used the following goana code for GO terms 21.05 for ch4
index <- which(output.shared$LogFC > 0 & output.shared$fdr < 0.05)
sig <- output.shared[index,]
topgo <- goana(de=sig$ENTREZID,universe=output.shared$ENTREZID,species="Hs")
topgo$adjP <- p.adjust(topgo$P.DE, method="BH")
topgo <- topgo %>% filter(topgo$adjP < 0.05 & topgo$Ont=="BP") 
#topgo <- topgo %>% filter(topgo$P.DE < 0.05) 
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))
de.by.go <- lapply(de.by.go, FUN=function(x) { x[x %in% sig$ENTREZID] })
result <- data.frame(matrix(NA,nrow = dim(topgo)[1],ncol = 5))
for (i in 1:dim(topgo)[1]){
  genename <- sig[sig$ENTREZID %in% de.by.go[[paste(rownames(topgo)[i],topgo$Ont[i],sep = ".")]],]
  result[i,1] <- topgo$Term[i]
  result[i,2] <- topgo$Ont[i]
  result[i,3] <- topgo$adjP[i]
  result[i,4] <- topgo$DE[i]
  result[i,5] <- topgo$N[i]
  result[i,6] <- topgo$DE[i]*100/topgo$N[i]
  vec <- c(genename$SYMBOL)
  fvec <- shQuote(vec, type = "cmd")
  comma_vec <- paste(fvec, collapse = ", ")
  result[i,7] <- comma_vec
}
colnames(result) <- c("Term","Ont","adjP","noDE","N","hits","DE")
rownames(result) <- rownames(topgo)
write.csv(result,"/group/card2/Neda/dn.csv")

#downregulated
result %>% 
  top_n(40, wt=-adjP) %>% 
  mutate(hitsPerc=noDE*100/N) %>%
  ggplot(aes(x=hitsPerc, 
             y=Term, 
             colour=adjP, 
             size=noDE)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="Adj.p value", size="Count")


first <- read.csv(file = "../../../../../../../results/CM-firstBatch.csv")
second <- read.csv(file = "../../../../../../../results/CM-secondBatch.csv")

shared.up <- intersect(first$ENSEMBL[first$ptvscorCM>0],second$ENSEMBL[second$ptvscorCM>0])
shared.dn <- intersect(first$ENSEMBL[first$ptvscorCM<0],second$ENSEMBL[second$ptvscorCM<0])

sig.up <- first[first$ENSEMBL %in% shared.up,]
sig.dn <- first[first$ENSEMBL %in% shared.dn,]


sig.merge <- c(shared.up,shared.dn)
sig <- first[first$ENSEMBL %in% sig.merge,]
DB <- paste("org", "Hs", "eg", "db", sep = ".")
require(DB, character.only = TRUE)
GO2ALLEGS <- paste("org", "Hs", "egGO2ALLEGS", sep = ".")
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]

sig <- sig %>% filter(sig$ptvscorCM < -1.5)
topgo <- topGO(goana(de=rankedlist$ENTREZID,species="Hs"),number = 1000)
topgo <- topgo %>% filter(topgo$P.DE < 0.05 & topgo$Ont=="BP") 
topgo <- topgo %>% filter(topgo$P.DE < 0.05) 
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))
de.by.go <- lapply(de.by.go, FUN=function(x) { x[x %in% sig$ENTREZID] })
result <- data.frame(matrix(NA,nrow = dim(topgo)[1],ncol = 5))
for (i in 1:dim(topgo)[1]){
  genename <- sig[sig$ENTREZID %in% de.by.go[[paste(rownames(topgo)[i],topgo$Ont[i],sep = ".")]],]
  result[i,1] <- topgo$Term[i]
  result[i,2] <- topgo$Ont[i]
  result[i,3] <- topgo$P.DE[i]
  result[i,4] <- topgo$DE[i]
  vec <- c(genename$SYMBOL)
  fvec <- shQuote(vec, type = "cmd")
  comma_vec <- paste(fvec, collapse = ", ")
  result[i,5] <- comma_vec
}
colnames(result) <- c("Term","Ont","pvalue","no. DE","DE")
rownames(result) <- rownames(topgo)
write.csv(result,"/group/card2/Neda/dn-shared-noTPM.csv")

#genename <- sig[sig$ENTREZID %in% de.by.go[["GO:0030217.BP"]],]
#index <- match(genename$SYMBOL, output$SYMBOL)
#index <- na.omit(index)
#View(output[index,])
output$SYMBOL[index]
#index <- match(genename$SYMBOL, shared)
#index <- na.omit(index)
#shared[index]


load("../../../../../../scRNAseq-ES/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.c2,rankedlist$ENTREZID)
reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
cardio.camera <- cameraPR(sig$ptvscorCM.1,reactome.id)
cardio.camera.up <- cardio.camera[cardio.camera[,2]=="Up",]
library(dplyr)
cardio.camera.up <- cardio.camera.up %>% filter(cardio.camera.up$FDR < 0.05)

cardio.camera.dn <- cardio.camera[cardio.camera[,2]=="Down",]
library(dplyr)
cardio.camera.dn <- cardio.camera.dn %>% filter(cardio.camera.dn$FDR < 0.05)

index <- match(rownames(cardio.camera.dn), names(reactome.id))
write.csv(cardio.camera.dn,"/group/card2/Neda/REACTOME-updn.csv")
genes <- first[reactome.id[[index[6]]],]
##link:https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#load-pathways
library(fgsea)

ranks <- output.shared$rank
names(ranks) <- output.shared$ENTREZID
head(ranks)
barplot(sort(ranks, decreasing = T))

load("../../../../../../scRNAseq-ES/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.H,output.shared$ENTREZID)
reactome.id <-c2.id[grep("HALLMARK",names(c2.id))]
#l <- list()
#for ( i in 1:length(reactome.id)){
#  if(length(reactome.id[i][[1]]) > 15 & length(reactome.id[i][[1]]) < 500)
#  l <- c(l,reactome.id[i])
#}
fgseaRes <- fgsea(reactome.id, ranks)
pval <- data.frame(fgseaRes$pval)
rownames(pval) <- rownames(fgseaRes$pathway)
fdr <- apply(pval, 2, function(x) p.adjust(x, method="BH"))
#head(fgseaRes[order(padj, -abs(NES)), ], n=10)
fgseaRes <- fgseaRes %>% filter(fgseaRes$pval < 0.05)

##make a jitter plot
library(ggplot2)
library(ggrepel)
sciplex.markers$sig <- ifelse((sciplex.markers$p_val_adj < 0.05 & sciplex.markers$avg_log2FC > 3), sciplex.markers$`ann$SYMBOL[index]`, "")

pdf(file="../../../../../../../ggplot_jitter_markers_iPSC.pdf",width = 10,height = 5)
ggplot(sciplex.markers, aes(x = cluster, y = avg_log2FC, label = sig)) + 
  geom_jitter(aes(color = ifelse(p_val_adj < 0.05, "red", "black")), position = position_jitter(seed = 1))+ scale_color_identity()+
  xlab("Cluster") + ylab("Average Log2FC") + geom_text(aes(label = sig), position = position_jitter(seed = 1))+
  theme_classic()
dev.off()

all$biorep <- paste(all$orig.ident,all$density,sep = "_")
cm.split <- SplitObject(all,split.by = "biorep" )
#Calculating statistics
libmode <- rep(NA,14)
numgmode <- rep(NA,14)
numcells <- rep(NA,14)

names(libmode) <- names(numgmode) <- names(numcells) <- names(cm.split)

mylibs <- mynumg <- vector("list", 14) 

for ( i in 1:length(cm.split)){
  mylibs[[i]] <- colSums(cm.split[[i]])
  mynumg[[i]] <- colSums(cm.split[[i]]@assays$RNA@counts!=0)
  numcells[i] <- ncol(cm.split[[i]])
  libmode[i] <- find_modes(mylibs[[i]])[1]
  numgmode[i] <- find_modes(mynumg[[i]])[1]
}

#Number of cells per well
par(mar=c(6,4,2,2))
mycols <- rep(ggplotColors(3),rep(16,144))
barplot(numcells,col=mycols,las=2,main="Number of nuclei per well", ylim = c(0,500), cex.names = 0.5)
abline(h=100,lty=2,lwd=2)
#Library size distributions (sequencing depth)
par(mfrow=c(1,1))
par(mar=c(8,4,3,2))
boxplot(mylibs,names=names(libmode),las=2, main="Library size distributions-1st batch",ylim=c(0,1800))

# no genes detected 
par(mfrow=c(1,1))
par(mar=c(8,4,3,2))
boxplot(mynumg,names=names(libmode),las=2, main="Number of genes detected-1st batch", ylim=c(1,1000))


####end of 08.04
heartexplorer <- read.csv("../../../../../../../pseudo_in vivo.csv")
output <- read.csv("../../../../../../../outputpseudo-in vitro.csv")
heartexplorer <- heartexplorer %>% filter(heartexplorer$fdr.DNfib < 0.05)
CM <- output %>% filter(output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM > 0)
Fib <- output %>% filter(output$fdr.ptvscorFib < 0.05 & output$LogFC.ptvscorFib > 0)
pseudo <- output %>% filter(output$fdr.ptvscor < 0.05 & output$LogFC.ptvscor > 0)
shared<- intersect(intersect(CM$SYMBOL, Fib$SYMBOL), pseudo$SYMBOL)
#"KIF26B" "MYBPC3"
View(Fib)
CM <- output %>% filter(output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM < 0)
Fib <- output %>% filter(output$fdr.ptvscorFib < 0.05 & output$LogFC.ptvscorFib < 0)
pseudo <- output %>% filter(output$fdr.ptvscor < 0.05 & output$LogFC.ptvscor < 0)
shared<- intersect(intersect(CM$SYMBOL, Fib$SYMBOL), pseudo$SYMBOL)

shared<- intersect(heartexplorer$SYMBOL, output$SYMBOL[which(output$fdr.ptvscorFib < 0.05)])

library(ggrepel)
output$sig <- ifelse((output$SYMBOL %in% shared) , output$SYMBOL, "")
output$col <- NA
index <- which(output$LogFC.ptvscor > 0 & output$fdr.ptvscor < 0.05)
output$col[index] <- "up"
index <- which(output$LogFC.ptvscor < 0 & output$fdr.ptvscor < 0.05)
output$col[index] <- "down"
ggplot(data=output, aes(x=LogFC.ptvscor, y=-log10(fdr.ptvscor), colour = col, label=sig)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20000) +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_classic()




## compare patient-CM with patient_corrected_CM
# Read in data (Seurat object) generated in "./Notebook3_preprocess_data.R"
sciplex <- readRDS("../../../output/sciplex-annotatedCellType.Rds")

pt.ptcorrected <- merge(x = sciplex[[1]],y = sciplex[[4]])
cardio <- subset(pt.ptcorrected, subset = Broad_celltype=="CM")
DefaultAssay(cardio) <- "RNA"
par(mar=c(4,4,2,1))
plot(density(cardio$nFeature_RNA),main="Number of genes detected")
abline(v=200,col=4, lty=3)
legend("topright",lty=2,col=4,legend="#genes = 200")

plot(density(cardio$nCount_RNA),main="Library size")
abline(v=300,col=4, lty=3)
legend("topright",lty=2,col=4,legend="library size = 300")

cardio.list <- SplitObject(cardio, split.by = "orig.ident")
min <- min(sapply(cardio.list, ncol))
for (i in 1:length(cardio.list)) {
  cardio.list[[i]] <- SCTransform(cardio.list[[i]], verbose = FALSE)
}
cardio.anchors <- FindIntegrationAnchors(object.list = cardio.list, dims=1:30,anchor.features = 3000,k.filter=min)
cardio.integrated <- IntegrateData(anchorset = cardio.anchors,dims=1:30)
DefaultAssay(object = cardio.integrated) <- "integrated"
cardio.integrated <- ScaleData(cardio.integrated, verbose = FALSE)
cardio.integrated <- RunPCA(cardio.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(cardio.integrated,ndims=50)
VizDimLoadings(cardio.integrated, dims = 1:4, reduction = "pca")
DimPlot(cardio.integrated, reduction = "pca",group.by="orig.ident")
DimHeatmap(cardio.integrated, dims = 1:15, cells = 500, balanced = TRUE)
cardio.integrated <- FindNeighbors(cardio.integrated, dims = 1:20)
cardio.integrated <- FindClusters(cardio.integrated, resolution = 0.9)
table(Idents(cardio.integrated))
set.seed(10)
cardio.integrated <- RunUMAP(cardio.integrated, reduction = "pca", dims = 1:10)
DimPlot(cardio.integrated, reduction = "umap",label=TRUE,label.size = 6,pt.size = 0.5)+NoLegend()


### DEG 
sciplex <- readRDS("../../../output/sciplex-annotatedCellType.Rds")
all <- list[[1]]@assays$RNA@counts
all <- cbind(all,list[[4]]@assays$RNA@counts)
#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(all),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(all))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 229
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 1366

# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
table(ann$ENSEMBL==rownames(all))
#TRUE 
#20633

#remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
#num_cells_expressed <- rowSums(all !=0)
numzero.genes <- rowSums(all ==0)
#fraction_cells_expressed <- (num_cells_expressed/ncol(all))
#keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
table(numzero.genes > (ncol(all)-20))
keep.genes <- numzero.genes < (ncol(all)-20)

table(keep.genes)
all <- all[keep.genes,]
dim(all)
#7530  6614
ann <- ann[keep.genes,]

#remove low quality nuclei
keep.nuclei <- colSums(all !=0 ) >= 200
all  <- all[,keep.nuclei]

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
xy
#FALSE  TRUE 
#16566   621
all<- all[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(all))
#TRUE 
#16566 
sciplex.sub <- merge(x = sciplex[[1]],y = sciplex[[4]])
sciplex.sub <- subset(sciplex.sub, subset=Broad_celltype != "NA")
index <- intersect(colnames(all), colnames(sciplex.sub))
sciplex.sub <- sciplex.sub[,index]
all <- all[,index]
broadct <- factor(sciplex.sub$Broad_celltype)
index <- which(sciplex.sub$orig.ident =="patient.corrected")
sciplex.sub$biorep <- "NA"
sciplex.sub$biorep[index] <- paste("cor",sciplex.sub$replicate[index],sep = "")
index <- which(sciplex.sub$orig.ident =="patient")
sciplex.sub$biorep[index] <- paste("pt",sciplex.sub$replicate[index],sep = "")
sam <- factor(sciplex.sub$biorep,levels=c("correp1","correp2","correp3","correp4","correp5","correp6","correp7","correp8","ptrep1","ptrep2","ptrep3","ptrep4","ptrep5","ptrep6","ptrep7","ptrep8"))
newgrp <- paste(broadct,sam,sep=".")
newgrp <- factor(newgrp,levels=paste(rep(levels(broadct),each=16),levels(sam),sep="."))
table(newgrp)

des <- model.matrix(~0+newgrp)
colnames(des) <- levels(newgrp)
dim(des)
#1660   32
dim(all)
#7297 1781
pb <- all %*% des
y.pb <- DGEList(pb)

saminfo <- matrix(unlist(strsplit(colnames(y.pb$counts),split="[.]")),ncol=2,byrow=TRUE)
bct <- factor(saminfo[,1])
indiv <- factor(saminfo[,2])
group <- rep(NA,ncol(y.pb))
group[grep("cor",indiv)] <- "cor"
group[grep("pt",indiv)] <- "pt"
group <- factor(group,levels=c("cor","pt"))
table(rownames(y.pb)==ann$ENSEMBL)
#TRUE 
#7297

y.pb$genes <- ann
bct2 <- as.character(bct)
bct2[bct2 == "CM"] <- "CM"
bct2[bct2 == "Fib"] <- "Fib"



newgrp <- paste(bct2,group,sep=".")
newgrp <- factor(newgrp)

design <- model.matrix(~0+newgrp)
colnames(design) <- levels(newgrp)

par(mfrow=c(2,2))
y.pb <- calcNormFactors(y.pb)
v <- voom(y.pb,design,plot=TRUE)

fit <- lmFit(v,design)
cont.cardio <- makeContrasts(ptvscorCM = CM.pt-CM.cor,
                             ptvscorFib = Fib.pt-Fib.cor, 
                             ptvscor = 0.5*(CM.pt+Fib.pt)-0.5*(CM.cor+Fib.cor),
                             levels=design)
fit.cardio <- contrasts.fit(fit,contrasts = cont.cardio)
fit.cardio <- eBayes(fit.cardio,robust=TRUE)

summary(decideTests(fit.cardio))
treat.cardio <- treat(fit.cardio,lfc=0)

dt.cardio<-decideTests(treat.cardio)

summary(dt.cardio)

options(digits=3)
topTreat(treat.cardio,coef=3,n=20,p.value=0.05)[,-c(1:3)]

fdr <- apply(treat.cardio$p.value, 2, function(x) p.adjust(x, method="BH"))
output <- data.frame(treat.cardio$genes,LogFC=treat.cardio$coefficients,AveExp=treat.cardio$Amean,tstat=treat.cardio$t, pvalue=treat.cardio$p.value, fdr=fdr)
library(dplyr)
#output <- output %>% filter(output$fdr.corvsptCM < 0.05)
write.csv(output[index,], "../../../../../../../Fib.Patient-PatientCor.csv")
hm <- c("HDAC9","FLNC","SORBS2","ANK2","WT1","FGF9","RBM24","SIX1","CDH2","RGS4","MYEF2","FLNB","TLL2","DDR2","THY1","COL1A1","COL1A2")

index <- match(hm,shared)
index <- na.omit(index)
hooray <- output[index,]


heartexplorer <- read.csv("../../../../../../../pseudo_in vivo.csv")
output <- read.csv("../../../../../../../outputpseudo-in vitro.csv")
heartexplorer <- heartexplorer %>% filter(heartexplorer$fdr.DNfib < 0.05)
CM <- output %>% filter(output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM > 0)
Fib <- output %>% filter(output$fdr.ptvscorFib < 0.05 & output$LogFC.ptvscorFib > 0)
pseudo <- output %>% filter(output$fdr.ptvscor < 0.05 & output$LogFC.ptvscor > 0)
shared<- intersect(intersect(CM$SYMBOL, Fib$SYMBOL), pseudo$SYMBOL)
#"KIF26B" "MYBPC3"
View(Fib)
CM <- output %>% filter(output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM < 0)
Fib <- output %>% filter(output$fdr.ptvscorFib < 0.05 & output$LogFC.ptvscorFib < 0)
pseudo <- output %>% filter(output$fdr.ptvscor < 0.05 & output$LogFC.ptvscor < 0)
shared<- intersect(intersect(CM$SYMBOL, Fib$SYMBOL), pseudo$SYMBOL)

shared<- intersect(heartexplorer$SYMBOL, output$SYMBOL[which(output$fdr.ptvscorFib < 0.05)])

library(ggrepel)
output$sig <- ifelse((output$SYMBOL %in% shared) , output$SYMBOL, "")
output$col <- NA
index <- which(output$LogFC.ptvscor > 0 & output$fdr.ptvscor < 0.05)
output$col[index] <- "up"
index <- which(output$LogFC.ptvscor < 0 & output$fdr.ptvscor < 0.05)
output$col[index] <- "down"
ggplot(data=output, aes(x=LogFC.ptvscor, y=-log10(fdr.ptvscor), colour = col, label=sig)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20000) +
  scale_color_manual(values=c("blue", "red", "black")) +
  #geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+ theme_classic()

load("../../../../../../scRNAseq-ES/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.c2,treat.cardio$genes$ENTREZID)
reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
cardio.camera <- cameraPR(treat.cardio$t[,c(1,2)],reactome.id)
cardio.camera.up <- cardio.camera[cardio.camera[,1]=="Up",]
library(dplyr)
cardio.camera.up <- cardio.camera.up %>% filter(cardio.camera.up$FDR < 0.05)

cardio.camera.dn <- cardio.camera[cardio.camera[,1]=="Down",]
library(dplyr)
cardio.camera.dn <- cardio.camera.dn %>% filter(cardio.camera.dn$FDR < 0.05)

fib.camera <- cameraPR(treat.cardio$t[,2],reactome.id)
fib.camera.up <- fib.camera[fib.camera[,2]=="Up",]
library(dplyr)
fib.camera.up <- fib.camera.up %>% filter(fib.camera.up$FDR < 0.05)

fib.camera.dn <- fib.camera[fib.camera[,2]=="Down",]
library(dplyr)
fib.camera.dn <- fib.camera.dn %>% filter(fib.camera.dn$FDR < 0.05)

DB <- paste("org", "Hs", "eg", "db", sep = ".")
require(DB, character.only = TRUE)
GO2ALLEGS <- paste("org", "Hs", "egGO2ALLEGS", sep = ".")
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]

#for CM index <- which((output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM > 1.5) | (output$fdr.ptvscorCM < 0.05 & output$LogFC.ptvscorCM < -1.5))
index <- which(output$fdr.ptvscor < 0.05)
sig <- output[index,]
topgo <- topGO(goana(de=sig$ENTREZID,universe=output$ENTREZID,species="Hs"),number = 1000)
topgo <- topgo %>% filter(topgo$P.DE < 0.05 & topgo$Ont=="BP") 
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))
de.by.go <- lapply(de.by.go, FUN=function(x) { x[x %in% sig$ENTREZID] })
genename <- sig[sig$ENTREZID %in% de.by.go[["GO:0070371.BP"]],]
index <- match(genename$SYMBOL, output$SYMBOL)
index <- na.omit(index)
View(output[index,])
output$SYMBOL[index]
#index <- match(genename$SYMBOL, shared)
#index <- na.omit(index)
#shared[index]


data <- output[index,]
data1 <- data[order(data$LogFC.ptvscor,decreasing = T),]
data1 <- data1 %>% filter(data1$LogFC.ptvscor > 1.5)
data2 <- data[order(data$LogFC.ptvscor,decreasing = F),]
data2 <- data2 %>% filter(data2$LogFC.ptvscor < -1.5)
data3<- rbind(data1,data2)
#data.merge <- data3
data.merge<- rbind(data.merge,data3)


index <- match(data.merge$ENSEMBL, rownames(v))
index <- na.omit(index)


sumexpr <- matrix(NA, nrow= length(index), ncol=16)
new <- v$E[index,c(1:32)]
for ( j in 1:length(index)){
  for (i in 1:16){
      sumexpr[j,i] <- sum(new[j,i],new[j,(i+16)])/2
  }
}
i <- data.merge$SYMBOL[match(rownames(v)[index], data.merge$ENSEMBL)]
library(NMF)
par(mfrow=c(1,1))
#aheatmap(v$E[index,c(1:32)],Rowv = NA,Colv = NA, labRow = i,
aheatmap(sumexpr,Rowv = NA,Colv = NA, labRow = i,
         fontsize=5,color="-RdYlBu",cexRow =1, cexCol = 1,
         scale="row")

## compare CM with CM-4F
sciplex.markers <- FindMarkers(sciplex[[2]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, ident.1="0", ident.2 = "1")
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(sciplex.markers),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(sciplex.markers),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(sciplex.markers))
index <- match(rownames(sciplex.markers),ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])
write.csv(sciplex.markers,"../../../../../../../up-CMvsCM4F.csv")



sub.sciplex <- list()
for (i in 1:2){
sub.sciplex[[i]] <- sciplex[[i+1]]
names(sub.sciplex)[[i]] <- names(sciplex)[[i+1]]
}
sub.sciplex[[1]] <- RunUMAP(sub.sciplex[[1]], reduction = "pca", dims = 1:15)
DimPlot(sub.sciplex[[1]], reduction = "umap",label=T,label.size = 5)
sub.sciplex[[2]]$Broad_celltype <- "NA"
index <- which(sub.sciplex[[2]]$RNA_snn_res.0.1=="2")
sub.sciplex[[2]]$Broad_celltype[index] <- "Unknown"

features <- SelectIntegrationFeatures(object.list = sub.sciplex)

sciplex.anchors <- FindIntegrationAnchors(object.list = sub.sciplex, dims=1:30,anchor.features = features)
sciplex.integrated <- IntegrateData(anchorset = sciplex.anchors,dims=1:30)
DefaultAssay(object = sciplex.integrated) <- "integrated"
sciplex.integrated <- ScaleData(sciplex.integrated, verbose = FALSE)
set.seed(10)
sciplex.integrated <- RunPCA(sciplex.integrated, npcs = 50, verbose = FALSE)
set.seed(10)
sciplex.integrated <- RunUMAP(sciplex.integrated, reduction = "pca", dims = 1:15)
ElbowPlot(sciplex.integrated,ndims=50)
set.seed(10)
sciplex.integrated <- FindNeighbors(sciplex.integrated, dims = 1:15)
sciplex.integrated <- FindClusters(sciplex.integrated, resolution = 0.2)

DimPlot(sciplex.integrated, reduction = "umap",label=T,label.size = 5, split.by = "treatment")
sciplex.integrated$cellType <- factor(sciplex.integrated$cellType, levels = c("CM","CM-4F","Fib","Prlf","pro-EpC","Unknown"))
Idents(sciplex.integrated) <- sciplex.integrated$cellType
DefaultAssay(sciplex.integrated) <- "RNA"
                                                  #MKI67           #TOP2A             #E2F1           #ANLN               #WEE1              #CDK1           #CDK4             #CCND1             #CCNB1             #DSP             #COL1A1             #COL1A2             #TTN               #MYBPC3         #MYH6              #RYR2          #MYL2         #WT1                #TBX18             #PDGFRA          #PDGFRB           #POSTN
#VlnPlot(sciplex.integrated, features = c("ENSG00000148773","ENSG00000131747","ENSG00000101412","ENSG00000011426","ENSG00000166483","ENSG00000170312","ENSG00000135446","ENSG00000110092","ENSG00000134057","ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000111245","ENSG00000184937","ENSG00000112837","ENSG00000134853","ENSG00000113721","ENSG00000133110"), ncol = 3)
FeaturePlot(sciplex.integrated, features = c("ENSG00000096696","ENSG00000108821","ENSG00000164692","ENSG00000155657","ENSG00000134571","ENSG00000197616","ENSG00000198626","ENSG00000111245"), ncol = 3, label = T, split.by = "treatment")


### DEG 

all.patient <- list[[1]]@assays$RNA@counts
all.patient <- cbind(all.patient,list[[4]]@assays$RNA@counts)
#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all.patient),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(all.patient),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(all.patient))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 229
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 1366
  
# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
table(ann$ENSEMBL==rownames(all.patient))
#TRUE 
#20633

#remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
num_cells_expressed <- rowSums(all.patient !=0)
fraction_cells_expressed <- (num_cells_expressed/ncol(all.patient))
keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
table(keep.genes)
all.patient <- all.patient[keep.genes,]
dim(all.patient)
#17187  6614
ann <- ann[keep.genes,]
  
#remove low quality nuclei
keep.nuclei <- colSums(all.patient !=0 ) >= 200
all.patient  <- all.patient[,keep.nuclei]
  
# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
xy
#FALSE  TRUE 
#16566   621
all.patient<- all.patient[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(all.patient))
#TRUE 
#16566 
sciplex.integrated <- subset(sciplex.integrated, subset=Broad_celltype!="Unknown")
index <- match(colnames(sciplex.integrated), colnames(all.patient))
all.patient <- all.patient[,index]
broadct <- factor(sciplex.integrated$Broad_celltype)
index <- which(sciplex.integrated$orig.ident =="patient_corrected")
sciplex.integrated$biorep <- "NA"
sciplex.integrated$biorep[index] <- paste("gfp",sciplex.integrated$replicate[index],sep = "")
sam <- factor(sciplex.integrated$biorep,levels=c("factorsrep1","factorsrep2","factorsrep3","factorsrep4","factorsrep5","factorsrep6","factorsrep7","factorsrep8","factorsrep9","factorsrep10","factorsrep11","factorsrep12","factorsrep13","factorsrep14","gfprep1","gfprep2","gfprep3","gfprep4","gfprep5","gfprep6","gfprep7","gfprep8","gfprep9","gfprep10","gfprep11","gfprep12","gfprep13","gfprep14"))
newgrp <- paste(broadct,sam,sep=".")
newgrp <- factor(newgrp,levels=paste(rep(levels(broadct),each=28),levels(sam),sep="."))
table(newgrp)

des <- model.matrix(~0+newgrp)
colnames(des) <- levels(newgrp)
dim(des)
#5847   112
dim(all.patient)
#16566  5847
pb <- all.patient %*% des
y.pb <- DGEList(pb)

saminfo <- matrix(unlist(strsplit(colnames(y.pb$counts),split="[.]")),ncol=2,byrow=TRUE)
bct <- factor(saminfo[,1])
indiv <- factor(saminfo[,2])
group <- rep(NA,ncol(y.pb))
group[grep("gfp",indiv)] <- "gfp"
group[grep("factors",indiv)] <- "factors"
group <- factor(group,levels=c("gfp","factors"))
table(rownames(y.pb)==ann$ENSEMBL)
#TRUE 
#16566

y.pb$genes <- ann
bct2 <- as.character(bct)
bct2[bct2 == "CM"] <- "CM"
bct2[bct2 == "pro-EpC"] <- "proEpC"
bct2[bct2 == "Fib"] <- "Fib"
bct2[bct2 == "Prlf"] <- "Prlf"


newgrp <- paste(bct2,group,sep=".")
newgrp <- factor(newgrp)

design <- model.matrix(~0+newgrp)
colnames(design) <- levels(newgrp)

v <- voom(y.pb,design,plot=TRUE,normalize.method = "cyclicloess")
fit <- lmFit(v,design)

cont.cardio <- makeContrasts(FvsGFP = Fib.factors-Fib.gfp,
                             levels=design)
fit.cardio <- contrasts.fit(fit,contrasts = cont.cardio)
fit.cardio <- eBayes(fit.cardio,robust=TRUE)

summary(decideTests(fit.cardio))
treat.cardio <- treat(fit.cardio,lfc=0)

dt.cardio<-decideTests(treat.cardio)

summary(dt.cardio)

options(digits=3)
topTreat(treat.cardio,coef=1,n=20,p.value=0.05)[,-c(1:3)]

fdr <- apply(treat.cardio$p.value, 2, function(x) p.adjust(x, method="BH"))
output <- data.frame(treat.cardio$genes,LogFC=treat.cardio$coefficients,AveExp=treat.cardio$Amean,tstat=treat.cardio$t, pvalue=treat.cardio$p.value, fdr=fdr)
library(dplyr)
hm <- c("HDAC9","FLNC","SORBS2","ANK2","WT1","FGF9","RBM24","SIX1","CDH2","RGS4","MYEF2","FLNB","TLL2","DDR2","THY1","COL1A1","COL1A2")

index <- match(hm,output$SYMBOL)
index <- na.omit(index)
hooray <- output[index,]



load("../../../../../../scRNAseq-ES/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.c2,treat.cardio$genes$ENTREZID)
reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
cardio.camera <- cameraPR(treat.cardio$t[,2],reactome.id)
cardio.camera.up <- cardio.camera[cardio.camera[,2]=="Up",]
library(dplyr)
cardio.camera.up <- cardio.camera.up %>% filter(cardio.camera.up$FDR < 0.05)

cardio.camera.dn <- cardio.camera[cardio.camera[,2]=="Down",]
library(dplyr)
cardio.camera.dn <- cardio.camera.dn %>% filter(cardio.camera.dn$FDR < 0.05)

write.csv(cardio.camera, "../../../../../../../Fib.PatientcorVSpatient.csv")

### merge heart tissue and 2D 
fetal.integrated <- readRDS(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/fetal-int.Rds")
load(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/fetalObjs.Rdata")
Idents(fetal.integrated) <- fetal.integrated$integrated_snn_res.0.3
DimPlot(fetal.integrated, reduction = "tsne",label=TRUE,label.size = 6)+NoLegend()

young.integrated <- readRDS(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/young-int.Rds")
load(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/youngObjs.Rdata")
Idents(young.integrated) <- young.integrated$integrated_snn_res.0.3
DimPlot(young.integrated, reduction = "tsne",label=TRUE,label.size = 6)+NoLegend()

dcm.integrated <- readRDS(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/dcm-int.Rds")
load(file="/group/card2/Neda/MCRI_LAB/Single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/dcmObjs.Rdata")
Idents(dcm.integrated) <- dcm.integrated$integrated_snn_res.0.3
DimPlot(dcm.integrated, reduction = "tsne",label=TRUE,label.size = 6)+NoLegend()


patient_Cor <- readRDS(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-2nd-22.12.2021/output/sciplex-patientCor.Rds")
patient_DSPmutmut <- readRDS(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-2nd-22.12.2021/output/sciplex-patient.Rds")

patient_DSPmutmut$biorep <- patient_DSPmutmut$replicate
patient_Cor$biorep <- patient_Cor$replicate

patient_DSPmutmut@meta.data <- patient_DSPmutmut@meta.data[, -which(colnames(patient_DSPmutmut@meta.data) %in% c('treatment','replicate'))]
patient_Cor@meta.data <- patient_Cor@meta.data[, -which(colnames(patient_Cor@meta.data) %in% c('treatment','replicate'))]
DefaultAssay(patient_DSPmutmut) <- "RNA"



heart <- merge(fetal.integrated, y=c(young.integrated,dcm.integrated,patient_Cor,patient_DSPmutmut), project = "heart")
heart <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-2nd-22.12.2021/output/heart-int-FYDvitro.RDS")

CM <- subset(heart, subset = Broad_celltype =="CM")

index <- which(CM$orig.ident=="patient")
CM$biorep[index] <- paste("p",CM$biorep[index],sep = "")


heart.list <- SplitObject(CM, split.by = "orig.ident")

heart.list <- SplitObject(sciplex.sub, split.by = "orig.ident")
min(sapply(heart.list, ncol))
for (i in 1:length(heart.list)) {
  heart.list[[i]] <- SCTransform(heart.list[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = heart.list)
heart.anchors <- FindIntegrationAnchors(object.list = heart.list, dims=1:30,anchor.features = features, k.filter =734 )
heart.integrated <- IntegrateData(anchorset = heart.anchors,dims=1:30)
DefaultAssay(object = heart.integrated) <- "integrated"
#cardio.integrated <- readRDS("../must-do-projects/cardio-int-FYAD-filtered.Rds")
heart.integrated <- ScaleData(heart.integrated, verbose = FALSE)
saveRDS(heart.integrated,file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/data/cm-int-FYDvitro-filtered.Rds")
heart.integrated <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/data/heart-int-YA-filtered.Rds")

set.seed(10)
heart.integrated <- RunPCA(heart.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(heart.integrated,ndims=50)
heart.integrated <- FindNeighbors(heart.integrated, dims = 1:20)
heart.integrated <- FindClusters(heart.integrated, resolution = 0.1)
set.seed(10)
heart.integrated <- RunUMAP(heart.integrated, reduction = "pca", dims = 1:20)
Idents(heart.integrated) <- heart.integrated$Broad_celltype


# Compare average expression in "High-throughput RNA sequencing of paraformaldehyde-fixed single cells"  
# link:https://www.nature.com/articles/s41467-021-25871-2    
# link:https://github.com/tay-lab/FD-seq/blob/main/technical_replicates/data_analysis.R

setwd("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-3rd-29.04.2022/code/sci-plex-master/large_screen")
library(Seurat)
library(dplyr)
library(Homo.sapiens)
library(edgeR)
source("../../../../../GENERAL_CODES/ggplotColors.R")
source("../../../../../GENERAL_CODES/normCounts.R")
source("../../../../../GENERAL_CODES/findModes.R")

# Read in data (seurat object) generated in "./Notebook3_preprocess_data.R"
data <- readRDS("../../../output/all29471-seurat.Rds")

data$replicate <- factor(data$replicate, levels = paste("rep",1:8, sep = ""))

##split cells into fresh and frozen

data.frozen <- subset(data,subset = treatment %in% c("DMSO","flash"))

`%notin%` <- Negate(`%in%`)
data.fresh  <- subset(data,subset = treatment %notin% c("DMSO","flash"))

sub <- subset(data, subset = treatment %in% c("","DMSO","flash"))
sub <- subset(sub, subset = orig.ident %in% c("patient")) 

list <- SplitObject(sub, split.by = "treatment")

all <- list()
for ( i in 1:length(list)){
  all[[i]] <- list[[i]]@assays$RNA@counts
  names(all)[[i]] <- names(list)[[i]]
}

for (i in 1:length(all)){
  #Take out mitochondrial, ribosomal and genes with no annotation
  ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all[[i]]),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
  m <- match(rownames(all[[i]]),ann$ENSEMBL)
  ann <- ann[m,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #20633
  mito <- grep("mitochondrial",ann$GENENAME)
  length(mito)
  #[1] 229
  ribo <- grep("ribosomal",ann$GENENAME)
  length(ribo)
  #[1] 201
  missingEZID <- which(is.na(ann$ENTREZID))
  length(missingEZID)
  #[1] 1366
  
  # Filter out non-informative genes
  chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
  length(chuck)
  all[[i]] <- all[[i]][-chuck,]
  ann <- ann[-chuck,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  #  num_cells_expressed <- rowSums(all[[i]] !=0)
  #  fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
  #  keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
  #  table(keep.genes)
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  #num_cells_expressed <- rowSums(all[[i]] !=0)
  numzero.genes <- rowSums(all[[i]] ==0)
  #fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
  #keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
  table(numzero.genes > (ncol(all[[i]])-20))
  keep.genes <- numzero.genes < (ncol(all[[i]])-20)
  all[[i]] <- all[[i]][keep.genes,]
  dim(all[[i]])
  ann <- ann[keep.genes,]
  
  #remove low quality nuclei
  keep.nuclei <- colSums(all[[i]] !=0 ) >= 200
  all[[i]]  <- all[[i]][,keep.nuclei]
  
  # Take out sex chromosome genes
  xy <- ann$TXCHROM %in% c("chrX","chrY")
  table(xy)
  xy
  #FALSE  TRUE 
  #13571   516
  all[[i]] <- all[[i]][!xy,]
  ann <- ann[!xy,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #13571
}  

orig.ident <- list()
for ( i in 1:length(all)){
  orig.ident[[i]] <- rep(levels(factor(list[[i]]@meta.data$orig.ident)),ncol(all[[i]]))
  names(orig.ident)[[i]] <- names(all)[[i]]
}

replicate <- list()
for ( i in 1:length(all)){
  index <- match(colnames(all[[i]]),colnames(list[[i]]))
  replicate[[i]] <- list[[i]]$replicate[index]
  names(replicate)[[i]] <- names(all)[[i]]
}

treatment <- list()
for ( i in 1:length(all)){
  treatment[[i]] <- rep(levels(factor(list[[i]]@meta.data$treatment)),ncol(all[[i]]))
  names(treatment)[[i]] <- names(all)[[i]]
}

sciplex <- list()
for (i in 1:length(all)){
  sciplex[[i]] <- CreateSeuratObject(counts = all[[i]], project = names(all)[[i]])
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = orig.ident[[i]], col.name="orig.ident")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = replicate[[i]], col.name="replicate")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = treatment[[i]], col.name="treatment")
  names(sciplex)[[i]] <- names(all)[[i]]
}

sample <- merge(sciplex[[1]], y=sciplex[[3]], add.cell.ids = c('DMSO',''))

sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample, selection.method = 'vst')
sample <- ScaleData(sample)
Idents(sample) <- sample$treatment
avg.sample <- data.frame(log1p(AverageExpression(sample)$RNA))
colnames(avg.sample)[2] <- "Fresh"
library(ggplot2)
my.theme <- theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
                 axis.text = element_text(size=12, color='black'), legend.position = "none")
par(mfrow=c(2,1))
par(mar=c(6,4,10,10))
pdf(file="../../../../../../../seurat_DMSO_vs_Fresh_expression.pdf",width = 5,height = 5)
ggplot(avg.sample, aes(DMSO,Fresh)) +
  geom_point(size=1, colour = "gray60") +
  geom_abline(slope = 1, intercept = 0, colour = "red", size=1) +
  xlab("Log-normalized\nexpression level (DMSO)") + ylab("Log-normalized\nexpression level (Fresh)") +
  theme_classic() + my.theme
dev.off()
#ggsave("../../../../../../../seurat_rep1_vs_rep2_expression.png", dpi=600, width=2.7, height=2.7, units="in")
cor(avg.sample, method="pearson")
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(avg.sample),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(avg.sample),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(avg.sample))
index <- match(rownames(avg.sample),ann$ENSEMBL)
avg.sample <- cbind(avg.sample,ann$SYMBOL[index])
write.csv(avg.sample,file="../../../../../../../avg.sample.csv")

#CVs vs mean
set.seed(1)
index <- sample(colnames(sciplex[[3]]), 278) 
Fresh.sub <- sciplex[[3]][,index]

DMSO <- sciplex[[1]]
DMSO <- NormalizeData(DMSO, normalization.method = "RC")
Fresh.sub <- NormalizeData(Fresh.sub,normalization.method = "RC")

avg.exp.DMSO <- apply(DMSO@assays$RNA@data,1,mean)
var.DMSO <- apply(DMSO@assays$RNA@data,1,var)
CV2.DMSO <- var.DMSO/avg.exp.DMSO^2

avg.exp.Fresh <- apply(Fresh.sub@assays$RNA@data,1,mean)
var.Fresh <- apply(Fresh.sub@assays$RNA@data,1,var)
CV2.Fresh <- var.Fresh/avg.exp.Fresh^2
#DMSO <- scran::modelGeneCV2(sciplex[[1]]@assays$RNA@counts, size.factors = 1000000)
#Fresh.sub <- scran::modelGeneCV2(Fresh.sub@assays$RNA@counts)

df_1 <- data.frame(cbind(avg.exp.DMSO,CV2.DMSO))
df_2 <- data.frame(cbind(avg.exp.Fresh,CV2.Fresh))

#find HKG and non-HKG  link:https://www.tau.ac.il/~elieis/HKG/   paper:https://www.sciencedirect.com/science/article/pii/S1097276518308803#sec4
HKG <- read.csv("../../../../../../../HK_genes",sep = "", header = F)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(df_1),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(df_1),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(df_1))
index <- match(rownames(df_1),ann$ENSEMBL)
df_1 <- cbind(df_1,ann$SYMBOL[index])

ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(df_2),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(df_2),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(df_2))
index <- match(rownames(df_2),ann$ENSEMBL)
df_2 <- cbind(df_2,ann$SYMBOL[index])

# Visualizing the mean-CV2 fit, coloring treatment
index <- order(df_1$avg.exp.DMSO,decreasing = T)[1:2000]
df_1o <- tibble(
  mean = df_1$avg.exp.DMSO[index], 
  CV2 = df_1$CV2.DMSO[index], 
  #trend = metadata(DMSO)$trend(mean),
  SYMBOL = df_1[index,3],
  treatment = "DMSO"
)

index <- match(HKG$V1,df_1o$SYMBOL)
index <- na.omit(index)
df_1o$HKgene <- 0
df_1o$HKgene[index] <- 1

index <- order(df_2$avg.exp.Fresh,decreasing = T)[1:2000]
df_2o <- tibble(
  mean = df_2$avg.exp.Fresh[index], 
  CV2 = df_2$CV2.Fresh[index], 
  SYMBOL = df_2[index,3],
  #trend = metadata(Fresh.sub)$trend(mean),
  treatment = "Fresh"
)

index <- match(HKG$V1,df_2o$SYMBOL)
index <- na.omit(index)
df_2o$HKgene <- 0
df_2o$HKgene[index] <- 1

df <- rbind(df_1o,df_2o)
df$treatment <- as.factor(df$treatment)

df$HKgene <- as.factor(df$HKgene)
df_2o$HKgene <- as.factor(df_2o$HKgene)
p <- ggplot(df) + 
  geom_point(aes(x = mean, y = CV2, col = HKgene), alpha = 0.4) + 
  #geom_line(aes(x = mean, y = trend), col = 'darkred') +
  scale_x_log10() +
  theme_minimal() +
  labs(x = 'Average of UMIs per million', y = 'CV2')+theme_classic() + my.theme
ggsave('../../../../../../../cv2~mean.pdf')


genes <- read.csv("../../../../../../../ENS",sep = "", header = F)
numZero <- data.frame(rowSums(all[[1]] !=0))
index <- match(genes$V1,rownames(numZero))
index <- na.omit(index)
DMSO.genes <- data.frame(c$rowSums.a....0.[index])

numZero.f <- data.frame(rowSums(all[[3]] !=0))
index <- match(genes$V1,rownames(numZero.f))
index <- na.omit(index)
fresh.genes <- data.frame(numZero.f$rowSums.all..3......0.[index])

ann <- AnnotationDbi:::select(Homo.sapiens,keys=g$`genes$V1`,columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(g$`genes$V1`,ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==g$`genes$V1`)
index <- match(g$`genes$V1`,ann$ENSEMBL)
g <- cbind(g,ann$SYMBOL[index])



#Check batch effect between sciplex 2&3 for patient line


# Read in data (seurat object) generated in "./Notebook3_preprocess_data.R"
#data is from sciplex3
data <- readRDS("../../../output/all29471-seurat.Rds")

data$replicate <- factor(data$replicate, levels = paste("rep",1:8, sep = ""))

##split cells into fresh and frozen

data.frozen <- subset(data,subset = treatment %in% c("DMSO","flash"))

`%notin%` <- Negate(`%in%`)
data.fresh  <- subset(data,subset = treatment %notin% c("DMSO","flash"))
### GO TO LINE 1031


list <- list()
for (cell_line in levels(factor(data.fresh$orig.ident))){
  for (treatment in levels(factor(data.fresh$treatment))){
    for (replicate in levels(factor(data.fresh[,which(data.fresh$orig.ident==cell_line & data.fresh$treatment==treatment)]$replicate))){
      if (treatment==""){
        list[[paste(cell_line,replicate,sep = "_")]] <- data.fresh[,which(data.fresh$orig.ident==cell_line & data.fresh$treatment == treatment & data.fresh$replicate == replicate)]
        list[[paste(cell_line,replicate,sep = "_")]]$replicate <- rep(replicate,ncol(list[[paste(cell_line,replicate,sep = "_")]]))
      }
      else{list[[paste(cell_line,treatment,replicate,sep = "_")]] <- data.fresh[,which(data.fresh$orig.ident==cell_line & data.fresh$treatment == treatment & data.fresh$replicate == replicate)] 
      list[[paste(cell_line,treatment,replicate,sep = "_")]]$replicate <-  rep(replicate,ncol(list[[paste(cell_line,treatment,replicate,sep = "_")]]))
      }
    }
  }
}

#data is from sciplex2
data.prev <- readRDS("../../../../sci-RNA-seq3-2nd-22.12.2021/output/all23294-seurat.Rds")

data.prev$replicate <- factor(data.prev$replicate, levels = paste("rep",1:16, sep = ""))


list.prev <- list()
for (cell_line in levels(factor(data.prev$orig.ident))){
  for (treatment in levels(factor(data.prev$treatment))){
    for (replicate in levels(factor(data.prev[,which(data.prev$orig.ident==cell_line & data.prev$treatment==treatment)]$replicate))){
      if (treatment==""){
        list.prev[[paste(cell_line,replicate,sep = "_")]] <- data.prev[,which(data.prev$orig.ident==cell_line & data.prev$treatment == treatment & data.prev$replicate == replicate)]
        list.prev[[paste(cell_line,replicate,sep = "_")]]$replicate <- rep(replicate,ncol(list.prev[[paste(cell_line,replicate,sep = "_")]]))
      }
      else{list.prev[[paste(cell_line,treatment,replicate,sep = "_")]] <- data.prev[,which(data.prev$orig.ident==cell_line & data.prev$treatment == treatment & data.prev$replicate == replicate)] 
      list.prev[[paste(cell_line,treatment,replicate,sep = "_")]]$replicate <-  rep(replicate,ncol(list.prev[[paste(cell_line,treatment,replicate,sep = "_")]]))
      }
    }
  }
}

list.final <- list()
list.final[1:8] <- list[1:8]
list.final[9:22] <- list.prev[1:14]
list.final[23:30] <- list[25:32]
list.final[31:44] <- list.prev[43:56]

names(list.final) <- c("patient_sciplex3_rep1","patient_sciplex3_rep2","patient_sciplex3_rep3","patient_sciplex3_rep4","patient_sciplex3_rep5","patient_sciplex3_rep6","patient_sciplex3_rep7","patient_sciplex3_rep8",
                        "patient_sciplex2_rep1","patient_sciplex2_rep2","patient_sciplex2_rep3","patient_sciplex2_rep4","patient_sciplex2_rep5","patient_sciplex2_rep6","patient_sciplex2_rep7","patient_sciplex2_rep8","patient_sciplex2_rep9",
                       "patient_sciplex2_rep10","patient_sciplex2_rep11","patient_sciplex2_rep12","patient_sciplex2_rep13","patient_sciplex2_rep14",
                       "patient.cor_sciplex3_rep1","patient.cor_sciplex3_rep2","patient.cor_sciplex3_rep3","patient.cor_sciplex3_rep4","patient.cor_sciplex3_rep5","patient.cor_sciplex3_rep6","patient.cort_sciplex3_rep7","patient.cor_sciplex3_rep8",
                       "patient.cor_sciplex2_rep1","patient.cor_sciplex2_rep2","patient.cor_sciplex2_rep3","patient.cor_sciplex2_rep4","patient.cor_sciplex2_rep5","patient.cor_sciplex2_rep6","patient.cor_sciplex2_rep7","patient.cor_sciplex2_rep8","patient.cor_sciplex2_rep9",
                       "patient.cor_sciplex2_rep10","patient.cor_sciplex2_rep11","patient.cor_sciplex2_rep12","patient.cor_sciplex2_rep13","patient.cor_sciplex2_rep14")

# MDS plot of all samples
#To get a high-level idea of the overall sources of variability in the dataset, I have summed the counts over all cells within a sample to obtain a “pseudobulk” sample and made MDS plots using functions in edgeR.
for (i in 1:length(list.final)){
  #remove low quality nuclei
  keep.nuclei <- colSums(list.final[[i]]@assays$RNA@counts !=0 ) >= 200
  list.final[[i]]  <- list.final[[i]][,keep.nuclei]
}  

pseudobulk <- matrix(NA,ncol=length(list.final),nrow=nrow(list.final[[1]]))
colnames(pseudobulk) <- names(list.final)
rownames(pseudobulk) <- rownames(list.final[[1]])
for (i in 1:length(list.final))
  pseudobulk[,i] <- rowSums(list.final[[i]])

#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(pseudobulk),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(pseudobulk),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(pseudobulk))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 231
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 775

# Filter out non-informative genes
#ENSG00000251562 stands for MALAT1
chuck <- unique(c(mito,ribo,missingEZID,which(rownames(pseudobulk)=="ENSG00000251562")))
length(chuck)
#[1] 1127

pseudobulk <- pseudobulk[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(pseudobulk))
#TRUE 
#19506

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
#xy
#FALSE  TRUE 
#18698   808
pseudobulk <- pseudobulk[!xy,]
ann <- ann[!xy,]

#Filter lowly expressed genes
y <- DGEList(pseudobulk)
keep <- rowSums(cpm.DGEList(y)>=35)>=8  #median(y$samples$lib.size) is ~159,988. A CPM of 62 is used as it corresponds to a count of 10 for this median library size (~159,988) in this data set. 8 = minimum number of replicates per condition.
table(keep)
#keep
#FALSE  TRUE 
#11758  6940 
y <- y[keep,]
y <- calcNormFactors(y)
dim(y)
#[1]6940   44

y$samples$group <- c(rep("patient_sciplex3",8),rep("patient_sciplex2",14),rep("patient.cor_sciplex3",8),rep("patient.cor_sciplex2",14))
y$samples$group <- factor(y$samples$group, levels = c("patient_sciplex3","patient_sciplex2","patient.cor_sciplex3","patient.cor_sciplex2"))
#,"tan4","blue","chocolate","violetred4","darkgoldenrod2","tomato","tan4","blue","chocolate"
y$samples$colour = rep(c("darkorchid","darkgoldenrod1","blue","chocolate"),4)[factor(y$samples$group)]
pdf(file="/group/card2/Neda/plotMDS-WT.mutationDMSOflash.pdf",width = 10,height = 5)
plotMDS(y, top = 500,dim=c(1,2), pch=rep(c(21,21,21,21),4)[factor(y$samples$group)], bg=y$samples$colour, cex = 2, gene.selection = "common")   #character, "pairwise" to choose the top genes separately for each pairwise comparison between the samples or "common" to select the same genes for all comparisons.
legend("bottomleft", title="Condition", # << THIS IS THE HACKISH PART
       legend=c("patient_sciplex3","patient_sciplex2","patient.cor_sciplex3","patient.cor_sciplex2"), ncol=1, box.lwd = 1)
abline(v=0,h=0,col="gray", lwd=1, lty=2)

#prcomp(y)
dev.off()

### generate a UMAP 12.08.2022

# Read in data (seurat object) generated in "./Notebook3_preprocess_data.R"
#data is from sciplex3
data <- readRDS("../../../output/all29471-seurat.Rds")

data$replicate <- factor(data$replicate, levels = paste("rep",1:8, sep = ""))


`%notin%` <- Negate(`%in%`)
data.fresh  <- subset(data,subset = treatment %notin% c("DMSO","flash"))
sub <- subset(data.fresh, subset = treatment %in% c(""))
sub <- subset(sub, subset = orig.ident %in% c("patient", "patient.corrected"))  ##"WT.mutation"
list <- SplitObject(sub, split.by = "orig.ident")


#data is from sciplex2
data.prev <- readRDS("../../../../sci-RNA-seq3-2nd-22.12.2021/output/all23294-seurat.Rds")
data.prev$replicate <- factor(data.prev$replicate, levels = paste("rep",1:16, sep = ""))
sub <- subset(data.prev, subset = treatment %in% c(""))
sub <- subset(sub, subset = orig.ident %in% c("patient", "patient_corrected"))  ##"WT.mutation"
list.prev <- SplitObject(sub, split.by = "orig.ident")


list.final <- list()
list.final[1] <- list[1]
list.final[2] <- list.prev[1]
list.final[3] <- list[2]
list.final[4] <- list.prev[2]

names(list.final) <- c("patient_sciplex3",
                       "patient_sciplex2",
                       "patient.cor_sciplex3",
                       "patient.cor_sciplex2")


all <- list()
for ( i in 1:length(list.final)){
  all[[i]] <- list.final[[i]]@assays$RNA@counts
  names(all)[[i]] <- names(list.final)[[i]]
}


for (i in 1:length(all)){
  #Take out mitochondrial, ribosomal and genes with no annotation
  ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all[[i]]),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
  m <- match(rownames(all[[i]]),ann$ENSEMBL)
  ann <- ann[m,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #20633
  mito <- grep("mitochondrial",ann$GENENAME)
  length(mito)
  #[1] 229
  ribo <- grep("ribosomal",ann$GENENAME)
  length(ribo)
  #[1] 201
  missingEZID <- which(is.na(ann$ENTREZID))
  length(missingEZID)
  #[1] 1366
  
  # Filter out non-informative genes
  chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
  length(chuck)
  all[[i]] <- all[[i]][-chuck,]
  ann <- ann[-chuck,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  #  num_cells_expressed <- rowSums(all[[i]] !=0)
  #  fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
  #  keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
  #  table(keep.genes)
  
  #remove low-expressed genes according to this tutorial https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  #num_cells_expressed <- rowSums(all[[i]] !=0)
  numzero.genes <- rowSums(all[[i]] ==0)
  #fraction_cells_expressed <- (num_cells_expressed/ncol(all[[i]]))
  #keep.genes <- fraction_cells_expressed >= (0.01/100) #0.01% :) 
  table(numzero.genes > (ncol(all[[i]])-20))
  keep.genes <- numzero.genes < (ncol(all[[i]])-20)
  all[[i]] <- all[[i]][keep.genes,]
  dim(all[[i]])
  ann <- ann[keep.genes,]
  
  #remove low quality nuclei
  keep.nuclei <- colSums(all[[i]] !=0 ) >= 200
  all[[i]]  <- all[[i]][,keep.nuclei]
  
  # Take out sex chromosome genes
  xy <- ann$TXCHROM %in% c("chrX","chrY")
  table(xy)
  xy
  #FALSE  TRUE 
  #13571   516
  all[[i]] <- all[[i]][!xy,]
  ann <- ann[!xy,]
  table(ann$ENSEMBL==rownames(all[[i]]))
  #TRUE 
  #13571
}  


orig.ident <- list()
for ( i in 1:length(all)){
  orig.ident[[i]] <- rep(levels(factor(list.final[[i]]@meta.data$orig.ident)),ncol(all[[i]]))
  names(orig.ident)[[i]] <- names(all)[[i]]
}

replicate <- list()
for ( i in 1:length(all)){
  index <- match(colnames(all[[i]]),colnames(list.final[[i]]))
  replicate[[i]] <- list.final[[i]]$replicate[index]
  names(replicate)[[i]] <- names(all)[[i]]
  replicate[[i]] <- as.character(replicate[[i]])
}

treatment <- list()
for ( i in 1:length(all)){
  treatment[[i]] <- rep(levels(factor(list.final[[i]]@meta.data$treatment)),ncol(all[[i]]))
  names(treatment)[[i]] <- names(all)[[i]]
}

sciplex <- list()
for (i in 1:length(all)){
  sciplex[[i]] <- CreateSeuratObject(counts = all[[i]], project = names(all)[[i]])
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = orig.ident[[i]], col.name="orig.ident")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = replicate[[i]], col.name="replicate")
  sciplex[[i]] <- AddMetaData(object=sciplex[[i]], metadata = treatment[[i]], col.name="treatment")
  names(sciplex)[[i]] <- names(all)[[i]]
}

for (i in 1:length(sciplex)){
  sciplex[[i]] <- NormalizeData(sciplex[[i]],scale.factor = 10000)
  sciplex[[i]] <- FindVariableFeatures(sciplex[[i]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(sciplex[[i]])
  sciplex[[i]] <- ScaleData(sciplex[[i]], features = all.genes)
  sciplex[[i]] <- RunPCA(sciplex[[i]], features = VariableFeatures(object = sciplex[[i]]))
  sciplex[[i]] <- RunUMAP(sciplex[[i]], reduction = "pca", dims = 1:15)
  sciplex[[i]] <- FindNeighbors(sciplex[[i]], dims = 1:15)
  sciplex[[i]] <- FindClusters(sciplex[[i]], resolution = 0.1)
}

sciplex[[1]] <- FindClusters(sciplex[[1]], resolution = 0.1)
sciplex[[1]] <- RunUMAP(sciplex[[1]], reduction = "pca", dims = 1:15)
DimPlot(sciplex[[1]], reduction = "umap",label=T,label.size = 5)

sciplex[[1]]$Broad_celltype <- "NA"
index <- which(sciplex[[1]]$RNA_snn_res.0.1=="1")
sciplex[[2]]$Broad_celltype[index] <- "Fib"
saveRDS(sciplex,file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/sci-RNA-seq3-3rd-29.04.2022/output/sciplex-annotatedCellType.Rds")
Idents(sciplex[[2]]) <- sciplex[[2]]$Broad_celltype


sciplex.markers <- FindAllMarkers(sciplex[[1]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(sciplex[[1]]),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(sciplex[[1]]),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(sciplex[[1]]))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])


sciplex[[1]]$sciplex <- NA
sciplex[[1]]$sciplex <- "patient_sciplex3"
sciplex[[2]]$sciplex <- NA
sciplex[[2]]$sciplex <- "patient_sciplex2"
sciplex.sub <- merge(x = sciplex[[1]],y = sciplex[[2]])
sciplex.sub <- subset(sciplex.sub, subset=Broad_celltype != "NA")

heart.list <- SplitObject(sciplex.sub, split.by = "sciplex")
min(sapply(heart.list, ncol))
for (i in 1:length(heart.list)) {
  heart.list[[i]] <- SCTransform(heart.list[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = heart.list)
heart.anchors <- FindIntegrationAnchors(object.list = heart.list, dims=1:30,anchor.features = features, k.filter = 734)
heart.integrated <- IntegrateData(anchorset = heart.anchors,dims=1:30)
DefaultAssay(object = heart.integrated) <- "integrated"
heart.integrated <- ScaleData(heart.integrated, verbose = FALSE)
saveRDS(heart.integrated,file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/data/cm-int-FYDvitro-filtered.Rds")
heart.integrated <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/data/heart-int-YA-filtered.Rds")

set.seed(10)
heart.integrated <- RunPCA(heart.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(heart.integrated,ndims=50)
heart.integrated <- FindNeighbors(heart.integrated, dims = 1:20)
heart.integrated <- FindClusters(heart.integrated, resolution = 0.1)
set.seed(10)
heart.integrated <- RunUMAP(heart.integrated, reduction = "pca", dims = 1:20)
Idents(heart.integrated) <- heart.integrated$Broad_celltype
DimPlot(heart.integrated, reduction = "umap",label=T,label.size = 5, split.by = "sciplex")
DimPlot(heart.integrated, reduction = "umap",label.size = 5)
Idents(heart.integrated) <- heart.integrated$sciplex
DefaultAssay(heart.integrated) <- "integrated"

DefaultAssay(heart.integrated) <- "RNA"
sciplex.markers <- FindAllMarkers(heart.integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
#sciplex.markers <- sciplex.markers %>% filter(sciplex.markers$p_val_adj < 0.05)
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(heart.integrated),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(heart.integrated),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(heart.integrated))
index <- match(sciplex.markers$gene,ann$ENSEMBL)
sciplex.markers <- cbind(sciplex.markers,ann$SYMBOL[index])

DefaultAssay(heart.integrated) <- "SCT"       #RYR2          #TNNT2                  #TTN         #COL1A2         #FN1                #MYBPC3
VlnPlot(heart.integrated, features = c("ENSG00000198626","ENSG00000118194","ENSG00000155657","ENSG00000164692","ENSG00000115414","ENSG00000134571","ENSG00000168542"),ncol = 3, split.by = "sciplex")
VlnPlot(heart.integrated, features = c("ENSG00000198626"), ncol = 3)
VlnPlot(heart.integrated, features = c("ENSG00000198626"), ncol = 3, split.by = "sciplex")        
        

#04.04.2023
#Mapping and annotating query datasets
#Tutorial link: https://satijalab.org/seurat/articles/integration_mapping.html

library(Seurat)

heart <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/EvangelynSim-snRNAseq/data/heart-int-FND-filtered.Rds")
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
saveRDS(FY.integrated, file="/group/card2/Neda/FY.integrated.Rds")
FY.integrated <- readRDS(file="/group/card2/Neda/MCRI_LAB/must-do-projects/JamesHudsonLab/COVID19_CytokineStorm_hCOs/output/Endo-Data-Cleaning/FY.integrated.Rds")

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


ctrl <- sciplex[[1]]
DefaultAssay(ctrl) <- "RNA"

FY.anchors <- FindTransferAnchors(reference = FY.integrated, query = ctrl,
                                  dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = FY.anchors, refdata = FY.integrated$bc,
                            dims = 1:30)
ctrl <- AddMetaData(ctrl, metadata = predictions)

#ctrl$prediction.match <- ctrl$predicted.id == ctrl$Broad_celltype
#table(ctrl$prediction.match)
table(ctrl$predicted.id)

FY.integrated <- RunUMAP(FY.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
ctrl <- MapQuery(anchorset = FY.anchors, reference = FY.integrated, query = ctrl,
                 refdata = list(celltype = "Broad_celltype"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(FY.integrated, reduction = "umap", group.by = "Broad_celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(ctrl, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)+ ggtitle("Query transferred labels")
p1 + p2


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

saveRDS(ctrl,file= "/group/card2/Neda/ctrl-filtered.Rds")
ctrl <- readRDS("/group/card2/Neda/MCRI_LAB/must-do-projects/JamesHudsonLab/COVID19_CytokineStorm_hCOs/output/Endo-Data-Cleaning/ctrl.Rds")
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


#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(all),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(all))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 229
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 1366

# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
all <- all[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(all))
numzero.genes <- rowSums(all@assays$RNA@counts ==0)
table(numzero.genes > (ncol(all)-20))
keep.genes <- numzero.genes < (ncol(all)-20)
all <- all[keep.genes,]
dim(all)
ann <- ann[keep.genes,]

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
xy

all <- all[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(all))


broadct <- "CM"
index <- which(cm.subset$orig.ident =="MCHTB11.5" | cm.subset$orig.ident =="MCHTB11.9"| cm.subset$orig.ident =="MCHTB11.4")
index <- na.omit(index)
cm.subset$biorep <- "NA"
cm.subset$biorep[index] <- paste("pt",cm.subset$replicate[index],sep = "")
index <- which(cm.subset$orig.ident =="MCHTB11.9 S2")
index <- na.omit(index)
cm.subset$biorep[index] <- paste("cor",cm.subset$replicate[index],sep = "")
sam <- factor(cm.subset$biorep,levels=c("correp1","correp2","correp3","correp4","correp5","correp6","correp7","correp8","correp9","correp10","correp11","correp12","correp13","correp14","correp15","correp16","ptrep1","ptrep2","ptrep3","ptrep4","ptrep5","ptrep6","ptrep7","ptrep8","ptrep9","ptrep10","ptrep11","ptrep12","ptrep13","ptrep14","ptrep15","ptrep16"))
#sam <- factor(cm.subset$biorep,levels=c("correp1","correp2","correp3","correp4","correp5","correp6","correp7","correp8","correp9","correp10","correp11","correp12","correp13","correp14","correp15","correp16","ptrep1","ptrep2","ptrep3","ptrep4","ptrep5","ptrep6","ptrep7","ptrep8"))
newgrp <- paste(broadct,sam,sep=".")
newgrp <- factor(newgrp,levels=paste("CM",levels(sam),sep="."))
table(newgrp)

des <- model.matrix(~0+newgrp)
colnames(des) <- levels(newgrp)
dim(des)
#33514    32
dim(all)
#12696 33514
pb <- all@assays$RNA@counts %*% des
y.pb <- DGEList(pb)

saminfo <- matrix(unlist(strsplit(colnames(y.pb$counts),split="[.]")),ncol=2,byrow=TRUE)
bct <- factor(saminfo[,1])
indiv <- factor(saminfo[,2])
group <- rep(NA,ncol(y.pb))
group[grep("cor",indiv)] <- "cor"
group[grep("pt",indiv)] <- "pt"
group <- factor(group,levels=c("cor","pt"))
table(rownames(y.pb)==ann$ENSEMBL)
#TRUE 
#12696 

y.pb$genes <- ann
bct2 <- as.character(bct)
bct2[bct2 == "CM"] <- "CM"
#bct2[bct2 == "Fib"] <- "Fib"



newgrp <- paste(bct2,group,sep=".")
newgrp <- factor(newgrp)

design <- model.matrix(~0+newgrp)
colnames(design) <- levels(newgrp)

par(mfrow=c(1,1))
y.pb <- calcNormFactors(y.pb)

######
all <- readRDS("../../../../sci-RNA-seq3-5th-DSP-TPM/output/first-batch/all52090-seurat.Rds")
all$replicate <- factor(all$replicate, levels = paste("rep",1:16, sep = ""))
index <- match(colnames(data.int), colnames(all))
index <- na.omit(index)
all <- all[,index]

#Take out mitochondrial, ribosomal and genes with no annotation
ann <- AnnotationDbi:::select(Homo.sapiens,keys=rownames(all),columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","TXCHROM"),keytype = "ENSEMBL")
m <- match(rownames(all),ann$ENSEMBL)
ann <- ann[m,]
table(ann$ENSEMBL==rownames(all))
#TRUE 
#20633
mito <- grep("mitochondrial",ann$GENENAME)
length(mito)
#[1] 229
ribo <- grep("ribosomal",ann$GENENAME)
length(ribo)
#[1] 201
missingEZID <- which(is.na(ann$ENTREZID))
length(missingEZID)
#[1] 1366

# Filter out non-informative genes
chuck <- unique(c(mito,ribo,missingEZID,which(ann$ENSEMBL=="ENSG00000251562")))
length(chuck)
all <- all[-chuck,]
ann <- ann[-chuck,]
table(ann$ENSEMBL==rownames(all))
numzero.genes <- rowSums(all@assays$RNA@counts ==0)
table(numzero.genes > (ncol(all)-20))
keep.genes <- numzero.genes < (ncol(all)-20)
all <- all[keep.genes,]
dim(all)
ann <- ann[keep.genes,]

# Take out sex chromosome genes
xy <- ann$TXCHROM %in% c("chrX","chrY")
table(xy)
xy

all <- all[!xy,]
ann <- ann[!xy,]
table(ann$ENSEMBL==rownames(all))

#hm <- read.delim("/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/Fetal-Gene-Program-snRNAseq/data/cellTypeMarkers.txt",stringsAsFactors = FALSE, sep="\t", header = T)
hm <- c("ANLN","MKI67","CENPF","ACTN2","MYH7","TNNT2","COL1A2","COL1A1","VCAN")
hgene <- toupper(hm)
hgene <- unique(hgene)
y <- DGEList(all@assays$RNA@counts)
log_counts <- normCounts(y,log=TRUE,prior.count=0.5)


m <- match(ann$ENSEMBL[ann$SYMBOL %in% hgene],rownames(log_counts))
m <- m[!is.na(m)]
broadct <- factor(data.int$Broad_celltype)
sumexpr <- matrix(NA,nrow=dim(log_counts)[1],ncol=length(levels(broadct)))
rownames(sumexpr) <- rownames(log_counts)
colnames(sumexpr) <- levels(broadct)

for(i in 1:nrow(sumexpr)){
  sumexpr[i,] <- tapply(log_counts[i,],broadct,mean)
}
library(NMF)
#pdf(file="/group/card2/Neda/MCRI_LAB/must-do-projects/EnzoPorrelloLab/dilated-cardiomyopathy/Fetal_Gene_Program_snRNAseq/output/Figures/cellTypeMarkers.pdf",width=30,height=25)
aheatmap(sumexpr[m,],Rowv = NA,Colv = NA, labRow = ann$SYMBOL[m],
         fontsize=5,color="-RdYlBu",cexRow =1, cexCol = 1,
         scale="row")
logcounts.all <- normCounts(all.keep,log=TRUE,prior.count=0.5)

all.bct <- factor(data.int$Broad_celltype,
                  levels = c("CM","CM(Prlf)","Fib"))
sample <- data.int$density
design <- model.matrix(~0+all.bct+sample)
colnames(design)[1:(length(levels(all.bct)))] <- levels(all.bct)

mycont <- matrix(0,ncol=length(levels(all.bct)),nrow=length(levels(all.bct)))
colnames(mycont)<-levels(all.bct)
diag(mycont)<-1
mycont[upper.tri(mycont)]<- -1/(length(levels(all.bct))-1)
mycont[lower.tri(mycont)]<- -1/(length(levels(all.bct))-1)

# Fill out remaining rows with 0s
zero.rows <- matrix(0,ncol=length(levels(all.bct)),nrow=(ncol(design)-length(levels(all.bct))))
test <- rbind(mycont,zero.rows)

fit <- lmFit(log_counts,design)
fit.cont <- contrasts.fit(fit,contrasts=test)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

fit.cont$genes <- ann$ENTREZID

treat.all <- treat(fit.cont,lfc=0.5)

#DefaultAssay(data.) <- "RNA"

sig.genes <- gene.label <- vector("list", ncol(treat.all))
for(i in 1:length(sig.genes)){
  top <- topTreat(treat.all,coef=i,n=Inf,sort.by="t")
  sig.genes[[i]] <- rownames(top)[top$logFC>0][1:10]
  gene.label[[i]] <- paste(rownames(top)[top$logFC>0][1:10],colnames(treat.all)[i],sep="-")
} 

csig <- unlist(sig.genes)
genes <- unlist(gene.label)

missing <- is.na(match(csig,rownames(data.int)))

csig2 <- csig[!missing]

gene.cols <- rep(c(ggplotColors(8),"grey"),each=10)
gene.cols <- gene.cols[!missing]

d <- duplicated(csig2)
csig2 <- csig2[!d]
gene.cols <- gene.cols[!d]
load("/group/card2/Neda/MCRI_LAB/single_cell_nuclei_rnaseq/Porello-heart-snRNAseq/output/RDataObjects/human_c2_v5p2.rdata")
c2.id <- ids2indices(Hs.c2,treat.all$genes)
reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
cm.camera <- cameraPR(treat.all$t[,1],reactome.id)
cm.camera.up <- cm.camera[cm.camera[,2]=="Up",]

cardioprlf.camera <- cameraPR(treat.all$t[,2],reactome.id)
cardioprlf.camera.up <- cardioprlf.camera[cardioprlf.camera[,2]=="Up",]

fib.camera <- cameraPR(treat.all$t[,3],reactome.id)
fib.camera.up <- fib.camera[fib.camera[,2]=="Up",]


nsets <- 5
all.cam <- rbind(cm.camera.up[1:nsets,], cardioprlf.camera.up[1:nsets,],
                 fib.camera.up[1:nsets,])

scores <- -log10(all.cam$PValue)
names(scores) <- rownames(all.cam)
names(scores) <- gsub("REACTOME_","",names(scores))

par(mfrow=c(1,1))
par(mar=c(5,20,3,2))
barplot(scores[length(scores):1],horiz = T,las=2,col=rev(rep(c(ggplotColors(8),"grey"),each=nsets)),cex.names=0.9,
        cex.axis = 1.5,xlab="-log10(PValue)",cex.lab=1.5)
abline(v= -log10(0.05),lty=2)


cm.invitro <- read.csv("/group/card2/Neda/CM-both-updnns.csv")
cm.invivo <- read.csv("/group/card2/Neda/invivo-cm-DvsN.csv")
table(cm.invivo$fdr.DNcm < 0.05 & cm.invivo$log2FC.Fetal.ND..in.CM < 0  )
table(cm.invitro$fdr < 0.05 & cm.invitro$LogFC<0 )
invivo.dn <- cm.invivo %>% filter(cm.invivo$fdr.DNcm < 0.05 & cm.invivo$LogFC.DNcm < 0 )
invivo.up <- cm.invivo %>% filter(cm.invivo$fdr.DNcm < 0.05 & cm.invivo$LogFC.DNcm > 0 )
invitro.up <- cm.invitro %>% filter(cm.invitro$fdr < 0.05 & cm.invitro$LogFC > 0 )
invitro.dn <- cm.invitro %>% filter(cm.invitro$fdr < 0.05 & cm.invitro$LogFC < 0 )
share.up <- intersect(invivo.up$ENSEMBL, invitro.up$ENSEMBL)
share.dn <- intersect(invivo.dn$ENSEMBL, invitro.dn$ENSEMBL)
up <- cm.invitro[cm.invitro$ENSEMBL %in% share.up,]
dn <- cm.invitro[cm.invitro$ENSEMBL %in% share.dn,]


#https://github.com/stephenturner/msigdf/blob/master/data-raw/msigdf.R
#used this code for fetal gene program Ch4 22.05
output.shared <- read.csv("../../../../../../../CM-both-updnns.csv")

DB <- paste("org", "Hs", "eg", "db", sep = ".")
require(DB, character.only = TRUE)
GO2ALLEGS <- paste("org", "Hs", "egGO2ALLEGS", sep = ".")
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]

topgo <- goana(de=dn$ENTREZID,universe=output.shared$ENTREZID,species="Hs")
topgo$adjP <- p.adjust(topgo$P.DE, method="BH")
topgo <- topgo %>% filter(topgo$P.DE < 0.05 & topgo$Ont=="BP") 
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))
de.by.go <- lapply(de.by.go, FUN=function(x) { x[x %in% dn$ENTREZID] })
result <- data.frame(matrix(NA,nrow = dim(topgo)[1],ncol = 5))
for (i in 1:dim(topgo)[1]){
  genename <- dn[dn$ENTREZID %in% de.by.go[[paste(rownames(topgo)[i],topgo$Ont[i],sep = ".")]],]
  result[i,1] <- topgo$Term[i]
  result[i,2] <- topgo$Ont[i]
  result[i,3] <- topgo$P.DE[i]
  result[i,4] <- topgo$adjP[i]
  result[i,5] <- topgo$DE[i]
  result[i,6] <- topgo$N[i]
  result[i,7] <- topgo$DE[i]*100/topgo$N[i]
  vec <- c(genename$SYMBOL)
  fvec <- shQuote(vec, type = "cmd")
  comma_vec <- paste(fvec, collapse = ", ")
  result[i,8] <- comma_vec
}
colnames(result) <- c("Term","Ont","P.DE","adjP","noDE","N","hits","DE")
rownames(result) <- rownames(topgo)
