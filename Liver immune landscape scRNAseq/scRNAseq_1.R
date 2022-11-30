

setwd("C:/Users/LKC/Desktop/Damien/scNASH/Datasets/")

#Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(clustree)
library(tibble)
library(DESeq2)
library(enrichR)

#### Loading data #### 
mat_wk8_ctrl <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk8_ctrl_filtered_feature_bc_matrix/")
mat_wk8_lidpad <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk8_LIDPAD_filtered_feature_bc_matrix/")
mat_wk8_lysm <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk8_LysM_filtered_feature_bc_matrix/")
mat_wk8_ab <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk8_Ab_filtered_feature_bc_matrix/")

mat_wk12_ctrl <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk12_ctrl_filtered_feature_bc_matrix/")
mat_wk12_lidpad <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk12_LIDPAD_filtered_feature_bc_matrix/")
mat_wk12_lysm <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk12_LysM_filtered_feature_bc_matrix/")
mat_wk12_ab <- Read10X(data.dir= "/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/Dataset_1/wk12_Ab_filtered_feature_bc_matrix/")


w8_ctrl<-CreateSeuratObject(counts=mat_wk8_ctrl,project="w8_ctrl",min.cells = 3,min.features = 200)
w8_lidpad<-CreateSeuratObject(counts=mat_wk8_lidpad,project="w8_lidpad",min.cells = 3,min.features = 200)
w8_lysm<-CreateSeuratObject(counts=mat_wk8_lysm,project="w8_lysm",min.cells = 3,min.features = 200)
w8_ab<-CreateSeuratObject(counts=mat_wk8_ab,project="w8_ab",min.cells = 3,min.features = 200)
w12_ctrl<-CreateSeuratObject(counts=mat_wk12_ctrl,project="w12_ctrl",min.cells = 3,min.features = 200)
w12_lidpad<-CreateSeuratObject(counts=mat_wk12_lidpad,project="w12_lidpad",min.cells = 3,min.features = 200)
w12_lysm<-CreateSeuratObject(counts=mat_wk12_lysm,project="w12_lysm",min.cells = 3,min.features = 200)
w12_ab<-CreateSeuratObject(counts=mat_wk12_ab,project="w12_ab",min.cells = 3,min.features = 200)

immune <- merge(w8_ctrl, y = c(w8_lidpad, w8_ab, w8_lysm,
                               w12_ctrl, w12_lidpad, w12_ab, w12_lysm), 
                add.cell.ids = c("w8_ctrl", "w8_lidpad", "w8_Ab","w8_LysM","w12_ctrl", "w12_lidpad", "w12_Ab","w12_LysM"),
                project = "HepaticImmune")

immune$orig.ident<-factor(immune$orig.ident,levels=c("w8_ctrl", "w8_lidpad", "w8_ab","w8_lysm","w12_ctrl", "w12_lidpad", "w12_ab","w12_lysm"))

immune.list<-SplitObject(immune,split.by = "orig.ident")




####Data filtering and QC####

# Add number of genes per UMI for each cell to metadata
DefaultAssay(immune) <- "RNA"
immune$log10GenesPerUMI <- log10(immune$nFeature_RNA) / log10(immune$nCount_RNA)

# Compute percent mito ratio
immune$mitoRatio <- PercentageFeatureSet(object = immune, pattern = "^MT-")
immune$mitoRatio <- immune@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- immune@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=seq_folder, x=nUMI, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=seq_folder, y=log10(nGene), fill=seq_folder)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~seq_folder)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=seq_folder, x=mitoRatio, fill=seq_folder)) + 
  geom_density(alpha = 0.2) + 
  #scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1) + 
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = seq_folder, fill=seq_folder)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# housekeeping genes
# Load the the list of house keeping genes
hkgenes <- read.table(paste0("/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/Datasets/tirosh_housekeepers.txt"), skip = 2)
hkgenes <- as.vector(hkgenes$V1)
# remove hkgenes that were not found
hkgenes.found <- which(toupper(rownames(immune@assays$RNA@data)) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(immune@assays$RNA@data[hkgenes.found, ] > 0)
immune <- AddMetaData(object = immune, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

#VlnPlot(object = immune, features.plot = c("nGene", "nUMI", "percent.mito","n.exp.hkgenes"), nCol = 4,  point.size.use = 0.1)
VlnPlot(object = immune, features ="n.exp.hkgenes") +
  geom_hline(yintercept = 55)



#### Filter #
f_immune <- subset(x = immune, 
                   subset= (nCount_RNA >= 500) & 
                     (nFeature_RNA >= 250) & 
                     (log10GenesPerUMI > 0.80) & 
                     (mitoRatio < 0.20)  &
                     (n.exp.hkgenes > 55) )

f_immune <- subset(x= immune,
                   cells = intersect(colnames(f_immune),colnames(immune)))


## Gene-level filtering ## 

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = f_immune, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
temp <- CreateSeuratObject(filtered_counts, meta.data = f_immune@meta.data)


f_immune<- subset(f_immune, features=c(paste0("prot_",rownames(immune@assays$prot)),paste0("rna_",intersect(rownames(temp), rownames(f_immune)))) )


# Create f_metadata dataframe
f_metadata <- f_immune@meta.data

# Add cell IDs to f_metadata
f_metadata$cells <- rownames(f_metadata)

# Rename columns
f_metadata <- f_metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Visualize the number UMIs/transcripts per cell
f_metadata %>% 
  ggplot(aes(color=seq_folder, x=nUMI, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
f_metadata %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250)

# Visualize the distribution of genes detected per cell via boxplot
f_metadata %>% 
  ggplot(aes(x=seq_folder, y=log10(nGene), fill=seq_folder)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
f_metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +  
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~seq_folder)

# Visualize the distribution of mitochondrial gene expression detected per cell
f_metadata %>% 
  ggplot(aes(color=seq_folder, x=mitoRatio, fill=seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
f_metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = seq_folder, fill=seq_folder)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)



VlnPlot(object = immune, features ="n.exp.hkgenes") +
  geom_hline(yintercept = 55)


immune<- f_immune
rm(f_immune)


### Cell cycle score 


# Load cell cycle markers
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Check Cell cycle phases
immune_phase <- NormalizeData(immune)

# Score cells for cell cycle
immune_phase <- CellCycleScoring(immune_phase, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# View cell cycle scores and phases assigned to cells                                 
View(immune_phase@meta.data)    

# Identify the most variable genes
immune_phase <- FindVariableFeatures(immune_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
immune_phase <- ScaleData(immune_phase)

# Perform PCA
immune_phase <- RunPCA(immune_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(immune_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
DimPlot(immune_phase,
        reduction = "pca",
        group.by= "Phase")



#### Integration, PCA, clustering ####

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_immune, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_immune <- PrepSCTIntegration(object.list = split_immune, anchor.features = integ_features)



# Find best matches
integ_anchors <- FindIntegrationAnchors(object.list = split_immune,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)

# Integrate across conditions
immune_integ <- IntegrateData(anchorset = integ_anchors, 
                              normalization.method = "SCT")
View(immune_integ)

immune_integ <- RunPCA(immune_integ, verbose = FALSE)


immune_integ$orig.ident<-factor(immune_integ$orig.ident,levels=c("w8_ctrl", "w8_lidpad", "w8_ab","w8_lysm","w12_ctrl", "w12_lidpad", "w12_ab","w12_lysm"))

PCAPlot(immune_integ,
        split.by = "orig.ident")  


#determine dimensionality 
#split_immune <- JackStraw(split_immune,num.replicate = 200)
#split_immune<- ScoreJackStraw(split_immune, dims = 1:50)
#JackStrawPlot(ds23.r.integ, dims = 1:50)
ElbowPlot(immune_integ, ndims=50)



# UMAP and Clustering
set.seed(44)
immune_integ <- RunUMAP(immune_integ, reduction = "pca", dims = 1:35)
immune_integ <- RunTSNE(immune_integ, reduction = "pca", dims = 1:35)

DimPlot(immune_integ, reduction = "tsne")
DimPlot(immune_integ, reduction = "umap")


#immune_integ <- FindNeighbors(immune_integ, reduction = "pca", dims = 1:35)
#immune_integ <- FindClusters(immune_integ, resolution = 0.5)




setwd("/acrc/jinmiao/CJM_lab/Damien/Projects/scMultiOmics_batchEffect/Script/scMultiOmics_project/archive/integ_cluster_res")
immune_integ<-FindNeighbors(immune_integ, dims = 1:35)
reses<-seq(from = 2.0, to = 4.0, by = 0.1)
for (res in reses){
  immune_integ<-FindClusters(immune_integ, resolution = res)
  nores<-gsub(pattern = ".", replacement = "", res, fixed = T)
  png(paste0("integ_umap_res", nores, ".png"), width = 5.5*1.4, 
      height = 5.5, units = "in", res = 300)
  print(DimPlot(immune_integ, reduction = "umap", label = T)+ NoLegend())
  dev.off()
  png(paste0("integ_tsne_res", nores, ".png"), width = 5.5*1.4, 
      height = 5.5, units = "in", res = 300)
  print(DimPlot(immune_integ, reduction = "tsne", label = T)+ NoLegend())
  dev.off()
}

png("integ_clustree.png", width = 11*2.6, height = 11*2.0, 
    units = "in", res = 300)
print(clustree(immune_integ, prefix = "integrated_snn_res."))
dev.off()



# Visualization

# Chose res 3.4
Idents(immune_integ) <- immune_integ$integrated_snn_res.3.4

p1 <- DimPlot(immune_integ, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(immune_integ, reduction = "umap", label = TRUE) + NoLegend()
library(cowplot)
plot_grid(p1, p2)
DimPlot(immune_integ, reduction = "umap", split.by = "orig.ident", label = F) +NoLegend()
DimPlot(immune_integ, reduction = "umap", split.by = "orig.ident", label = T) +NoLegend()


saveRDS(immune_integ, file = "clust_immune_integ.rds")


#### Cell Annotation ####


clust_immune_integ <- readRDS("C:/Users/LKC/Desktop/Damien/scNASH/Datasets/clust_immune_integ.rds")


library(celldex)
library(SingleR)
library(pheatmap)
library(multtest)
library(metap)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

ref_immgen <- celldex::ImmGenData()

{
  immune_sce<- clust_immune_integ
  DefaultAssay(immune_sce)<-"integrated"
  immune_sce[['SCT']]<-NULL
  immune_sce[['RNA']]<-NULL
  immune_sce<-as.SingleCellExperiment(immune_sce)
  
  pred_immgen <- SingleR(test= immune_sce, ref = ref_immgen, labels = ref_immgen$label.main)
  
  #annotation
  table(pred_immgen$labels)
  
  plotScoreHeatmap(pred_immgen)
  
  #annotation VS clusters 
  tab <- table(Assigned=pred_immgen$pruned.labels, Cluster=immune_sce$integrated_snn_res.3.4)
  pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
}

### Using rna 


immune_sce.rna <- clust_immune_integ
DefaultAssay(immune_sce.rna)<-"RNA"
immune_sce.rna[['SCT']]<-NULL
immune_sce.rna[['integrated']]<-NULL
immune_sce.rna<-as.SingleCellExperiment(immune_sce.rna)

set.seed(888)
pred_immgen.rna <- SingleR(test= immune_sce.rna, ref = ref_immgen, labels = ref_immgen$label.main)
pred_immgen.rna.fine <- SingleR(test= immune_sce.rna, ref = ref_immgen, labels = ref_immgen$label.fine)
# annotations 
table(pred_immgen.rna$labels)

plotScoreHeatmap(pred_immgen.rna)

#annotation VS clusters 
tab <- table(Assigned=pred_immgen.rna$pruned.labels, Cluster=immune_sce.rna$integrated_snn_res.3.4)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


### rna vs integrated 
tab <- table(Integrated=pred_immgen$pruned.labels, 
             RNA=pred_immgen.rna$pruned.labels)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))


#Add singleR-ImmGen annotation 

clust_immune_integ <- AddMetaData(clust_immune_integ, pred_immgen.rna$labels, col.name = "SingleR_RNA_ImmGen")
clust_immune_integ <- AddMetaData(clust_immune_integ, pred_immgen.rna.fine$labels, col.name = "SingleR_RNA_ImmGen_fine")

plot1 <- DimPlot(clust_immune_integ, reduction = "umap", group.by = "integrated_snn_res.3.4", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend() + ggtitle('Leiden Clusters')
plot2 <- DimPlot(clust_immune_integ, reduction = "umap", group.by = 'SingleR_RNA_ImmGen', label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend() + ggtitle('SingleR Annotations')
plot1 + plot2



features <- c("Cd3d", "Cd3e", "Cd28", "Trac", "Sell", "Gimap5", "Crem", "Cd4", "Cd8a", "Cd8b1",
              "Ms4a1", "Cd79a", "Ccr7", "Ighm", "Ighd",  "H2-Aa","Cd27", "Cxcr4", "Fcrl4","Cxcr3","Fcer1g", "Zeb2", "Itgax", "Cd86", "Mzb1", "Cd38", "CD138",  "Xbp1", "Iglc2",
              "Fcgr1", "Ly6c2", "Fcgr3a",
              "Cd68", "Mnda", "Marco", "Adgre", "Clec4f", "Timd4",
              "Irf8", "Lilra", 
              "Nkg7", "Gnly",
              "Mki67", "Stmn1", "Pclaf"
)
table(Idents(clust_immune_integ))
DotPlot(clust_immune_integ, features = features, cluster.idents= TRUE)

#FindAllMarkers(clust_immune_integ)

test<- FindMarkers(clust_immune_integ, ident.1 = 0, ident.2 = 22)

#B cell heterogeineity trials
#- Find Plasma cell marker Ighg
test1<- FindMarkers(clust_immune_integ, ident.1 = 10, ident.2 = 8)
test2<- FindMarkers(clust_immune_integ, ident.1 = 10, ident.2 = 15)
test3<- FindMarkers(clust_immune_integ, ident.1 = 10, ident.2 = 50)
View(test)

test1[grep("Ighg*",rownames(test1)),]
DotPlot(clust_immune_integ, features = rownames(test1[grep("Ighg*",rownames(test1)),]), cluster.idents= TRUE)
test1[which(rownames(test1) == c("Zeb2", "Ighv1-81", "Igha")),]

rownames(clust_immune_integ[grep("Ighg*",rownames(clust_immune_integ)),])


rownames(clust_immune_integ[grep("Il4r*",rownames(clust_immune_integ)),])
rownames(clust_immune_integ[grep("Tcl1*",rownames(clust_immune_integ)),])

#checking histocompat transcripts
rownames(clust_immune_integ[grep("H2-a*",rownames(clust_immune_integ)),])
FeaturePlot(clust_immune_integ, features = rownames(clust_immune_integ[grep("H2-a*",rownames(clust_immune_integ)),]))

#checking Ig transcripts 
rownames(clust_immune_integ[grep("Igh*",rownames(clust_immune_integ)),]) #Ighm, Ighd, Igha, Ighvs...
FeaturePlot(clust_immune_integ, features = c("Igkc", "Iglc1", "Iglc2", "Iglc3","Iglv1", "Ighm", "Ighd", "Igha", "H2-Aa"))

VlnPlot(clust_immune_integ, features = c("Igkc", "Iglc1", "Iglc2", "Iglc3","Iglv1", "Ighm", "Ighd", "Igha", "H2-Aa"))

FeaturePlot(clust_immune_integ, features = c("Ms4a1", "Cd79a", "Cd38", "Cd5", "Cd80", "Cd86", "Xbp1","Mzb1", "Ighm", "Ighd", "Igha", "H2-Aa"))


saveRDS(clust_immune_integ, file = "C:/Users/LKC/Desktop/Damien/scNASH/Datasets/clust_immune_integ.rds")



#05/04/22 - trying to label 43,44,48,58 
#use different ref - DICE 
#http://bioconductor.org/books/release/SingleRBook/advanced-options.html


library(celldex)
dice <- DatabaseImmuneCellExpressionData(ensembl=FALSE)

immune_sce.rna <- clust_immune_integ
DefaultAssay(immune_sce.rna)<-"RNA"
immune_sce.rna[['SCT']]<-NULL
immune_sce.rna[['integrated']]<-NULL
immune_sce.rna<-as.SingleCellExperiment(immune_sce.rna)
set.seed(888)


common <- intersect(rownames(immune_sce.rna), toupper(rownames(dice)))
library(SingleR)
set.seed(888)
trained <- trainSingleR(dice[common,], labels=dice$label.fine, aggr.ref=TRUE)
pred_dice.rna <- classifySingleR(immune_sce.rna, trained, assay.type=1)
table(pred$labels)

clust_immune_integ <- AddMetaData(clust_immune_integ, pred_immgen.rna$labels, col.name = "SingleR_RNA_DICE")
#clust_immune_integ <- AddMetaData(clust_immune_integ, pred_immgen.rna.fine$labels, col.name = "SingleR_RNA_DICE_fine")

tab <- table(Assigned=pred_dice.rna$pruned.labels, Cluster=immune_sce.rna$integrated_snn_res.3.4)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))

tab1 <- table(Assigned=pred_dice.rna$pruned.labels, Cluster=immune_sce.rna$SingleR_RNA_ImmGen)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))



# + Manual curation

#overall 

Idents(clust_immune_integ) <- clust_immune_integ$integrated_snn_res.3.4
clust_immune_integ <- RenameIdents(clust_immune_integ, 
                                   
                                   '0' = "Mature B cells",
                                   '1' = "Mature B cells",
                                   '2' = "Mature B cells",
                                   '3' = "Mature B cells",
                                   '4' = "Mature B cells",
                                   '6' = "Mature B cells",
                                   '8' = "Mature B cells",
                                   '9' = "Mature B cells",
                                   '19' = "Mature B cells",
                                   '20' = "Mature B cells",
                                   '32' = "Mature B cells",
                                   '37' = "Mature B cells",
                                   '15' = "TrB cells",
                                   '50' = "TrB cells",
                                   '10' = "B_c10",
                                   '29' = "B_c29",
                                   '52' = "B_c52",
                                   
                                   
                                   "5" = "Naive CD8 T" ,
                                   "7" = "Naive CD8 T",
                                   "12" = "Naive CD8 T",
                                   "13" = "_13 Foxp3mid Treg",
                                   "14" = "Naive CD8 T",
                                   "16" = "Naive CD4 T",
                                   "17" = "_17",
                                   "18" = "Naive CD4 T",
                                   "21" = "_21 CD4- CD8- T",
                                   "22" = "_22",
                                   "25" = "gdT",
                                   "26" = "CCR6+ Tgd",
                                   "27" = "_27Foxp3+ Treg?",
                                   "28" = "Naive CD4 T",
                                   "30" = "_30 doublen neg NKG7+",
                                   "35" = "_35",
                                   "38" = "Th1",
                                   "39" = "Naive CD4 T",
                                   "40" = "_40 doubleneg",
                                   "41" = "Th17",
                                   "43" = "_43 nothing.",
                                   "46" = "Prolif. CD8 T",
                                   "48" = "_48 same as 43. nth",
                                   "49" = "_49 T-innate like ",
                                   "51" = "_51 CD8 mid",
                                   "55" = "_55 NKT",
                                   "56" = "_56 Ccl4+ Cd8a+ T ",
                                   "60" = "_60 T-innate Ccl4+",
                                   
                                   "11" = "NK cells" ,
                                   "23" = "NK cells",
                                   "24" = "NK cells",
                                   "33" = "NK cells",
                                   "45" = "NK cells",
                                   "42" = "NK cells",
                                   "34" = "ILC_c34",
                                   "36" = "ILC_c36",
                                   "44" = "_Unk_profl_c44",
                                   "47" = "Neutrophils",
                                   "31" = "Monocytes",
                                   "57" = "Kupffer cells",
                                   "58" = "_Unk_c58",
                                   "59" = "cDC",
                                   "61" = "pDC",
                                   "53" = "Mo-macs",
                                   "54" = "plasma cells")  


clust_immune_integ <- RenameIdents(clust_immune_integ,
                                   "5" = "Naive CD8 T" ,
                                   "7" = "Naive CD8 T",
                                   "12" = "Naive CD8 T",
                                   "_13 Foxp3mid Treg" = "Treg_c13",
                                   "14" = "Naive CD8 T",
                                   "16" = "Naive CD4 T",
                                   "_17" = "NKT_c17",
                                   "18" = "Naive CD4 T",
                                   "_21 CD4- CD8- T" = "NKT_c21",
                                   "_22" = "Memory CD8+ Tem",
                                   "25" = "gdT",
                                   "CCR6+ Tgd" = "CCR6+ gdT",
                                   "_27Foxp3+ Treg?" = "Treg_c27",
                                   "28" = "Naive CD4 T",
                                   "_30 doublen neg NKG7+" = "NKT_c30",
                                   "_35" = "Activated dnT_35",
                                   "38" = "Tnf+ CD4 CTL",
                                   "39" = "Naive CD4 T",
                                   "_40 doubleneg" = "Naive dnT_c40",
                                   "41" = "Th17",
                                   "_43 nothing." = "_43 Bmarkers",
                                   "Prolif. CD8 T" = "Memory CD8 T_c46",
                                   "_48 same as 43. nth" = "_48 Unk",
                                   "_49 T-innate like " = "T-innate-like_c49",
                                   "_51 CD8 mid" = "Activated dnT_c51",
                                   "_55 NKT" = "Activated NKT_c55",
                                   "_56 Ccl4+ Cd8a+ T " = "Memory CD8+ Tem",
                                   "_60 T-innate Ccl4+" = "Memory CD8+ Trm"
)  




DimPlot(clust_immune_integ, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) +  ggtitle('Annotations')

clust_immune_integ$annotv1 <- Idents(clust_immune_integ) 



all.markers <- FindAllMarkers(clust_immune_integ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers1<- all.markers
all.markers1 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
all.markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(subset(clust_immune_integ, downsample = 50), features = top20$gene) + NoLegend()

write.csv(all.markers, "C:/Users/LKC/Desktop/Damien/scNASH/Datasets/all_annotv1_DEGs.csv")

names(table(clust_immune_integ$annotv1))

Idents(clust_immune_integ) <- clust_immune_integ$annotv1
clust_immune_integ <- RenameIdents(clust_immune_integ,
                                   "Treg_c13" = "Treg",
                                   "NKT_c17" = "NKT",          
                                   "NKT_c21" = "NKT",            
                                   "Memory CD8+ Tem" = "CD8 T cells",   
                                   "CCR6+ gdT" = "gdT",         
                                   "Treg_c27" =  "Treg",         
                                   "NKT_c30" = "NKT",
                                   "Activated dnT_35" = "dnT",  
                                   "Naive dnT_c40" = "Naive dnT",     
                                   "_43 Bmarkers" = "43",      
                                   "Memory CD8 T_c46"  = "CD8 T cells",
                                   "_48 Unk"  = "48",
                                   "T-innate-like_c49" = "T-innate-like",
                                   "Activated dnT_c51" = "dnT",
                                   "Activated NKT_c55" = "NKT",
                                   "Memory CD8+ Trm" = "CD8 T cells",   
                                   "Mature B cells" = "B cells",    
                                   "TrB cells" = "B cells",         
                                   "B_c10" = "B cells",             
                                   "B_c29" = "B cells",            
                                   "B_c52" = "B cells",
                                   "Naive CD8 T" = "Naive CD8 T",      
                                   "Naive CD4 T" = "Naive CD4 T",      
                                   "gdT" = "gdT",              
                                   "Th1"  = "Th1",    
                                   "Th17" = "Th17", 
                                   
                                   "NK cells" = "NK cells",         
                                   "ILC_c34"= "ILC",
                                   "ILC_c36" = "ILC",          
                                   "_Unk_profl_c44" = "44",   
                                   "Neutrophils"  = "Neutrophils",     
                                   "Monocytes"  = "Monocytes",       
                                   "Kupffer cells" = "Kupffer cells",    
                                   "_Unk_c58"   = "58",       
                                   "cDC"= "cDC",
                                   "pDC"= "pDC",              
                                   "Mo-macs" = "Mo-macs",          
                                   "plasma cells" = "plasma cells"
                                   
)  

clust_immune_integ$annotv1_maj1 <- Idents(clust_immune_integ)


Idents(clust_immune_integ) <- clust_immune_integ$annotv1
clust_immune_integ <- RenameIdents(clust_immune_integ,
                                   "Treg_c13" = "CD4 T cells",
                                   "NKT_c17" = "NKT",          
                                   "NKT_c21" = "NKT",            
                                   "Memory CD8+ Tem" = "CD8 T cells",   
                                   "CCR6+ gdT" = "gdT",         
                                   "Treg_c27" =  "CD4 T cells",         
                                   "NKT_c30" = "NKT",
                                   "Activated dnT_35" = "dnT",  
                                   "Naive dnT_c40" = "dnT",     
                                   "_43 Bmarkers" = "43",      
                                   "Memory CD8 T_c46"  = "CD8 T cells",
                                   "_48 Unk"  = "48",
                                   "T-innate-like_c49" = "T-innate-like",
                                   "Activated dnT_c51" = "dnT",
                                   "Activated NKT_c55" = "NKT",
                                   "Memory CD8+ Trm" = "CD8 T cells",   
                                   "Mature B cells" = "B cells",    
                                   "TrB cells" = "B cells",         
                                   "B_c10" = "B cells",             
                                   "B_c29" = "B cells",            
                                   "B_c52" = "B cells",
                                   "Naive CD8 T" = "Naive CD8 T",      
                                   "Naive CD4 T" = "Naive CD4 T",      
                                   "gdT" = "gdT",              
                                   "Th1"  = "CD4 T cells",    
                                   "Th17" = "CD4 T cells", 
                                   
                                   "NK cells" = "NK cells",         
                                   "ILC_c34"= "ILC",
                                   "ILC_c36" = "ILC",          
                                   "_Unk_profl_c44" = "44",   
                                   "Neutrophils"  = "Neutrophils",     
                                   "Monocytes"  = "Monocytes",       
                                   "Kupffer cells" = "Kupffer cells",    
                                   "_Unk_c58"   = "58",       
                                   "cDC"= "cDC",
                                   "pDC"= "pDC",              
                                   "Mo-macs" = "Mo-macs",          
                                   "plasma cells" = "plasma cells"
                                   
)  


clust_immune_integ$annotv1_maj2 <- Idents(clust_immune_integ)

mycol <- DiscretePalette(n=25, palette ="polychrome")
names(mycol)= c( "Naive CD8 T",
                 "CD8 T cells",
                 "Naive CD4 T",
                 "CD4 T cells",
                 "Naive dnT",
                 "dnT",
                 "T-innate-like",
                 "gdT",
                 "NKT",
                 "B cells",
                 "plasma cells",
                 "ILC",
                 "NK cells",
                 "Neutrophils",
                 "Monocytes",
                 "Kupffer cells",
                 "Mo-macs",
                 "cDC",
                 "pDC",
                 "43",
                 "44",
                 "48",
                 "58")



#### Cell abundances ####

Idents(clust_immune_integ) <- clust_immune_integ$annotv1_maj2
clust_immune_integ <- RenameIdents(
  clust_immune_integ, "Naive CD8 T" = "T lymphocytes",
  "CD8 T cells"= "T lymphocytes",
  "Naive CD4 T"= "T lymphocytes",
  "CD4 T cells"= "T lymphocytes",
  "Naive dnT"= "T lymphocytes",
  "dnT"= "T lymphocytes",
  "T-innate-like"= "T lymphocytes",
  "gdT"= "T lymphocytes",
  "NKT"= "T lymphocytes",
  "B cells"= "B lymphocytes",
  "plasma cells"= "B lymphocytes",
  "ILC" = "myeloid-derived",
  "NK cells" = "myeloid-derived",
  "Neutrophils"= "myeloid-derived",
  "Monocytes"= "myeloid-derived",
  "Kupffer cells"= "myeloid-derived",
  "Mo-macs"= "myeloid-derived",
  "cDC"= "myeloid-derived",
  "pDC"= "myeloid-derived") 


new <- subset(clust_immune_integ, idents= c("T lymphocytes", "B lymphocytes", "myeloid-derived") )

new$major <- Idents(new)



Idents(new) <- new$orig.ident
split.all.clust<- SplitObject(new , split.by = "orig.ident")
split.all.wk8 <- split.all.clust[1:4]
split.all.wk12 <- split.all.clust[5:8]


##### Wk 8 all #####
all.prop.wk8 <- c()
for (i in 1:length(split.all.wk8) ){
  all.prop.wk8<- cbind(all.prop.wk8, as.vector(summary(split.all.wk8[[i]]$annotv1_maj2)))
}

all.prop.wk8 <- as.data.frame(all.prop.wk8)
rownames(all.prop.wk8) <- names(summary(new$annotv1_maj2))
colnames(all.prop.wk8) <- names(summary(new$orig.ident))[1:4]

all.prop.wk8[all.prop.wk8 == 0] <-NA
all.prop.wk8 <- na.omit(all.prop.wk8)

immune.total <- as.data.frame(summary(new$orig.ident))

all.prop.wk8.pct<-c()
for (i in 1:ncol(all.prop.wk8)){
  all.prop.wk8.pct <- cbind(all.prop.wk8.pct, all.prop.wk8[,i] / immune.total[,1][i])
}

all.prop.wk8.pct <- as.data.frame(all.prop.wk8.pct)
rownames(all.prop.wk8.pct) <- rownames(all.prop.wk8)
colnames(all.prop.wk8.pct) <- names(summary(new$orig.ident))[1:4]

all.prop.wk8.pct$Celltype <- rownames(all.prop.wk8.pct)


all.prop.wk8.pct.log2fc <- log2(all.prop.wk8.pct[1:4]/all.prop.wk8.pct$w8_ctrl)



all.prop.pct <- cbind(all.prop.wk8.pct[,1:4], all.prop.wk12.pct)

all.prop.pct$Celltype <- rownames(all.prop.pct)
all.prop.pct <- all.prop.pct[!(rownames(all.prop.pct) %in% c(43,44,48,58)), ]

all.prop.pct.gg<-all.prop.pct%>%tidyr::gather(Treatment,Count,1:8)
all.prop.pct.gg$Treatment<-factor(all.prop.pct.gg$Treatment, 
                                  levels= c('w8_ctrl','w8_lidpad','w8_ab','w8_lysm',
                                            'w12_ctrl','w12_lidpad','w12_ab','w12_lysm'))


all.prop.pct.gg$Celltype<-factor(all.prop.pct.gg$Celltype, 
                                 levels= c("T lymphocytes", "B lymphocytes", "myeloid-derived") )
all.prop.pct.gg$Celltype <- forcats::fct_rev(all.prop.pct.gg$Celltype)
all.prop.pct.gg$Count <- round(all.prop.pct.gg$Count, digits = 5)


all.prop.pct.gg <- all.prop.pct.gg[which(!is.na(all.prop.pct.gg$Celltype)),]
all.prop.w8.pct.ggplot<-ggplot(all.prop.pct.gg[1:72,],aes(fill=Treatment,y=Count,x=Celltype))+
  geom_bar(position="stack", stat = 'identity')+
  guides(fill=guide_legend(reverse=F))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Celltype")+
  ylab("Relative proportion") 



##### WK 12 all#####

all.prop.wk12 <- c()
for (i in 1:length(split.all.wk12) ){
  all.prop.wk12<- cbind(all.prop.wk12, as.vector(summary(split.all.wk12[[i]]$annotv1_maj2)))
}

all.prop.wk12 <- as.data.frame(all.prop.wk12)
rownames(all.prop.wk12) <- names(summary(new$annotv1_maj2))
colnames(all.prop.wk12) <- names(summary(new$orig.ident))[5:8]
all.prop.wk12[all.prop.wk12 == 0] <-NA
all.prop.wk12 <- na.omit(all.prop.wk12)

immune.total <- as.data.frame(summary(new$orig.ident))

all.prop.wk12.pct<-c()
for (i in 1:ncol(all.prop.wk12)){
  all.prop.wk12.pct <- cbind(all.prop.wk12.pct, all.prop.wk12[,i] / immune.total[,1][i+4])
}
all.prop.wk12.pct <- as.data.frame(all.prop.wk12.pct)
rownames(all.prop.wk12.pct) <- rownames(all.prop.wk12)
colnames(all.prop.wk12.pct) <- names(summary(new$orig.ident))[5:8]

#sanity check 
colSums(all.prop.wk12.pct) #good.

all.prop.wk12.pct$Celltype <- rownames(all.prop.wk12.pct)
all.prop.wk12.pct.log2fc <- log2(all.prop.wk12.pct[1:4]/all.prop.wk12.pct$w12_ctrl)







#### Figure 4 - global prop #### 
all.prop.w8.pct.ggplot + ggtitle("Week 8") + coord_flip() 
all.prop.w8.pct.ggplot + coord_flip(ylim = c(0,0.01), xlim = c(0.95,3.9)) + theme_bw() #last 4
all.prop.w8.pct.ggplot + coord_flip(ylim = c(0,0.015), xlim = c(9.1,9.1)) + theme_bw() # plasma cells
write.csv(all.prop.pct.gg, "C:/Users/LKC/Desktop/Damien/scNASH/Output/Fig4_global_proportion.csv")